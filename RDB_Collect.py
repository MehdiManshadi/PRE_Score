import os
import re
import requests
from bs4 import BeautifulSoup
from flask import Flask, jsonify, request

app = Flask(__name__)

REGULOME_BASE_URL = "http://legacy.regulomedb.org/index"
HEADERS = {
    "User-Agent": "Mozilla/5.0 (compatible; PRE_Score/1.0; +https://example.com)"
}


def fetch_regulome_data(rsid):
    params = {"rsid": rsid}
    resp = requests.get(REGULOME_BASE_URL, params=params, headers=HEADERS, timeout=15)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, "html.parser")

    gene = None
    for label in soup.find_all(string=re.compile(r"\bGene\b", re.I)):
        parent = label.parent
        if parent and parent.name in ("td", "th"):
            sibling = parent.find_next_sibling("td")
            if sibling:
                gene = sibling.get_text(strip=True)
                break

    if not gene:
        dl = soup.find("dl")
        if dl:
            for dt, dd in zip(dl.find_all("dt"), dl.find_all("dd")):
                if dt.get_text(strip=True).lower().startswith("gene"):
                    gene = dd.get_text(strip=True)
                    break

    genome_browser_url = None
    for a in soup.find_all("a", href=True):
        text = a.get_text(separator=" ", strip=True).lower()
        if "genome" in text and "browser" in text:
            genome_browser_url = requests.compat.urljoin(REGULOME_BASE_URL, a["href"])
            break
        if "browser" in text and "genome" in text:
            genome_browser_url = requests.compat.urljoin(REGULOME_BASE_URL, a["href"])
            break

    snp = rsid
    if genome_browser_url:
        m = re.search(r"(chr[0-9XY]+[:]\d+)", genome_browser_url, re.I)
        if m:
            snp = m.group(1)

    return {
        "rsid": rsid,
        "query_url": resp.url,
        "gene": gene,
        "snp": snp,
        "genome_browser_url": genome_browser_url,
    }


@app.route("/api/v1/regulome", methods=["GET"])
def regulome_api():
    rsid = request.args.get("rsid")
    if not rsid:
        return jsonify({"error": "missing rsid query parameter"}), 400

    try:
        data = fetch_regulome_data(rsid)
    except requests.HTTPError as exc:
        return jsonify({"error": "failed to fetch RegulomeDB page", "details": str(exc)}), 502
    except requests.RequestException as exc:
        return jsonify({"error": "request to RegulomeDB failed", "details": str(exc)}), 502
    except Exception as exc:
        return jsonify({"error": "failed to parse RegulomeDB response", "details": str(exc)}), 500

    return jsonify(data)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.getenv("PORT", "5000")), debug=False)
