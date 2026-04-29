import requests
from bs4 import BeautifulSoup

def query_legacy_regulomedb(rsid):
    url = "http://legacy.regulomedb.org/results"

    payload = {
        "query": rsid
    }

    headers = {
        "User-Agent": "Mozilla/5.0"
    }

    response = requests.post(url, data=payload, headers=headers)
    response.raise_for_status()

    return response.text


def extract_bcell_targets(html):
    soup = BeautifulSoup(html, "html.parser")

    results = []

    # Example: iterate over table rows
    rows = soup.find_all("tr")

    for row in rows:
        text = row.get_text(separator=" ").strip()

        if "B cell" in text:
            results.append(text)

    return results


if __name__ == "__main__":
    rsid = "rs8062446"

    html = query_legacy_regulomedb(rsid)
    results = extract_bcell_targets(html)

    for r in results:
        print(r)