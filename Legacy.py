import os
import requests

LDLINK_BASE_URL = "https://ldlink.nih.gov/LDlinkRest/ldproxy"


def get_ldproxy_rs(
    rsid,
    population="CEU",
    window=50000,
    collapse_transcripts=True,
    annotation="RegulomeDB",
    ld_measure="R2",
    token="f46027d576ea",
):
    if token is None:
        token = os.getenv("LDLINK_TOKEN")
    if not token:
        raise ValueError(
            "LDlink API token required. Set LDLINK_TOKEN env var or pass token argument."
        )

    params = {
        "var": rsid,
        "pop": population,
        "r2_d": ld_measure,
        "window": str(window),
        "collapsed": "yes" if collapse_transcripts else "no",
        "annot": annotation,
        "token": token,
    }

    response = requests.get(LDLINK_BASE_URL, params=params, timeout=30)
    response.raise_for_status()
    return response.text


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Query LDlink LDproxy with a variant rsID.")
    parser.add_argument("rsid", help="Variant rsID")
    parser.add_argument("--token", help="LDlink API token")
    parser.add_argument("--population", default="CEU")
    parser.add_argument("--window", type=int, default=50000, help="Window size in base pairs")
    parser.add_argument("--ld-measure", default="R2", help="LD measure, e.g. R2 or D")
    parser.add_argument("--no-collapse", dest="collapse_transcripts", action="store_false", help="Do not collapse transcripts")
    parser.add_argument("--annotation", default="RegulomeDB", help="Annotation source")
    args = parser.parse_args()

    result = get_ldproxy_rs(
        args.rsid,
        population=args.population,
        window=args.window,
        collapse_transcripts=args.collapse_transcripts,
        annotation=args.annotation,
        ld_measure=args.ld_measure,
        token=args.token,
    )
    print(result)


if __name__ == "__main__":
    main()
