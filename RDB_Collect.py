import requests

def query_regulomedb(rsid, assembly="GRCh38"):
    url = "https://regulomedb.org/regulome-summary/"

    params = {
        "regions": rsid,
        "genome": assembly,
        "maf": "0.01"
    }

    headers = {
        "Accept": "application/json",
        "User-Agent": "Mozilla/5.0"
    }

    response = requests.get(url, params=params, headers=headers)
    response.raise_for_status()

    return response.json()

if __name__ == "__main__":
    rsid = "rs6589939"
    assembly = "GRCh38"
    data = query_regulomedb(rsid, assembly)

    for Ex in data["@graph"]:
        if Ex['biosample_ontology']['term_name'] == "B cell" and (Ex['target_label'] is not None):
            print(Ex["target_label"])
    
    for record in data['@graph']:
        if record['method'] == 'chromatin state':
            cell_type = record['biosample_ontology']['term_name']
            if cell_type == "B cell":
                chrom_state = record['value']
                print(cell_type, chrom_state)


from Mixer import LDlinkClient as Lead
if __name__ == "__main__":
    client = Lead(token="f46027d576ea")
    rsid = "rs6589939"
    df = client.get_ldproxy(rsid)
    print(df.head())