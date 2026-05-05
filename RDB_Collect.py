import requests

from Mixer import Region


client = Region(token="f46027d576ea")
rsid = "rs8062446"
Coord, _ = client.get_ldproxy(rsid)

data = client.query_regulomedb(Coord[0])

gene = client.query_ensembl(Coord[0])

for experiments in data["@graph"]:
    if experiments['biosample_ontology']['term_name'] == "B cell" and (experiments['target_label'] is not None):
        print(experiments["target_label"])

for record in data['@graph']:
    if record['method'] == 'chromatin state':
        cell_type = record['biosample_ontology']['term_name']
        if cell_type == "B cell":
            chrom_state = record['value']
            print(cell_type, chrom_state)


