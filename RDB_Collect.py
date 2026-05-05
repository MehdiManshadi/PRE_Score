from Mixer import Region


client = Region(token="f46027d576ea")
rsid = "rs6589939"
RegDB_Coord, _ = client.get_ldproxy(rsid)

data = client.query_regulomedb(RegDB_Coord[0])

for Ex in data["@graph"]:
    if Ex['biosample_ontology']['term_name'] == "B cell" and (Ex['target_label'] is not None):
        print(Ex["target_label"])

for record in data['@graph']:
    if record['method'] == 'chromatin state':
        cell_type = record['biosample_ontology']['term_name']
        if cell_type == "B cell":
            chrom_state = record['value']
            print(cell_type, chrom_state)