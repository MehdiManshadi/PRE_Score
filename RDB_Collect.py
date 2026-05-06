from Mixer import Region


client = Region(token="f46027d576ea")
rsid = "rs8062446"
Coord, _ = client.get_ldproxy(rsid)

for c in Coord:
    print(c, "-------------------------------------")
    data = client.query_regulomedb(c)
    gene = client.query_ensembl(c)
    client.count_experiments(data, cell_type="B cell")


