from Mixer import Region


client = Region(token="f46027d576ea")
rsid = "rs8062446"
Coord, _ = client.get_ldproxy(rsid)

for c in Coord:
    print(c, "-------------------------------------")
    data = client.query_regulomedb(c)
    gene = client.query_ensembl(c)
    score, n_experiments, ChIP_Seq_targets = client.count_experiments(data, cell_type="B cell")
    print(f"RegulomeDB Score: {score}, Number of Experiments: {n_experiments}")
    print(f"ChIP-Seq Targets: {', '.join(ChIP_Seq_targets)}")

#print(f"RegulomeDB Score: {score}, Number of Experiments: {n_experiments}")
#print(f"ChIP-Seq Targets: {', '.join(ChIP_Seq_targets)}")


