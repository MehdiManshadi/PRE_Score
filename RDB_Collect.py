from Mixer import Region


client = Region(token="f46027d576ea", assembly="GRCh38")
rsid = ["rs10801908", "rs11256593", "rs438613"]
for r in rsid:
    Coord, _ = client.get_ldproxy(r)

    for c in Coord:
        print("RSID:", r, "Coord:", c, "-------------------------------------")
        client.coord = c  # Set the coordinate for the client instance
        score, n_experiments, ChIP_Seq_targets, Dataset = client.count_experiments(cell_type="B cell")
        #print(f"RegulomeDB Score: {score}, Number of Experiments: {n_experiments}")
        #print(f"ChIP-Seq Targets: {', '.join(ChIP_Seq_targets)}")
        #print(f"Datasets: {', '.join(Dataset)}")

    #print(f"RegulomeDB Score: {score}, Number of Experiments: {n_experiments}")
    #print(f"ChIP-Seq Targets: {', '.join(ChIP_Seq_targets)}")


