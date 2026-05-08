import requests
import os
import pandas as pd
from io import StringIO
import requests
from datetime import datetime

class Region:
    BASE_URL = "https://ldlink.nih.gov/LDlinkRest/ldproxy"

    def __init__(self, token=None, assembly="GRCh38"):
        self.token = token or os.getenv("LDLINK_TOKEN")
        self.assembly = assembly
        if not self.token:
            raise ValueError("LDlink API token required")
        
    def ldlink_to_dataframe(self, response_text):
        df = pd.read_csv(StringIO(response_text), sep="\t")
        return df
    
    def get_ldproxy(self, rsid, population="CEU", window=50000,
                    collapse_transcripts=True, annotation="RegulomeDB",
                    ld_measure="R2"):

        self.rsid = rsid

        params = {
            "var": rsid,
            "pop": population,
            "r2_d": ld_measure,
            "window": str(window),
            "collapsed": "yes" if collapse_transcripts else "no",
            "annot": annotation,
            "token": self.token,
            "genome_build": self.assembly.lower()
        }

        response = requests.get(self.BASE_URL, params=params, timeout=30)
        response.raise_for_status()
        
        df = self.ldlink_to_dataframe(response.text)
        Coord = df[df["R2"] >= 0.5]["Coord"]
        RS_Number = df[df["R2"] >= 0.5]["RS_Number"]
        return Coord, RS_Number

    def query_regulomedb(self):
        
        url = "https://regulomedb.org/regulome-summary/"
        
        RegDB_Coord = f"{self.coord.split(':')[0]}:{int(self.coord.split(':')[1])-1}-{self.coord.split(':')[1]}"

        params = {
            "regions": RegDB_Coord,
            "genome": self.assembly,
            "maf": "0.01"
        }

        headers = {
            "Accept": "application/json",
            "User-Agent": "Mozilla/5.0"
        }

        response = requests.get(url, params=params, headers=headers)
        response.raise_for_status()

        return response.json()
    
    def query_ensembl(self, coord):
        chrom, pos = coord.replace("chr", "").split(":")
        if self.assembly == "GRCh37":
            server = "https://grch37.rest.ensembl.org"
        elif self.assembly == "GRCh38": 
            server = "https://rest.ensembl.org"

        ext = f"/overlap/region/human/{chrom}:{pos}-{pos}?feature=gene"
        headers = {"Content-Type": "application/json"}

        response = requests.get(server + ext, headers=headers)
        response.raise_for_status()
        genes = response.json()

        #for gene in genes:
        #   print(gene["id"], gene["external_name"], gene["biotype"])

        return [gene["external_name"] for gene in genes if gene.get("external_name")] if genes else None
    
    def count_experiments(self, cell_type):
        score = 0
        ChIP_Seq_targets = []
        Dataset = []
        n_experiments = 0
        genes = self.query_ensembl(self.coord)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"regulomedb_experiments_{timestamp}.txt"
        
        if not genes:
            return [], [], [], []
        
        data = self.query_regulomedb()

        histone_marks = {
            "H2az": +1,  # associated with regulatory elements
            "H3K4me1": +1,  # enhancer associated
            "H3K04me1": +1,  # enhancer associated
            "H3K4me2": +1,  # promoter/enhancer
            "H3K04me2": +1,  # promoter/enhancer
            "H3K4me3": +1,  # promoters/transcription start
            "H3K04me3": +1,  # promoters/transcription start
            "H3K04me3Ohtam": +1,  # assume promoter/activating
            "H3K9ac": +1,  # active regulatory elements, promoters
            "H3K09ac": +1,  # active regulatory elements, promoters
            "H3K9me1": 0,  # preference for 5' end, neutral
            "H3K09me1": 0,  # preference for 5' end, neutral
            "H3K9me3": -1,  # repressive heterochromatin
            "H3K09me3": -1,  # repressive heterochromatin
            "H3K27ac": +1,  # active enhancers/promoters
            "H3K27me3": -1,  # polycomb repression
            "H3K36me3": 0,  # elongation mark, neutral
            "H3K79me2": 0,  # transcription-associated, neutral
            "H4K20me1": 0   # preference for 5' end, neutral
        }

        chromatin_state_scores = {
            'TssA': 1,           # Active TSS
            'TssFlnk': 1,        # Flanking TSS (15-state model)
            'TssFlnkU': 1,       # Flanking TSS Upstream
            'TssFlnkD': 1,       # Flanking TSS Downstream
            'Tx': 1,             # Strong transcription
            'TxWk': 1,           # Weak transcription
            'EnhG1': 1,          # Genic enhancer 1
            'EnhG2': 1,          # Genic enhancer 2
            'EnhA1': 1,          # Active Enhancer 1
            'EnhA2': 1,          # Active Enhancer 2
            'EnhWk': 1,          # Weak Enhancer
            'ZNF/Rpts': 0,       # ZNF genes & repeats
            'Het': -1,           # Heterochromatin
            'TssBiv': 0,         # Bivalent/Poised TSS
            'EnhBiv': 0,         # Bivalent Enhancer
            'ReprPC': -1,        # Repressed Polycomb
            'ReprPCWk': -1,      # Weak Repressed Polycomb
            'Quies': 0           # Quiescent/Low signal
        }
        with open(filename, "a") as f:
            for experiment in data["@graph"]:

                biosample = experiment.get("biosample_ontology", {})
                term_name = biosample.get("term_name")    

                if term_name == ["ChIP-seq"] or term_name == ["Histone ChIP-seq"] or term_name == ["chromatin state"] or term_name == ["DNase-seq"] or term_name == ["ATAC-seq"]:
                    print("No cell type information for this experiment, skipping...")
                    print("================================================================================-====================")
                    ValueError("No cell type information for this experiment, skipping...")

                if term_name == cell_type:
                    dataset_rel = experiment.get('dataset_rel')
                    dataset_rel = dataset_rel.strip("/").split("/")[-1]
                    Dataset.append(dataset_rel)   
                

                    if experiment.get("method") == "ChIP-seq":
                        target = experiment.get("target_label")
                        ChIP_Seq_targets.append(target)
                        #print(cell_type, experiment.get("method"),experiment["target_label"])
                        for gene in genes:
                            f.write(f"{self.rsid},{gene}, {dataset_rel}, {self.coord}, {cell_type},{experiment.get('method')},{experiment.get('target_label')}, {""}\n")
                    
                    if experiment.get("method") == "Histone ChIP-seq":
                        target = experiment.get("target_label")
                        score += histone_marks[target]
                        n_experiments += 1  # Count the number of histone ChIP-seq experiments for this cell type
                        #print(cell_type, experiment.get("method"),experiment["target_label"])
                        for gene in genes:
                            f.write(f"{self.rsid},{gene}, {dataset_rel}, {self.coord},{cell_type},{experiment.get('method')},{experiment.get('target_label')}, {histone_marks[target]}\n")

                    if experiment.get("method") == "chromatin state":
                        chrom_state = experiment.get("value")
                        score = chromatin_state_scores.get(chrom_state, 0) + score
                        n_experiments += 1  # Count the number of chromatin state experiments for this cell type
                        #print(cell_type, experiment.get("method"), chrom_state)
                        for gene in genes:
                            f.write(f"{self.rsid},{gene}, {dataset_rel}, {self.coord}, {cell_type},{experiment.get('method')}, {chrom_state}, {chromatin_state_scores.get(chrom_state, 0)}\n")

                    if experiment.get("method") == "DNase-seq":
                        value = experiment.get("value")
                        # n_experiments += 1  # Count the number of chromatin state experiments for this cell type
                        #print(cell_type, experiment.get("method"), value)
                        for gene in genes:
                            f.write(f"{self.rsid},{gene}, {dataset_rel}, {self.coord}, {cell_type},{experiment.get('method')}, {value}, {''}\n")

                    if experiment.get("method") == "ATAC-seq":
                        value = experiment.get("value")
                        # n_experiments += 1  # Count the number of chromatin state experiments for this cell type
                        #print(cell_type, "\t", experiment.get("method"), "\t", value)
                        for gene in genes:
                            f.write(f"{self.rsid},{gene}, {dataset_rel}, {self.coord}, {cell_type},{experiment.get('method')}, {value}, {''}\n")


        return score, n_experiments, ChIP_Seq_targets, set(Dataset)
