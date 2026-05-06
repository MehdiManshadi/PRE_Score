import requests
import os
import pandas as pd
from io import StringIO
import requests

class Region:
    BASE_URL = "https://ldlink.nih.gov/LDlinkRest/ldproxy"

    def __init__(self, token=None):
        self.token = token or os.getenv("LDLINK_TOKEN")
        if not self.token:
            raise ValueError("LDlink API token required")
        
    def ldlink_to_dataframe(self,response_text):
        df = pd.read_csv(StringIO(response_text), sep="\t")
        return df
    
    def get_ldproxy(self, rsid, population="CEU", window=50000,
                    collapse_transcripts=True, annotation="RegulomeDB",
                    ld_measure="R2", genome_build="grch38"):

        params = {
            "var": rsid,
            "pop": population,
            "r2_d": ld_measure,
            "window": str(window),
            "collapsed": "yes" if collapse_transcripts else "no",
            "annot": annotation,
            "token": self.token,
            "genome_build": genome_build
        }

        response = requests.get(self.BASE_URL, params=params, timeout=30)
        response.raise_for_status()
        
        df = self.ldlink_to_dataframe(response.text)
        Coord = df[df["R2"] >= 0.5]["Coord"]
        RS_Number = df[df["R2"] >= 0.5]["RS_Number"]
        return Coord, RS_Number

    def query_regulomedb(self, coord, assembly="GRCh38"):
        url = "https://regulomedb.org/regulome-summary/"
        
        RegDB_Coord = f"{coord.split(':')[0]}:{int(coord.split(':')[1])-1}-{coord.split(':')[1]}"

        params = {
            "regions": RegDB_Coord,
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
    
    def query_ensembl(self, coord):
        chrom, pos = coord.replace("chr", "").split(":")

        server = "https://rest.ensembl.org"
        ext = f"/overlap/region/human/{chrom}:{pos}-{pos}?feature=gene"

        headers = {"Content-Type": "application/json"}

        response = requests.get(server + ext, headers=headers)
        response.raise_for_status()
        genes = response.json()

        for gene in genes:
            print(gene["id"], gene["external_name"], gene["biotype"])

        return [gene["external_name"] for gene in genes] if genes else None
    
    def count_experiments(self, data, cell_type):

        histone_marks = {
            "H2az": +1,  # associated with regulatory elements
            "H3k4me1": +1,  # enhancer associated
            "H3k04me1": +1,  # enhancer associated
            "H3k4me2": +1,  # promoter/enhancer
            "H3k04me2": +1,  # promoter/enhancer
            "H3k4me3": +1,  # promoters/transcription start
            "H3k04me3": +1,  # promoters/transcription start
            "H3k04me3Ohtam": +1,  # assume promoter/activating
            "H3k9ac": +1,  # active regulatory elements, promoters
            "H3k09ac": +1,  # active regulatory elements, promoters
            "H3k9me1": 0,  # preference for 5' end, neutral
            "H3k09me1": 0,  # preference for 5' end, neutral
            "H3k9me3": -1,  # repressive heterochromatin
            "H3k09me3": -1,  # repressive heterochromatin
            "H3k27ac": +1,  # active enhancers/promoters
            "H3k27me3": -1,  # polycomb repression
            "H3k36me3": 0,  # elongation mark, neutral
            "H3k79me2": 0,  # transcription-associated, neutral
            "H4k20me1": 0   # preference for 5' end, neutral
        }

        for experiment in data["@graph"]:
            biosample = experiment.get("biosample_ontology", {})
            term_name = biosample.get("term_name")

            if term_name == cell_type and experiment.get("target_label") is not None:
                print(cell_type, experiment.get("method"),experiment["target_label"])

            if experiment.get("method") == "chromatin state":
                if term_name == cell_type:
                    chrom_state = experiment.get("value")
                    method = experiment.get("method")
                    print(cell_type, method, chrom_state)
            
