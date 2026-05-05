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
        
        RegDB_Coord = [f"{c.split(':')[0]}:{int(c.split(':')[1])-1}-{c.split(':')[1]}" for c in coord]

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

        return genes["id"] if genes else None