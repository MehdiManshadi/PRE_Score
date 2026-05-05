import requests
import os
import pandas as pd
from io import StringIO

class LDlinkClient:
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
                    ld_measure="R2"):

        params = {
            "var": rsid,
            "pop": population,
            "r2_d": ld_measure,
            "window": str(window),
            "collapsed": "yes" if collapse_transcripts else "no",
            "annot": annotation,
            "token": self.token,
        }

        response = requests.get(self.BASE_URL, params=params, timeout=30)
        response.raise_for_status()
        
        df = self.ldlink_to_dataframe(response.text)
        LDLinked = df[df["R2"] <= 0.5]["Coord"]

        return LDLinked