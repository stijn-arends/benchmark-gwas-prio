"""
Perform multiple fisher exact tests for a number of prioritization methods.
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats

from arg_parser import ArgumentParser, CLIArgValidator
from results.prioritization_methods import Downstreamer, Magma, Depict, PoPs, NetWAS

__author__ = "Stijn Arends"
__version__ = "v0.1"
__data__ = "9-8-2022"


@dataclass
class HPO:
    database: Path
    hpo_info: Path

    def __post_init__(self) -> None:
        self.hpo_data = pd.read_csv(self.database, compression='gzip', sep="\t")
        self.hpo_data.set_index('-', inplace=True)
        self.hpo_info_data = pd.read_csv(self.hpo_info, sep=",")


class FisherTest:

    def __init__(self):
        pass


    def create_fisher_table(self, overlap_genes, gwas_genes, hpo_genes):
        """
        Create a 2x2 contingency table needed for fisher exact test
        """
        tl = overlap_genes[~overlap_genes.isin(gwas_genes) & ~overlap_genes.isin(hpo_genes)].shape[0]
        bl = overlap_genes[overlap_genes.isin(gwas_genes) & ~overlap_genes.isin(hpo_genes)].shape[0]
        tr = overlap_genes[~overlap_genes.isin(gwas_genes) & overlap_genes.isin(hpo_genes)].shape[0]
        br = overlap_genes[overlap_genes.isin(gwas_genes) & overlap_genes.isin(hpo_genes)].shape[0]
        
        total = tl + bl + tr + br
        
        metrix = pd.DataFrame({"No HPO": [tl, bl, tl + bl], "Yes HPO": [tr, br, tr + br],
                        "sum": [tl + tr, bl + br, total]})
        metrix.index = ["No GWAS", "Yes GWAS", "sum"]
        return metrix

    def perform_fisher_exact_tests(self, hpo_data, gene_data, hpo_info):
        """
        Perform fisher's exact test on the intersect of the HPO genes, and
        genes produced by a gene prioritization method.
        """
        hpo_scores = hpo_info.copy()
        odds_ratios = []
        p_values = []
        for hpo in hpo_info["HPO ID"]:

            try: 
                hpo_genes = hpo_data.loc[hpo_data[hpo] == 1, hpo].index
                fisher_data = self.create_fisher_table(hpo_data.index, gene_data, hpo_genes).iloc[0:2, 0:2].values

                odds_ratio, p_val = stats.fisher_exact(fisher_data)
            except KeyError:
                odds_ratio, p_val = np.nan, np.nan
                
            odds_ratios.append(odds_ratio)
            p_values.append(p_val)
            
        hpo_scores["OR"] = odds_ratios
        hpo_scores["pvalues"] = p_values
        zscores = stats.norm.ppf(p_values)
        hpo_scores["zscores"] = np.where(zscores == np.inf, 4, zscores)
    
        return hpo_scores


def main():
    
    arg_parse = ArgumentParser()

    config = arg_parse.get_argument("c")
    method = arg_parse.get_argument("m")
    output_dir = arg_parse.get_argument("o")

    cli_validator = CLIArgValidator()

    cli_validator.validate_input_file(config)

    hpo_data = r"C:\Users\stijn\Documents\Master_DSLS\Semester_two\project\HPO\phenotype_to_genes_V1268_OMIMandORPHA.txt_matrix.txt.gz"
    hpo_info = r"C:\Users\stijn\Documents\Master_DSLS\Semester_two\project\comparing_methods\HPO_table.csv"

    methods = {"NetWAS": NetWAS, "PoPs": PoPs, "DEPICT": Depict, "Downstreamer": Downstreamer, "MAGMA":Magma}



if __name__ == "__main__":
    main()