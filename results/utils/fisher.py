from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats

@dataclass
class HPO:
    database: Path
    hpo_info: Path

    def __post_init__(self) -> None:
        """
        Read in the HPO database and the HPO info.
        """
        self.hpo_data = pd.read_csv(self.database, compression='gzip', sep="\t")
        self.hpo_data.set_index('-', inplace=True)
        self.hpo_info_data = pd.read_csv(self.hpo_info, sep=",")


class FisherTest:

    def create_fisher_table(self, overlap_genes: pd.Series, gwas_genes: pd.Series, hpo_genes: pd.Series) -> pd.DataFrame:
        """
        Create a 2x2 contingency table needed for fisher exact test.

        :parameters
        -----------
        overlap_genes - pd.Series
            List of genes that overlap with the HPO database and results of a prioritization method.
        gwas_genes - pd.Series
            List of significant genes
        hpo_genes - pd.Series
            List of HPO genes

        :returns
        --------
        metrix - pd.DataFrame
            2x2 contingency table
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

    def perform_fisher_exact_tests(self, hpo_data: pd.DataFrame, gene_data: pd.Series, hpo_info: pd.DataFrame) -> pd.DataFrame:
        """
        Perform fisher's exact test on the intersect of the HPO genes, and
        genes produced by a gene prioritization method.

        :parameters
        -----------
        hpo_data - pd.DataFrame
            HPO metric inside a pandas data frame
        gene_data - pd.Series
            List of genes
        hpo_info - pd.DataFrame
            A data frame containing the name of the GWAS trait, Related HPO term, and HPO ID

        :returns
        --------
        hpo_scores - pd.DataFrame
            The original hpo_info data frame with some additional information: OR and p values from the fisher exact test
            and zscores from the p values.
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