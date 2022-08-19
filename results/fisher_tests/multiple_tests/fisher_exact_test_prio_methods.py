"""
Perform multiple fisher exact tests for a number of prioritization methods.
"""

import sys
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats
import yaml
import os

root_dir = os.path.abspath(os.path.join(
                  os.path.dirname(__file__), 
                  os.pardir,
                  os.pardir))

sys.path.insert(0, root_dir)

from arg_parser import ArgumentParser, CLIArgValidator
from utils.prioritization_methods import Downstreamer, Magma, Depict, PoPs, NetWAS
from utils.fisher import HPO, FisherTest

__author__ = "Stijn Arends"
__version__ = "v0.1"
__data__ = "9-8-2022"


# @dataclass
# class HPO:
#     database: Path
#     hpo_info: Path

#     def __post_init__(self) -> None:
#         """
#         Read in the HPO database and the HPO info.
#         """
#         self.hpo_data = pd.read_csv(self.database, compression='gzip', sep="\t")
#         self.hpo_data.set_index('-', inplace=True)
#         self.hpo_info_data = pd.read_csv(self.hpo_info, sep=",")


# class FisherTest:

#     def create_fisher_table(self, overlap_genes: pd.Series, gwas_genes: pd.Series, hpo_genes: pd.Series) -> pd.DataFrame:
#         """
#         Create a 2x2 contingency table needed for fisher exact test.

#         :parameters
#         -----------
#         overlap_genes - pd.Series
#             List of genes that overlap with the HPO database and results of a prioritization method.
#         gwas_genes - pd.Series
#             List of significant genes
#         hpo_genes - pd.Series
#             List of HPO genes

#         :returns
#         --------
#         metrix - pd.DataFrame
#             2x2 contingency table
#         """
#         tl = overlap_genes[~overlap_genes.isin(gwas_genes) & ~overlap_genes.isin(hpo_genes)].shape[0]
#         bl = overlap_genes[overlap_genes.isin(gwas_genes) & ~overlap_genes.isin(hpo_genes)].shape[0]
#         tr = overlap_genes[~overlap_genes.isin(gwas_genes) & overlap_genes.isin(hpo_genes)].shape[0]
#         br = overlap_genes[overlap_genes.isin(gwas_genes) & overlap_genes.isin(hpo_genes)].shape[0]
        
#         total = tl + bl + tr + br
        
#         metrix = pd.DataFrame({"No HPO": [tl, bl, tl + bl], "Yes HPO": [tr, br, tr + br],
#                         "sum": [tl + tr, bl + br, total]})
#         metrix.index = ["No GWAS", "Yes GWAS", "sum"]
#         return metrix

#     def perform_fisher_exact_tests(self, hpo_data: pd.DataFrame, gene_data: pd.Series, hpo_info: pd.DataFrame) -> pd.DataFrame:
#         """
#         Perform fisher's exact test on the intersect of the HPO genes, and
#         genes produced by a gene prioritization method.

#         :parameters
#         -----------
#         hpo_data - pd.DataFrame
#             HPO metric inside a pandas data frame
#         gene_data - pd.Series
#             List of genes
#         hpo_info - pd.DataFrame
#             A data frame containing the name of the GWAS trait, Related HPO term, and HPO ID

#         :returns
#         --------
#         hpo_scores - pd.DataFrame
#             The original hpo_info data frame with some additional information: OR and p values from the fisher exact test
#             and zscores from the p values.
#         """
#         hpo_scores = hpo_info.copy()
#         odds_ratios = []
#         p_values = []
#         for hpo in hpo_info["HPO ID"]:

#             try: 
#                 hpo_genes = hpo_data.loc[hpo_data[hpo] == 1, hpo].index
#                 fisher_data = self.create_fisher_table(hpo_data.index, gene_data, hpo_genes).iloc[0:2, 0:2].values

#                 odds_ratio, p_val = stats.fisher_exact(fisher_data)
#             except KeyError:
#                 odds_ratio, p_val = np.nan, np.nan
                
#             odds_ratios.append(odds_ratio)
#             p_values.append(p_val)
            
#         hpo_scores["OR"] = odds_ratios
#         hpo_scores["pvalues"] = p_values
#         zscores = stats.norm.ppf(p_values)
#         hpo_scores["zscores"] = np.where(zscores == np.inf, 4, zscores)
    
#         return hpo_scores


def get_config(file: Path) -> dict:
    """
    Read in config file and return it as a dictionary.
    
    :parameter
    ----------
    file - str
        Configuration file in yaml format

    :returns
    --------
    config - dict
        Configuration file in dictionary form.
    """
    if not file.exists():
        raise FileExistsError(f"The file that was supplied does not exists: {file}")

    with open(file, 'r') as stream:
        config = yaml.safe_load(stream)

    return config


def read_hpo_info(hpo_info):
    hpo_info_data = pd.read_csv(hpo_info, sep=",")
    return hpo_info_data


def make_out_dir(path: Path) -> None:
    """
    Create a directory (if it does not exsit yet) to store the 
    data.

    :parameter
    ----------
    path - Path
        Location of directory

    :Excepts
    --------
    FileExistsError
        The directory already exists
    """
    try:
        path.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(f"[{make_out_dir.__name__}] {path} already exists.")


def write_out_data(data: pd.DataFrame, file: Path) -> None:
    """
    Write out a pandas data frame to a file.

    :parameters
    -----------
    data - pd.DataFrame
        Data
    file - Path
        Name and location of output file
    """
    data.to_csv(file, sep="\t")


def main():
    
    arg_parse = ArgumentParser()

    config_file = arg_parse.get_argument("c")
    method = arg_parse.get_argument("m")
    output_dir = arg_parse.get_argument("o")

    cli_validator = CLIArgValidator()
    cli_validator.validate_input_file(config_file)

    config = get_config(Path(config_file))

    out_dir = Path(output_dir) / method

    make_out_dir(out_dir)

    hpo_data = config["hpo_data"]
    hpo_info_data = read_hpo_info(Path(config["hpo_info"]))

    methods = {"NetWAS": NetWAS, "PoPs": PoPs, "DEPICT": Depict, "Downstreamer": Downstreamer, "MAGMA":Magma}

    hpo = HPO(database=hpo_data)
    fisher = FisherTest()

    method_instance = methods[method](hpo=hpo, fisher=fisher)

    for trait, file in config["traits"].items():
        print(f"Processing trait: {trait}")
        file = Path(file)

        method_data, genes = method_instance.read_data(file)

        overlap_hpo, overlap_genes, total_overlap = method_instance.get_overlap(hpo.hpo_data, genes)

        overlap_method = method_instance.get_overlap_genes(method_data, overlap_genes)

        sig_data_method, sig_genes = method_instance.filter_data(overlap_method)

        fish_results = method_instance.fisher.perform_fisher_exact_tests(overlap_hpo, sig_genes, hpo_info_data)
        
        out_file = out_dir / ("fisher_result_" + file.stem + ".csv")

        write_out_data(fish_results, out_file)


if __name__ == "__main__":
    main()
