"""
Perform fisher's exact tests on the results from a number of prioritization methods.

Idea:
create a function that can write out the fisher results to a file and a function that will print the results based on which mode the user selected.

to-do:
Add documentation to the methods/functions
"""

import sys
from pathlib import Path
# import numpy as np
# import pandas as pd
# import scipy.stats as stats
import yaml
import os

from matplotlib_venn import venn2, venn2_circles
import matplotlib.pyplot as plt

root_dir = os.path.abspath(os.path.join(
                  os.path.dirname(__file__), 
                  os.pardir,
                  os.pardir))

sys.path.insert(0, root_dir)

from arg_parser import ArgumentParser, CLIArgValidator
from utils.prioritization_methods import Downstreamer, Magma, Depict, PoPs, NetWAS
from utils.fisher import HPO, FisherTest

class VennDiagram:

    def __init__(self, n_hpo_genes, n_method_genes, total_overlap) -> None:
        self.n_hpo_genes = n_hpo_genes
        self.n_method_genes = n_method_genes
        self.total_overlap = total_overlap

    def __plot_venn_diagram(self):
        """
        Perhaps create two seperate functions, one for saving the plot and the other for showing.
        Or add an output dir to the function, or make a class for this and save the output dir there.
        """

        # depict venn diagram
        venn2(subsets=(self.n_hpo_genes, self.n_method_genes, self.total_overlap), 
            set_labels=('HPO genes', 'GWAS genes'),
            set_colors=("silver", "lightsteelblue"), alpha=0.7)

        venn2_circles(subsets=(self.n_hpo_genes, self.n_method_genes, self.total_overlap),
                    linewidth=1)

    def show(self):
        self.__plot_venn_diagram()
        plt.show()

    def save(self, output_file):
        self.__plot_venn_diagram()
        plt.savefig(output_file)


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


def save_fisher_results(file, fisher_table, odds_ratio, pval):
    with open(file, 'w') as f:
        # fish_table = fisher_data.to_string()
        f.write("2x2 contingency table:\n")
        f.write(fisher_table)
        f.write(f"\nFisher's exact test results:")
        f.write(f"odds ratio: {odds_ratio}, pvalue: {pval}\n")


def print_fisher_results(fisher_table, odds_ratio, pval):
    print("2x2 contingency table:\n")
    print(fisher_table)
    print(f"\nFisher's exact test results:")
    print(f"odds ratio: {odds_ratio}, pvalue: {pval}\n")


def main():
    arg_parse = ArgumentParser()

    config_file = arg_parse.get_argument("c")
    method = arg_parse.get_argument("m")
    output_dir = arg_parse.get_argument("o")

    # Save or plot mode
    save_mode = arg_parse.get_argument('s')
    mode= "save" if save_mode else "plot"

    cli_validator = CLIArgValidator()
    cli_validator.validate_input_file(config_file)
    cli_validator._check_arg_combination(output_dir, save_mode)

    config = get_config(Path(config_file))

    methods = {"NetWAS": NetWAS, "PoPs": PoPs, "DEPICT": Depict, "Downstreamer": Downstreamer, "MAGMA":Magma}

    hpo_data = config["hpo_data"]
    hpo = HPO(database=hpo_data)
    fisher = FisherTest()

    method_instance = methods[method](hpo=hpo, fisher=fisher)

    for trait, info in config["traits"].items():
        print(f"Processing trait: {trait}")
        if trait == "IBD":
            file = Path(info["file"])
            hpo_term = info["hpo_term"]
            method_data, genes = method_instance.read_data(file)

            overlap_hpo, overlap_genes, total_overlap = method_instance.get_overlap(hpo.hpo_data, genes)

            overlap_method = method_instance.get_overlap_genes(method_data, overlap_genes)

            sig_data_method, sig_genes = method_instance.filter_data(overlap_method)

            data_hpo_term, genes_hpo_term = hpo.get_data_hpo_term(hpo.hpo_data, hpo_term)

            n_hpo_genes = len(hpo.hpo_data)
            # print(f"Number of hpo genes: {n_hpo_genes}")

            # print(f"Number of total netwas IBD genes: {len(genes)}")

            # print(f"Genes overlapping with HPO: {total_overlap}")

            fisher_data = fisher.create_fisher_table(overlap_genes, sig_genes, genes_hpo_term)# .iloc[0:2, 0:2].values
            fisher_table = fisher_data.to_string()

            OR, pval = fisher.fishers_exact_test(fisher_data.iloc[0:2, 0:2].values)
            # print(f"OR: {OR}, pval: {pval}")

            # with open("test.txt", 'w') as f:
            #     fish_table = fisher_data.to_string()
            #     f.write("2x2 contingency table:\n")
            #     f.write(fish_table)
            #     f.write(f"\nFisher's exact test results:\n")
            #     f.write(f"\n\nodds ratio: {OR}, pvalue: {pval}")

            venn_diagram = VennDiagram(n_hpo_genes, len(genes), total_overlap)

            if mode == "save": 
                venn_diagram.save(Path(output_dir) / (trait + "_venn.png"))
                save_fisher_results(Path(output_dir) / (trait + "_fisher_results.txt"), fisher_table, OR, pval)
            else: 
                venn_diagram.show()
                print_fisher_results(fisher_table, OR, pval)
            # plot_venn_diagram(n_hpo_genes, len(genes), total_overlap, mode)


if __name__ == "__main__":
    main()
