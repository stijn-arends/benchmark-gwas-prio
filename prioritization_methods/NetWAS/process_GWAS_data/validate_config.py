"""
This module is designed to validate the contents of the configuration file associated 
with the script processing the GWAS summstats files for VEGAS. 
"""

import yaml
import linecache
from pathlib import Path


class ConfigValidator:

    def __init__(self, config):
        self.config = self.get_config(config)

    def get_config(self, file):
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
        with open(Path(file), 'r') as stream:
            config = yaml.safe_load(stream)
        return config

    def validate_config_file(self) -> None:
        self.validate_input_exists()
        self.validate_column_names()

    def validate_input_exists(self) -> None:
        """
        Check if a file exists. 
        """
        for info in self.config["traits"].values():
            if not Path(info["file"]).is_file():
                raise FileNotFoundError(f'Input file does not exist:\n{info["file"]}')

    def validate_output_exists(self) -> None:
        """
        Check if the output location exists. 
        """
        file = Path(self.config["output"])
        if not file.exists():
            file.mkdir(parents=True)

    def validate_column_names(self) -> None:
        """
        Validate that the given names for the SNP and pvalue columns 
        are actually inside of the corresponding file. 
        """
        for trait, info in self.config["traits"].items():
            print(trait)
            snp = info["columns"]["snp"]
            pval = info["columns"]["pval"]

            if snp == "None" or pval == "None":
                continue

            columns = linecache.getline(info["file"], 1).rstrip().split("\t")

            if not set([snp, pval]).issubset(columns):
                print(f"WARNING! The set value for the SNP and/or p value ({snp}, {pval}) column are not found amongst the column names of the corresponding file." \
                    " The program will try to automatically detect the correct column.")



def main():
    validator = ConfigValidator("config.yaml")
    file = "C:\\Users\\stijn\\Documents\\Master_DSLS\\Semester_two\\project\\data\\prostate_cancer\\meta_v3_onco_euro_overall_ChrAll_1_release_full.txt"
    # print(linecache.getline(file, 1))


    validator.validate_config_file()


if __name__ == "__main__":
    main()