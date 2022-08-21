#!/bin/bash
: '
A bash script to download and manipulate features needed to run PoPs.

Author - Stijn Arends, Raana Roohanitaziani
Date - 15-6-2022
contact - s.arends@st.hanze.nl, r.roohanitaziani@st.hanze.nl
'

repo=https://github.com/FinucaneLab/gene_features.git


# Check if the repository is already cloned
if [ ! -d "gene_features" ]; then
   git clone ${repo}
fi


cd gene_features

# Remove all files and folders except for the features dir
ls | grep -v "features" | xargs rm -r

cd features

# check if there is data to extract (i.e. are there any directories with data)
vals=($(find . -maxdepth 1 -type d))

# Check if there are no directories
# If there are no directories, then there is no data to extract, therefore we exit the program
if [ -z "${vals}" ]; then
    echo "No directories to extract data from, exiting..."
    exit
fi

# Extract all the data and add a prefix
for dir in ./*/     #
do
    dir=${dir%*/}      # remove the trailing "/"
    cd ${dir##*/}       # everything after the final "/"

    for f in * ; do mv "$f" ../${dir##*/}_"$f" ; done
    cd ..
    rm -r ${dir##*/}
done

gunzip -f *.gz

cd ../../

# Rename the files and columns by adding a prefix.
python3 rename_file_contents.py -d gene_features/features/ --prefix pops_