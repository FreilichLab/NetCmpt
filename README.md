<img src="dep_sign.png" width=120, height=120 align="left" />

# NetCmpt - beta version

This is a reimplementation of NetCmpt: a network-based tool for calculating the metabolic competition between bacterial species based on:

Anat Kreimer, Adi Doron-Faigenboim, Elhanan Borenstein, Shiri Freilich, **NetCmpt: a network-based tool for calculating the metabolic competition between bacterial species**, Bioinformatics, Volume **28**, Issue 16, August 2012, Pages 2195–2197, **https://doi.org/10.1093/bioinformatics/bts323**

This is a Unpublished **beta version**

## Dependencies

* [python (3.8.18)]
* [pandas (1.5.1)]
* [networkx (3.1.0)]
* [numpy (1.24.4)]

## Installation

### Download and install NetCom2 on linux (Ubuntu 20.04)

1. clone using git

``` shell

git clone https://github.com/FreilichLab/NetCmpt.git

```
2. press on Code pulldown menue and download as zip and extract


### Create virtual environment and install snakemake 

```shell


# mamba instaltion #

conda install -n base -c conda-forge mamba

cd NetCmpt

#Install virtual enviroment 

mamba env create -f NetCmpt_env.yml

#Or:

cd NetCmpt

conda env create -f NetCmpt_env.yml

```

## Input files for NetCmpt

Input file isa tab delimited lists of EC numbers separted by space

Example file: genomes.tsv

```shell

Spc_1	1.-.-.- 1.1.1.100 1.1.1.193 1.1.1.205
Spc_2	2.6.1.17 2.6.1.9 3.6.1.31 3.5.1.18

```

## Output File for NetCmpt

Effictive metabolic overlap (EMO) is as defined in original paper of NetCmpt.

Essential compounds are listed in DB['color_index.txt'].

An additional measure is calculated, Metabolic Overlap (MO),

Which is the same as EMO yet claculated from all compounds generted by the simulation procees.

MO is the fraction of all compounds Instead of fraction of essential compounds (EMO)

Output is a tab-delimited file listing for each pair of spices emo and mo calculations.

the format is (species_x)(tab)(species_y)(tab)(emo)(tab)(mo)

Example compete.tsv:

```shell

index	species_x	species_y	emo	mo
0	Spc_1	Spc_1	1.0	1.0
1	Spc_1	Spc_2	0.11864406779661019	0.09708737864077666
2	Spc_2	Spc_1	0.33333333333333337	0.35357142857142854
3	Spc_2	Spc_2	1.0	1.0

```

## Usage of NetCmpt

python NetCmpt.py (BaseFolder path) (genomes file) (analysisMode) (output file)

## Execute NetCmpt

```shell

# activate NetCmpt enviroment #

conda activate NetCmpt_env

# change directory to NetCom2_snakemake (working direcotry) #

cd NetCmpt 

python NetCmpt.py ./ ./genomes.tsv single-species compete.tsv

```

## Contributors

[Gon Carmi](https://www.freilich-lab.com/members) \
[Shiri Freilich](https://www.freilich-lab.com/) 

## References

Anat Kreimer, Adi Doron-Faigenboim, Elhanan Borenstein, Shiri Freilich, **NetCmpt: a network-based tool for calculating the metabolic competition between bacterial species**, Bioinformatics, Volume **28**, Issue 16, August 2012, Pages 2195–2197, **https://doi.org/10.1093/bioinformatics/bts323**

Tal O, Bartuv R, Vetcos M, Medina S, Jiang J, Freilich S: **NetCom: A Network-Based Tool for Predicting Metabolic Activities of Microbial Communities Based on Interpretation of Metagenomics Data**. Microorganisms 2021, 9(9):1838.
