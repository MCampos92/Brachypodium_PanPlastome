# Indels and Heteroplasmy

In this section, you can find the scripts to detect insertions/deletions in the alignment and parse it with the heteroplasmic positions.

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium* plastoems data included in the paper:
>
> "The *Brachypodium* panplastome reveals extended heteroplasmy and positive selection in the chloroplast genomes of its three annual model grasses" by Miguel Campos Cáceres.
>
> And co-authored by Ernesto Pérez-Collazos, John Vogel and Pilar Catalán. 

## Indels detection
To detect the indels belonging to the plastomes of all individuals we need:
- `Plastomes.fasta` that contains the plastome sequences aligned.
- `Indels.R` is the main script with all the commands used in R.

## Parsing heteroplasmy with the indels
To detect if the heteroplasmic positions belongs or not to the indels detected we used:
- `Heteroplasmy.csv` heteroplasmic positions detected using NovoPlasty.
- `Indels.csv` main output of the previous script.
- `Parse.py` is the python script used to parse the positions.

## Inferring the origin of the indels
To inferr the origin of the heteroplasmic positions in B. hybrdidum we need:
- `Origin.R` main script to split into one parental or another using R.
- `Heteroplasmy.csv` discarding the heteroplasmic positions of the parentals.
- Consensus sequence of both parentals.
