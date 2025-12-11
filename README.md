This is the project folder containing my custom snakemake and R scripts for the Streptococcus mutans project, to facilitate transparency and reproduceability 

Main workflows:
- Reference Mapping v2: Workflow which maps samples against a reference panel and helps the user ascertain the best fit through mapping statistics, plots and tables
- Strep_anvio: Workflow based on Anv'io software suite which assembles and bins samples and creates an interactive window for manual refinement 
- Strep_Phylogenetics 05 23: Workflow for phylogenetic tree construction.

Two smaller scripts for finding a representative subset of genomes based on ANI are here:
- fastANIdb: Workflow for calculating ANI between genomes
- Strep_Phylogenetics_09-21: Workflow that calculates a representative subset and does a phylogenetic plot of a distance matrix (using fastANIdb output)

Other scripts:
- Strep_pangenome: pangenome analysis of ancient MAGs based on An'vio
- misc_scripts: Miscellaneous R scripts used for analysis. paper_plots.R is the script which collects all the code used in creating figures and tables included in the manuscript
