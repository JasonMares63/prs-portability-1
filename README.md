# PRS Portability

This project attempts to probe the portability of polygenic risk scores (PRS) across and within populations.
Previous work has largely focused on manually drawn binary or categorical boundaries such as 1000 Genomes "Super population codes", inferences of "wealthy status" from UK Biobank, or manually defined geographical regions.
This work attempts to change that paradigm and represent replicability in terms of a continuum dependent on loci and overall genetic relatedness.
In a simple example, we could examine how helpful a PRS computed in a relatively homogeneous population A is for other people, arranged according to their genetic distance from the centroid of population A. 

## Data

This project uses data from four sources.
* UK Biobank
* 1000 Genomes Project
* Biobank Japan
* Neale Lab UK Biobank GWAS
    * Included for comparisons of GWAS results

## Running

This project was largely conducted on non-public computational resources and data.
Because of this, it is unlikely that the workflow could be replicated exactly.
However, each script's purpose is documented both in its name and comments within each file.
For questions, please open an issue or contact @zietzm.
