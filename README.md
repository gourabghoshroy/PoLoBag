# PoLoBag 
POLYNOMIAL LASSO BAGGING ALGORITHM FOR SIGNED GENE REGULATORY NETWORK INFERENCE

Author : Gourab Ghosh Roy

Maintainer : Gourab Ghosh Roy <g.ghoshroy@pgr.bham.ac.uk>

<b>Code </b>: 
The Python file PoLoBag.py can be run in Python 3.x (for Python 2.x just change the part for writing output to file at the end).
Here degree d = 2 in the algorithm.
Before running change the infilename parameter as per your static gene expression filename.
The output will be a sorted edge list. The first column denotes regulator, the second column denotes target gene and the third column is a signed edge weight. 
The weights are normalized to have a maximum absolute value of 1. Edges with non-zero weights are sorted in the output list based on the absolute value of these weights. 

<b>Test Datasets </b>: 

The following test datasets are provided 

| Dataset | Type (Simulated/Real) | Organism | Genes | Directed | Measurements/gene |
| --- | --- | --- | --- |  --- | --- | 
| A | Simulated | Yeast | 200  | Yes | 200 |
| B | Simulated | Yeast | 400  | Yes | 400  |
| C | Simulated | Yeast | 500  | Yes | 500  |
| D | Simulated | *E.coli* | 500  | Yes | 500 |
| E | Simulated | *E.coli* | 650  | Yes | 650 |
| F | Real | *E.coli* | 1419  | No | 907 |
| G | Real | Human | 522  | No | 171 |
| H | Real | *E.coli* | 99  | Yes | 24 |
| I | Real | Human | 408  | Yes | 200 |

For details and references of these datasets and the corresponding networks, please see the PoLoBag paper.

The data files with the gene expression are of the format 'input/Dataset_matrix.txt'. 
Every row contains the gene name followed by measurment values. There is a header row. There are no missing values. 
If the expression dataset contains logarithmic values, taking exponential of the values before applying the algorithm can give better inference.

The ground truth networks are provided in tab separated files 'gold_standard/Dataset_goldstandard_signed.tsv'.
For the datasets with directed networks, the first column denotes the regulator and the second column denotes the target gene.
For the undirected cases these are interactions without direction. 
The third column gives the sign of the edge (+/-). Some values can be 0 (unknown signs), as in dataset G.



