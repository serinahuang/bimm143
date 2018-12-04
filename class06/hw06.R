library(bio3d)
# Input: the four letter PDB identifier of the protein.

# Function: protein_drug() trims the protein data, keeping residues with an
# alpha chain and the CA element. It then plots the temperature factor along
# each amino acid residue in a line scatter plot, with an optional secondary
# structure estimation in the margins. To turn off margins, set margins = NULL
# in function argument.

# Output: none.

protein_drug <- function(x, margins = TRUE) {
  
  # load data of user specified protein
  protein <- read.pdb(x)
  
  # only keep amino acids in the alpha chain and of the CA element type
  chainA <- trim.pdb(protein, chain = "A", elety = "CA")
  
  # among amino acids chosen, select the column for temperature factor
  temp_fact <- chainA$atom$b
  
  # scatter plot temperature factor and secondary structure estimation
  # in the margins as default
  if (margins == TRUE) {
    plotb3(temp_fact, sse = chainA, type = "l", ylab = "Temperature Factor")
  }
  else {
    plotb3(temp_fact, sse = NULL, type = "l", ylab = "Temperature Factor")
  }
  
}

# Example:
x = "1E4Y"
protein_drug(x)