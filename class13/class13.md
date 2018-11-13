Class 13: Structural Bioinformatics II
================
Serina Huang
November 13, 2018

*in silico* Docking of Drugs to HIV-1 Protease
----------------------------------------------

### Section 1. Preparation of Protein and Ligand

``` r
library(bio3d)
file.name <- get.pdb("1HSG")
```

    ## Warning in get.pdb("1HSG"): ./1HSG.pdb exists. Skipping download

``` r
hiv <- read.pdb(file.name)
```

Q1. What are the 2 non-protein **resid** values in this structure? What does **resid** correspond to and how would you get a listing of all residue values in this structure?

``` r
# Quick summary of contents of this pdb structure
hiv
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

A1. HOH, water, and MK1, a Merck drug compound.

##### Split into separate protein and ligand files

``` r
# The protein-only portion is our receptor
prot <- trim.pdb(hiv, "protein")
prot
```

    ## 
    ##  Call:  trim.pdb(pdb = hiv, "protein")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
# Write this into a file
write.pdb(prot, file = "1hsg_protein.pdb")
# Inspect the resulting file that all are proteins
```

Do the same for ligand:

``` r
lig <- trim.pdb(hiv, "ligand")
lig
```

    ## 
    ##  Call:  trim.pdb(pdb = hiv, "ligand")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
# The only non-protein residue is MK1, and no HOH
write.pdb(lig, file = "1hsg_ligand.pdb")
# Inspect again
```

Note: If there is more than 1 ligand in a structure, you can get the ligand at a certain residue, e.g. at 902.

##### Prepare the protein

After loading in the protein molecule in AutoDockTools (ADT), add hydrogen atoms, which are smaller than 2 Angstroms and thus not displayed.

Grid keeps the polar hydrogens, i.e. only the hydrogens that can form hydrogen bonds. Then, ADT saves this file in .pdbqt format, which is akin to a .pdb file but with charge (e.g. -0.062) and type (e.g. Nitrogen). Grid Box lets you select where you want the ligand to fit. Usually, people box the whole moelcule.

The following code, executed in Terminal, can move the .pdbqt file to your desired directory.

    # Tilde means home area, back slack means the separator
    ls ~/
    # . means put it right here in this folder
    mv ~/1hsg_protein.pdbqt

Likewise, prepare the ligand in ADT.

### Section 2. Docking

    # Specify the protein and ligand files to fit
    receptor = 1hsg_protein.pdbqt
    ligand = 1hsg_ligand.pdbqt

    # Try 50 times (poses)
    num_modes = 50
    out = all.pdbqt

    # The coordinates restrain our ligand to dock at the center of the receptor
    # Part of the ligand can be outside the box
    # But restrain makes the program run faster
    center_x = 16
    center_y = 25
    center_z = 4
    size_x = 30
    size_y = 30
    size_z = 30

    # Random seed (optional)
    seed = 2009

    # The program will now calculate the configuration with the lowest energy

After downloading Autodock Vina, you can do the actual docking! This might take a while and your computer might start whirring.

    ~/Downloads/autodock_vina_1_1_2_mac/bin/vina --config config.txt --log log.txt

Food for thought: Would the same output be generated if a different seed were used?

To visualize the results of `all.pdbqt`, we must convert it to a regular .pdb format.

``` r
# This file has multiple poses, so we need an extra argument
res <- read.pdb("all.pdbqt", multi = TRUE)
write.pdb(res, "results.pdb")
res
```

    ## 
    ##  Call:  read.pdb(file = "all.pdbqt", multi = TRUE)
    ## 
    ##    Total Models#: 14
    ##      Total Atoms#: 50,  XYZs#: 2100  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 50  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, xyz, calpha, call

Now, we want to visualize the docking on VMD and compare to Merck's drug.

Open `1HSG.pdb` -&gt; Graphics -&gt; Representations:

Default Drawing Method is Lines.

    Visualize protein: NewCartoon, Secondary Structure, all
    Visualize ligand: Licorice, ColorID 8, resname MK1

Open `results.pdb` -&gt; Graphics -&gt; Representations:

Note: You can double click T to change which molecule is on top.

There are 14 different results, so you can select the pose on the slider. E.g. the one with -11.5 affinity is index 0 and the best. You can also play through the results on a low speed. Notice that index 1 is a mere 180º mirror image of index 0.

Create a new representation under `1HSG.pdb`:

    Surf, Secondary Structure, resname protein

Now, we see that our index 0 matches pretty well with Merck's compound, with the exception of the functional groups at the end, which are not interacting with the protein (receptor). This suggests that we can improve the binding further by changing the functional group to one that binds to the receptor.
