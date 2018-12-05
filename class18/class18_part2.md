Class 17: Cancer Genomics Part 2
================
Serina Huang
12/04/18

### Section 1. Protein Sequences from Healthy and Tumor Tissues

The following sequences are the results of an NGS analysis of a healthy vs. tumor tissue in one patient.

Q1. Identify sequence regions that contain all 9-mer peptides that are only found in the tumor.

A1. First, find the mutation sites by alignment (either by the web-based MUSCLE or the R `bio3d` package).

``` r
library(bio3d)
p53seq <- read.fasta("https://bioboot.github.io/bimm143_F18//class-material/lecture18_sequences.fa")
p53seq
```

    ##              1        .         .         .         .         .         60 
    ## P53_wt       MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
    ## P53_mutant   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMLDLMLSPDDIEQWFTEDPGP
    ##              **************************************** ******************* 
    ##              1        .         .         .         .         .         60 
    ## 
    ##             61        .         .         .         .         .         120 
    ## P53_wt       DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ## P53_mutant   DEAPWMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ##              **** ******************************************************* 
    ##             61        .         .         .         .         .         120 
    ## 
    ##            121        .         .         .         .         .         180 
    ## P53_wt       SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ## P53_mutant   SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ##              ************************************************************ 
    ##            121        .         .         .         .         .         180 
    ## 
    ##            181        .         .         .         .         .         240 
    ## P53_wt       RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ## P53_mutant   RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFVHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ##              ******************************** *************************** 
    ##            181        .         .         .         .         .         240 
    ## 
    ##            241        .         .         .         .         .         300 
    ## P53_wt       SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
    ## P53_mutant   SCMGGMNRRPILTIITLEV-----------------------------------------
    ##              ******************                                           
    ##            241        .         .         .         .         .         300 
    ## 
    ##            301        .         .         .         .         .         360 
    ## P53_wt       PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
    ## P53_mutant   ------------------------------------------------------------
    ##                                                                           
    ##            301        .         .         .         .         .         360 
    ## 
    ##            361        .         .         .  393 
    ## P53_wt       GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
    ## P53_mutant   ---------------------------------
    ##                                                
    ##            361        .         .         .  393 
    ## 
    ## Call:
    ##   read.fasta(file = "https://bioboot.github.io/bimm143_F18//class-material/lecture18_sequences.fa")
    ## 
    ## Class:
    ##   fasta
    ## 
    ## Alignment dimensions:
    ##   2 sequence rows; 393 position columns (259 non-gap, 134 gap) 
    ## 
    ## + attr: id, ali, call

We want to find a 9-mer sequence that fits into the pocket in the MLC molecule; and within that 9-mer, there must be a mutation site that's specific to the tumor.

Each asterisk denotes a match. We see that the end of the mutant sequence is truncated. We can also eyeball a few other mutation sites, but this approach is error-prone and unproductive. It's also problematic if there are more than one mutation site and/or more than one patient.

Computationally, we can extract all possible 9-mer sequences containing a mutation site. In other words, the mutation can either be the first or the last of the 9-mer.

The `conserv()` function in `bio3d` that searches for conserved sequences. What we want, however, is the "non-conserved" sequences (of the tumor).

``` r
score <- conserv(p53seq, method = "identity")
# Find the scores that are less than 1 and use which to get their indices
mut.ind <- which( (score < 1) )
# Verify the indices from the original fasta file
p53seq$ali[,1]
```

    ##     P53_wt P53_mutant 
    ##        "M"        "M"

``` r
p53seq$ali[,2]
```

    ##     P53_wt P53_mutant 
    ##        "E"        "E"

``` r
p53seq$ali[,3]
```

    ##     P53_wt P53_mutant 
    ##        "E"        "E"

``` r
p53seq$ali[,41]
```

    ##     P53_wt P53_mutant 
    ##        "D"        "L"

``` r
p53seq$ali[,mut.ind]
```

    ##            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## P53_wt     "D"  "R"  "R"  "D"  "S"  "S"  "G"  "N"  "L"  "L"   "G"   "R"  
    ## P53_mutant "L"  "W"  "V"  "V"  "-"  "-"  "-"  "-"  "-"  "-"   "-"   "-"  
    ##            [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
    ## P53_wt     "N"   "S"   "F"   "E"   "V"   "R"   "V"   "C"   "A"   "C"  
    ## P53_mutant "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"  
    ##            [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32]
    ## P53_wt     "P"   "G"   "R"   "D"   "R"   "R"   "T"   "E"   "E"   "E"  
    ## P53_mutant "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"  
    ##            [,33] [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42]
    ## P53_wt     "N"   "L"   "R"   "K"   "K"   "G"   "E"   "P"   "H"   "H"  
    ## P53_mutant "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"  
    ##            [,43] [,44] [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52]
    ## P53_wt     "E"   "L"   "P"   "P"   "G"   "S"   "T"   "K"   "R"   "A"  
    ## P53_mutant "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"  
    ##            [,53] [,54] [,55] [,56] [,57] [,58] [,59] [,60] [,61] [,62]
    ## P53_wt     "L"   "P"   "N"   "N"   "T"   "S"   "S"   "S"   "P"   "Q"  
    ## P53_mutant "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"  
    ##            [,63] [,64] [,65] [,66] [,67] [,68] [,69] [,70] [,71] [,72]
    ## P53_wt     "P"   "K"   "K"   "K"   "P"   "L"   "D"   "G"   "E"   "Y"  
    ## P53_mutant "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"  
    ##            [,73] [,74] [,75] [,76] [,77] [,78] [,79] [,80] [,81] [,82]
    ## P53_wt     "F"   "T"   "L"   "Q"   "I"   "R"   "G"   "R"   "E"   "R"  
    ## P53_mutant "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"  
    ##            [,83] [,84] [,85] [,86] [,87] [,88] [,89] [,90] [,91] [,92]
    ## P53_wt     "F"   "E"   "M"   "F"   "R"   "E"   "L"   "N"   "E"   "A"  
    ## P53_mutant "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"  
    ##            [,93] [,94] [,95] [,96] [,97] [,98] [,99] [,100] [,101] [,102]
    ## P53_wt     "L"   "E"   "L"   "K"   "D"   "A"   "Q"   "A"    "G"    "K"   
    ## P53_mutant "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"    "-"    "-"   
    ##            [,103] [,104] [,105] [,106] [,107] [,108] [,109] [,110] [,111]
    ## P53_wt     "E"    "P"    "G"    "G"    "S"    "R"    "A"    "H"    "S"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,112] [,113] [,114] [,115] [,116] [,117] [,118] [,119] [,120]
    ## P53_wt     "S"    "H"    "L"    "K"    "S"    "K"    "K"    "G"    "Q"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,121] [,122] [,123] [,124] [,125] [,126] [,127] [,128] [,129]
    ## P53_wt     "S"    "T"    "S"    "R"    "H"    "K"    "K"    "L"    "M"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,130] [,131] [,132] [,133] [,134] [,135] [,136] [,137] [,138]
    ## P53_wt     "F"    "K"    "T"    "E"    "G"    "P"    "D"    "S"    "D"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"

We just want the first four tumor-specific variants because we can't present "gaps" to the MHC molecule. We can take -8 and +8 positions of the mutant to get the series of 9-mers.

However, to do so, we need to get rid of the gaps first. Thankfully, there's already a `bio3d` function for it.

``` r
gaps <- gap.inspect(p53seq)
# Indices for true gaps
notgaps.ind <- gaps$f.inds

# Now we have mut.ind and notgaps.ind
# Intersect the mutant sites vs. non-gaps
# The intersect function is %in%, i.e. in which of the index is mut.ind in notgaps.ind?
# Try in console: c(1:5) %in% c(4:10)
# Does dimension need to equal? Nope!
mut.ind_notgaps <- mut.ind[mut.ind %in% notgaps.ind]
mut.ind_notgaps
```

    ## [1]  41  65 213 259

Now we can get the 8 before and 8 after positions for the 9-mers.

``` r
start.position <- mut.ind_notgaps - 8
end.position <- mut.ind_notgaps + 8

start.position
```

    ## [1]  33  57 205 251

``` r
end.position
```

    ## [1]  49  73 221 267

Let's try to get our first 9-mer.

``` r
p53seq$ali[ ,start.position[1]:end.position[1] ]
```

    ##            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## P53_wt     "S"  "P"  "L"  "P"  "S"  "Q"  "A"  "M"  "D"  "D"   "L"   "M"  
    ## P53_mutant "S"  "P"  "L"  "P"  "S"  "Q"  "A"  "M"  "L"  "D"   "L"   "M"  
    ##            [,13] [,14] [,15] [,16] [,17]
    ## P53_wt     "L"   "S"   "P"   "D"   "D"  
    ## P53_mutant "L"   "S"   "P"   "D"   "D"

Aha! We confirm that at position 9, there's one of our mutation sites.
