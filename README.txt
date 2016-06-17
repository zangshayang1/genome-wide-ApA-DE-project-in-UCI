# When I worked with Dr. Xie and Dr. Shi on a bioinformatics project about "Alternative PolyA Differential Splicing",
# I wrote quite a few scripts for data-analysis purpose. 
# Then I put them together into a command line tool that can be utilized to streamline the data process in similar projects.
# This is my first python project. 

# Please check --help page of each of the follow functions for detailed usage

usage: apap [-h]
            {IMfilter,PARcounter,MergeAPAL,MLCleaner,PARCSummer,rc2pct,calcRED,select2,hister}
            ...

positional arguments:
  {IMfilter,PARcounter,MergeAPAL,MLCleaner,PARCSummer,rc2pct,calcRED,select2,hister}
    IMfilter            IMfilter removes internally mis-primed reads
    
    PARcounter          PARcounter increment counts of reads at given APA
                        positions.
                        
    MergeAPAL           MergeAPAL merges a list of results from PARcounter
                        into one read counts (RC) table.
                        
    MLCleaner           MLCleaner removes the ApAs that: 1. are
                        "upstreamAntisense"; 2. are not associated with any
                        geneNames; 3. are within a certain range of each
                        other; 4. are repeated.
                        
    PARCSummer          PARCSummer calculates the sums of RC values of the
                        same gene.
                        
    rc2pct              rc2pct converts read counts to percentages using the
                        sums of RC values of the same gene.
                        
    calcRED             calcRED defines the approximal and distal pAs from the
                        top2 selected read counts table and generates a RED
                        score table.
                        
    select2             This function selects the top2 enriched ApAs per gene
                        from a standard RCtable
                        
    hister              This function takes a labeled ApA logFC lists and
                        outputs three types of data: 1. sliced fasta sequences
                        of a certain length centered at ApA sites for the
                        purpose of generating weblogos. 2. the counting data
                        of a certain short sequence pattern appearing at the
                        each position of the slicing window, for further
                        plotting purpose 3. the preliminary plots of the
                        counting data.

optional arguments:
  -h, --help            show this help message and exit


# Install
Download this package.
$python setup.py install

# Requirement
matplotlib (1.5.1)
numpy (1.10.4)
pathlib (1.0.1)
pysam (0.9.0)
# and other standarad libraries Anaconda-Python2.7 will suffice. 
