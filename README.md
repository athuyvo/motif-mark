# Motif-Mark

The program ```motif_mark_oop.py``` reads in a fasta file and motif text file to visualize exons and motifs on a given gene.   
Output a ```.png``` file with unique colors for each motif. 

# Capability 
- This program is capable of handing motifs with ambiguous bases. 
- Multiple sequences with a max of 10 
- Max 5 motif sequences 
- Will display overlapping motifs 
- Motifs are drawn to scale 


# Requirements 

- Must have pycairo installed 
- Input fasta file ```-f``` and input motif text file ```m``` 
- Motif file must be on a single line  
- Must have bioinfo.py module in the same directory 

# Example Usage 

```
./motif-mark-oop.py -f Figure_1.fasta -m Fig_1_motifs.txt 
```