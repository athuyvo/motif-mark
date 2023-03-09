#!/usr/bin/env python

import cairo
import math 
import bioinfo
import argparse
import re

class Motif:

    def __init__(self, sequence, exon, start_pos, stop_pos, color):
        '''Create a motif object with the original sequence, if exon, start position, stop position, and unique color.'''

        ## Data ##
        self.seq = sequence
        self.is_exon = exon
        self.start = start_pos
        self.stop = stop_pos
        self.colors = color
        self._owner = None
    
    # Draw motifs on gene 
    def draw(self, gene_count, context, motif_positions):
        start = self.start + 20
        end = self.stop + 20
        context.move_to(20, gene_count)
        
        # Change motif size if intron
        if self.is_exon is True: 
            rec_y0 = 25
            rec_y1 = 50
        else: 
            rec_y0 = 20
            rec_y1 = 40
        
        # Change motif size if overlap 
        if self.start in motif_positions or self.stop in motif_positions: 
            rec_y0 = rec_y0 - 5
            rec_y1 = rec_y1 - 5
                
        context.rectangle(start, gene_count-rec_y0, start-end, rec_y1) #(x0,y0,x1,y1)
        context.set_source_rgb(self.colors[0],self.colors[1],self.colors[2])
        context.fill()
        motif_positions.append(self.start)
        motif_positions.append(self.stop)

    
def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--fasta", required=True, type=str, help="absolute filepath to fasta file")
    parser.add_argument("-m", "--motif", required=True, type=str, help="absolute filepath to motif file")

    return(parser.parse_args())

args = get_args()
input_file = args.fasta
motif_file = args.motif
output_name = input_file.split(".")[0]

# Motif dictionary from motif text file 
motif_dict = {}

# List of current motif positions on gene
motif_positions = [] 

# dictionary of ambiguous bases as keys and its potential base as values 
ambig_bases_dict = {
    "Y": "[YCT]", "M": "[MAC]", "K": "[KGT]",
    "W": "[WAT]", "S": "[SCG]", "R": "[RAG]", 
    "B": "[BCGT]", "D": "[DAGT]", "H": "[HACT]",
    "V": "[VACG]" , "N": "[NACGTU]", "U": "[UA]"
}

# List of RGB colors
color_list = [(255,0,0), (0,255,0),(0,0,255), (255,255,0), (0,255,255), 
           (255,0,255), (192,192,192), (128,128,128), (128,0,0), 
           (128,128,0), (0,128,0), (128,0,128), (0,128,128), (0,0,128)]


# Create the coordinates to display your graphic, desginate output
width, height = 2000, 3000
surface = cairo.PDFSurface(output_name + ".pdf", width, height)
context = cairo.Context(surface)

# Read motif file and create a dictionary of motifs (considers ambiguous bases) as keys 
# and a motif object as values 

def read_motif(file):
    ''' This function reads a given motif file and finds ambiguous bases then creates a 
    motif dictionary with potential motifs as keys and a motif object as values.'''

    for line in file:
        key = line.strip()
        temp_line = key
        upper = True 

        # Check if string is lowercase (introns), if yes, convert to upper case 
        if temp_line.islower() is True:
            temp_line = temp_line.upper()
            upper = False  

        # Check if string is in ambiguous bases dictionary, if yes, 
        # add potential bases as a tuple in motif dictionary as keys 
        for base1 in ambig_bases_dict:
            
            # Found ambigous base
            if base1 in temp_line:

                if upper is False: 
                    key = key.replace(base1.lower(), ambig_bases_dict[base1].lower())
                else: 
                    key = key.replace(base1, ambig_bases_dict[base1])
            
        # Create motif object as dictionary value 
        # Assign a color to each unique motif 
        motif_dict[key] = Motif(line, upper, 0, 0, color_list[len(color_list)-1]) 
        del color_list[-1]


def find_motif(line):
    '''Find motif using motif dictionary in a given fasta sequence. Then draws motif if motif is found '''
   
    # Find motif in sequence line and draw motif on gene 
    for key in motif_dict: 
        motif_found = re.finditer(f"(?=({key}))", line, re.IGNORECASE)      

        for motif in motif_found:

            motif_dict[key].start = motif.start()
            motif_dict[key].stop = motif.start() + len(motif_dict[key].seq)
            motif_dict[key].draw(gene_count, context, motif_positions)

# Find exons in the fasta file
def find_exon(line, gene_count):
    exons = re.finditer("[A-Z]+", line)
   
    for e in exons: 
        motif = Motif(e, True, e.start(), e.end(), (64/255,64/255,64/255))
        motif.draw(gene_count, context, motif_positions)

# Create blank canvas 
def create_canvas():
    
    # Set background to white then reset source color
    context.set_source_rgb(1,1,1)
    context.paint()
    context.set_source_rgb(0,0,0)

# Draw scaled lines for each gene(sequence) encountered in fasta file 
def draw_gene(name, gene_count, length): 

    # Need to tell cairo where to put the brush, the color and width, and the shape you want it to draw   
    context.set_source_rgb(0,0,0)

    # Write out gene name as headers
    context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    context.set_font_size(15)
    context.move_to(10, gene_count-50)
    context.show_text(name)

    # Draw a line for gene 
    context.set_line_width(5)
    context.move_to(20, gene_count)        #(x,y)
    context.line_to(length+20, gene_count)
    context.stroke()

# Draw a legend for each motif found 
def draw_legend(gene_count):

    context.set_source_rgb(0,0,0)
    context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    context.set_font_size(20)

    # Write out sequence and draw a key for each motif 
    for key in motif_dict: 

        # Draw rectangle with motif color 
        context.set_source_rgb(motif_dict[key].colors[0], motif_dict[key].colors[1], motif_dict[key].colors[2])
        context.rectangle(20, gene_count, 20, 20)
        context.fill()

        # Write out motif sequence for key 
        context.move_to(50, gene_count+15)
        context.set_source_rgb(0,0,0)
        context.show_text((motif_dict[key].seq).strip())
        gene_count += 50
    
      
    
# Convert fasta to one line fasta 
new_fasta = bioinfo.oneline_fasta(input_file)     

# Open motif file and new one lined fasta file 
with open(motif_file, "r") as motif_f, \
open("oneLineFasta.fa", "r") as fasta:
   
    # Read motif file and create a dictionary with motifs as keys 
    read_motif(motif_f)  
    create_canvas()

    # keep track of each gene and create space between each gene 
    gene_count = 0 
    header_name = ""

        
    # Parse through fasta file and find sequence 
    for line in fasta: 

        line = line.strip()

        # Grab gene name from header line 
        if line.startswith(">"): 
            header_name= line
        
        else: 
            # Get gene length to draw motifs to scale 
            gene_count += 150    
            draw_gene(header_name, gene_count, len(line))
            find_exon(line, gene_count)    
            find_motif(line)  

        # Reset motif positions for next gene 
        motif_positions = []
            
    draw_legend(gene_count+150)

surface.write_to_png(output_name + ".png") # Output to PNG

#write out drawing
surface.finish()

               
