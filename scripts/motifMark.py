#!/usr/bin/env python

#--------------------------------------------------------------------------------------------------------------
# Motif Marking tool
# Author: Cameron Watson
# Last Updated: 2 Feb 2021
#--------------------------------------------------------------------------------------------------------------

import argparse
import re
import cairo
import numpy as np

#--------------------------------------------------------------------------------------------------------------
# USER INPUT
#--------------------------------------------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(description = "A tool for marking specified motifs in FASTA files. This program returns \
        an image (SVG) of genes with exons and introns denoted and specified motifs marked.",
    add_help = True)

    parser.add_argument("-f", "--fasta", 
    help = "A FASTA file containing sequences with exons capitalized and all other nucleotides in lower-case", required = True)

    parser.add_argument("-m", "--motif",
    help = "A plain-text file with one motif (<=10 base nucleotide sequence) per line", required = True)

    return parser.parse_args()

args = get_args()

#--------------------------------------------------------------------------------------------------------------
# FUNCTIONS
#--------------------------------------------------------------------------------------------------------------

def expandMotifs(motif_list):
    '''
    returns a dict of regex versions of all motifs with the motif length as value,
    also returns a dict with original motif:regex motif
    '''

    # RNA/DNA indifferent dictionary for IUPAC degenerate bases
    iupac_degens = {
        "w" : "[atuATU]", "b" : "[cgtuCGTU]", "r" : "[agAG]", "t" : "[tuTU]",
        "s" : "[cgCG]", "d" : "[agtuAGTU]", "n" : "[acgtuACGTU]", "u" : "[tuTU]",
        "m" : "[acAC]", "h" : "[actuACTU]", "y" : "[ctuCTU]", "a" : "[aA]",
        "k" : "[gtuGTU]", "v" : "[acgACG]", "z" : "-", "c" : "[cC]", "g" : "[gG]"
    }

    final_motif = {}
    maintainInput_motif = {}
    for motif in motif_list: # iterate thru all input motifs 
        motif = motif.lower() # looks for a motif in both exons and introns, case insensitive
        new_motif = ""
        for i in motif:
            new_motif += iupac_degens[i] # rebuild input motif with IUPAC regex from dict
        final_motif[new_motif] = len(motif)
        maintainInput_motif[motif] = new_motif # Ex: YYY:[ctuCTU][ctuCTU][ctuCTU]

    return final_motif, maintainInput_motif


class gene:

    def __init__(self, record, motifs):
        '''
        initializes an instance of gene. Takes a fasta record in list form: [header,seq]
        Isolates the record name, obtains coordinates for introns and exons
        '''
        self.longName = re.split(">", record[0])[1] # ditch the greater than sign
        self.name = re.split(" ", self.longName)[0] # isolate gene name
        self.seq = record[1]
        self.seqLen = len(self.seq)
        self.exons = [] # initializing, to be updated with exon coords
        self.introns = [] 
        self.motifs = motifs
        # calling functions to automatically update initialized attributes
        self.findRegions()
        self.findMotifs()

    def findRegions(self):
        '''
        called by init, returns all exon and intron coordinates using regex
        '''
        exon_matches = re.finditer("[A-Z]+", self.seq)
        intron_matches = re.finditer("[a-z]+", self.seq)
        # extracting coordinates (first char index: last char index) from re.iter objects
        for i in exon_matches:
            self.exons.append(i.span())
        for i in intron_matches:
            self.introns.append(i.span())
        return self

    def findMotifs(self):
        '''
        called by init, locates motifs within DNA/RNA sequence, returns
        dict of coordinates:motif
        '''
        tempMotif_dict = {}
        counter = 0
        for nuc in range(len(self.seq)): # iterates through each base in seq
            for mot in self.motifs: # for each base, iterate through each motif (O(n*m))
                mLen = self.motifs[mot]
                if (len(self.seq) - counter) > (mLen -1): # prevents indexing past length of seq
                    if re.fullmatch(mot, self.seq[nuc:(nuc+mLen)]):
                        # dict: (start pos, stop pos): [motif1, any perfectly overlapping motifs]
                        if (nuc,(nuc+mLen)) in tempMotif_dict:
                            tempMotif_dict[(nuc,(nuc+mLen))].append(mot)
                        else:
                            tempMotif_dict[(nuc,(nuc+mLen))] = [mot]
            counter += 1
        self.motifs = tempMotif_dict # writing over motif attribute w/new coordinate dict
        return self

    def sequencePurge(self):
        '''
        allows for convenient purging of self.seq to reduce memory usage,
        only name and feature coordinates need to be kept for image drawing
        '''
        self.seq = None
        return self

# Drawing functions

def addRegions(ctx, recs, dims):
    '''
    called by draw function; adds sequence lines, exon denotation, and gene labels
    '''
    width, height, lnSpacing = dims[0], dims[1], dims[2]
    spaceIter = 1
    prev_loc = 0
    for rec in recs: # adding introns (full seq length)
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(0.008)
        ctx.move_to(0.01,lnSpacing * spaceIter)
        ctx.line_to((rec.seqLen / width), lnSpacing * spaceIter)
        ctx.stroke()
        # finding halfway between two genes for the labels (labels above gene)
        label_location = np.mean([(lnSpacing*spaceIter),prev_loc])
        ctx.move_to(0.01,label_location)
        ctx.select_font_face("Calibri")
        ctx.set_font_size(.02)
        ctx.show_text(rec.name) # insert gene name
        prev_loc = lnSpacing * spaceIter
        for exon in rec.exons: # add exons over introns
            ctx.set_line_width(0.05)
            start, stop = exon[0], exon[1] # exon coordinates, to be scaled by image width
            ctx.move_to((start/width + 0.01), lnSpacing * spaceIter)
            ctx.line_to((stop/width + 0.01), lnSpacing * spaceIter)
            ctx.stroke()
        spaceIter += 1
    return ctx


def getColors(motDict):
    '''
    called by draw function; contains hardcoded RGB values for seven hopefully
    colorblind-friendly colors. once exhausted, colors are randomly generated
    '''

    # best approximation of the Darjeeling limited color pallete
    cols = [(.008,.47,.69),
            (.68,.22,0.09),
            (.117,.75,.803),
            (0.82,.61,.18),
            (.29,.6,0.49),
            (.35,.7,.9),
            (0,.6,.5)]
    
    colors = {}

    ln = 0
    for i in motDict:
        if ln < len(cols):
            # use the preset colors until they run out
            colors[motDict[i]] = cols[ln]
        else:
            # randomly generate colors for excessive motifs
            colors[motDict[i]] = np.random.uniform(size = 3)
        ln += 1
    return colors


def addMotifs(ctx, recs, dims, colors):
    '''
    called by draw function; adds size-scaled motif markers to gene graphics
    '''
    width, height, lnSpacing = dims[0], dims[1], dims[2]
    spaceIter = 1
    for rec in recs:
        for key in rec.motifs:
            # new color for each motif
            color = colors[rec.motifs[key][0]]
            ctx.set_source_rgba(color[0], color[1], color[2], 0.7)
            # places motifs on gene same way as exons
            start, stop = key[0], key[1]
            ctx.move_to((start/width + 0.01), lnSpacing * spaceIter)
            ctx.line_to((stop/width + 0.01), lnSpacing * spaceIter)
            ctx.stroke()
        spaceIter += 1
    return ctx


def addLegend(ctx, motDict, colors):
    '''
    called by draw function; adds motif/color legend in bottom right corner of figure
    '''
    ctx.move_to(.777,.56) # Legend title
    ctx.set_font_size(.03)
    ctx.set_source_rgb(0,0,0)
    ctx.select_font_face("Calibri")
    ctx.show_text("Legend")
    
    ctx.set_line_width(0.008) # intron labels
    ctx.move_to(0.78, 0.58)
    ctx.line_to(0.8, 0.58)
    ctx.stroke()
    ctx.set_font_size(.02)
    ctx.move_to(0.82, 0.59)
    ctx.show_text("Introns")
    
    ctx.set_line_width(0.02) # exon labels
    ctx.move_to(0.78, 0.606)
    ctx.line_to(0.8, 0.606)
    ctx.stroke()
    ctx.move_to(0.82, 0.616)
    ctx.show_text("Exons")

    jitter = 0.03 # dictates block spacing
    for i in motDict: # motif labels with block indicating motif color
        # block colors
        color = colors[motDict[i]]
        ctx.set_source_rgb(color[0], color[1], color[2])
        # block size and location
        ctx.set_line_width(0.02)
        ctx.move_to(0.78, 0.606+jitter)
        ctx.line_to(0.8, 0.606+jitter)
        ctx.stroke()
        # label text
        ctx.set_source_rgb(0,0,0)
        ctx.set_font_size(.02)
        ctx.move_to(0.82, 0.616+jitter)
        ctx.show_text(i.upper())
        jitter += 0.03
    return ctx


def draw(recList, name, motDict):
    '''
    takes a list, elements are class gene, output file name, and dict of non-regex:regex motifs. 
    generates motif mark graphic using pycario and calling addRegions, addMotifs, addLegend, getColors functions
    '''
    # order gene objects by seqence length, longest to shortest
    recList.sort(key = lambda rec: rec.seqLen, reverse = True)
    # determine spacing of gene graphics on surface
    lnSpacing = 1 / (len(recList) + 1)
    # set the drawing surface
    width = (recList[0].seqLen * 1.1)
    height = (width*0.75)
    surface = cairo.SVGSurface(name, width, height)
    context = cairo.Context(surface)
    dimensions = [width, height, lnSpacing] # var to pass to context updating functions
    # scaling so everything is 0-1
    context.scale(width, height)
    # background color
    context.rectangle(0,0,width,height)
    context.set_source_rgb(.95,.95,.95)
    context.fill()
    # calling all other drawing functions which return updated context
    context = addRegions(context, recList, dimensions)
    colors = getColors(motDict)
    context = addMotifs(context, recList, dimensions, colors)
    context = addLegend(context, motDict, colors)
    surface.finish() # export svg with fully rendered context

#--------------------------------------------------------------------------------------------------------------
# MAIN
#--------------------------------------------------------------------------------------------------------------

# open input motif file, store motifs in list
input_motifs = []
with open(args.motif, "r") as fh:
    for line in fh:
        input_motifs.append(line.strip("\n"))

# call function to generate regex motif dicts (see expandMotifs)
reg_motifs, saved_motifs = expandMotifs(input_motifs)

# open input fasta file
fasta_fh = open(args.fasta, "r")

# iterate through fasta file 
ln = 0
current_record = [None,""] # initialize list to store records [header,sequence]
record_export = []
for line in fasta_fh:
    if ln == 0: # captures first line
        current_record[0] = line.strip("\n")
    elif re.match("^>", line): # captures all subsequent header lines
        current_record = gene(current_record, reg_motifs)
        current_record = current_record.sequencePurge() # remove seq from memory, no longer needed
        record_export.append(current_record)
        current_record = [None,""] # reset and restart for next record
        current_record[0] = line.strip("\n")
    else: # sequence lines
        current_record[1] += line.strip("\n")
    ln += 1

# handle EOF
current_record = gene(current_record, reg_motifs)
current_record = current_record.sequencePurge()
record_export.append(current_record)

fasta_fh.close()

# output naming and calling draw functions

prefix = re.split("\.", args.fasta)[0]
outName = prefix + ".svg"

# drawing output, automatically exported as svg

draw(record_export, outName, saved_motifs)

