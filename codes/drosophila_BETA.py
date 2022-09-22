#!/bin/python

#########################
####    Peiyu Shi     ###
####    2021-09-16    ###
#########################

"""Script Description:

    This script is modified version of BETA-Binding and Expression Targets Analysis (BETA),which is suitable for analyzing datastes from Drosophila melanogaster (fly).It integrates ChIP-seq/CUT&Tag binding data and RNA-seq differentail expression data to predict a factor's targets.

    Part1:
    use binding data to calculate the regulatory potential score of each gene to be regulated by factor.
    1. For each gene in genome, input a distance(d). All binding sites near the transcription start site of the gene within the specified range (10 kb as default) are considered.
    2. calculate a sum of regulatory potential score for each gene use this formula: 
    Sg = lambda ldx: sum([math.exp(-0.5-4*t) for t in ldx])
    3. output is in bed format. the 5th column is score.
"""

import sys, os, time
import math
import re
from subprocess import call as subpcall
from optparse import OptionParser

# specify analyzed chromosomes
chroms = ["chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY"]

# regulatory potential score calculation function
Sg = lambda ldx: sum([math.exp(-0.5-4*t) for t in ldx])

# print current time and information on screen
def Info(infoStr):
        print "[%s] %s" %(time.strftime('%H:%M:%S'), infoStr)

def prepare_optparser():
        usage = """usage: %prog <-p binding_file> <-n name> <-g genome> <-d distance> [options]
            example : python drosophila_BETA.py -p factor1.narrowpeak -n factor1 -g /public/home/shipy3/DB/dm6/annotation/BETA_input.txt"""
        description = "BETA for drosophila --- Binding Expression Target Analysis for drosophila"
        # create a OptionParser
        optparser = OptionParser(version="%prog v1.00", description=description, usage=usage, add_help_option=False)
        # add options
        optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
        optparser.add_option("-p","--peakfile",dest="peakfile",type="string",help="Input the bed/narrowpeak format peak file of the factor")
        optparser.add_option("-n","--name",dest="name",type="string",
                                         help="this argument is used to name the result file")
        optparser.add_option("-d","--distance", dest="distance", type="int",
                                         help="Set a number which unit is 'base'. It will get peaks within this distance from gene TSS. default:10000 (10kb)", default=10000)
        optparser.add_option("-g","--genome",dest="genome",type="string",
                                         help="Select an annotaion file with the UCSC bed format file to search genes.")
        (options,args) = optparser.parse_args()

        if not options.peakfile and not options.genome:
            optparser.print_help()
            sys.exit(1)
        if not os.path.isfile(options.peakfile):
            Info('ERROR: Cannot find peak file, a tab-peak file must be given through -p(--peakfile).')
            sys.exit(1)
        if not os.path.isfile(options.genome):
            Info("ERROR: Genome file not found! An annottion file with the UCSC bed format file must be given through -g (--genome).")
            sys.exit(1)
        if not options.name:
            Info("ERROR: Name of the result file is not specificed. You can specify it by -n (--name).")
            sys.exit(1)

        Info("Argument List: ")
        Info("Name = " + options.name)
        Info("Peak File = " + options.peakfile)
        Info("Distance = %d bp" %options.distance)
        Info("Genome = %s" %options.genome)

        return options


##########################
##        PartI         ##
##########################

class PScore:
    #Caculate every gene's regulatory potential score, peaks within 10kb will be cosidered.
    def __init__(self, options):
        self.peakfile = options.peakfile
        self.genome = options.genome

        self.opts_string = "# Argument List:\n" +\
                           "# Name = %s\n" %options.name +\
                           "# peak file = %s\n" %options.peakfile +\
                           "# distance = %d bp\n" %options.distance 

        self.peaklist = {}

    def readfile(self):#reads the file and returns a vector: each element is a bed_row. 
        peakf = open(self.peakfile)
        peakf = peakf.readlines()
        peaks = []

        for peak in peakf:
            if peak.startswith("#") or not peak.strip():
                continue
            peak = peak.strip()
            peak = peak.split('\t')
            peaks.append(peak)

        count = 0
        self.peaklist = {}

        for line in peaks:
            #.bed-> 0:chrom 1:pStart 2:pEnd
            line = [line[0], int(line[1]), int(line[2])]
            try:
                self.peaklist[line[0]].append(line)
            except KeyError:
                self.peaklist[line[0]] = [line]
            count += 1
        
        for i in self.peaklist.keys():
            self.peaklist[i].sort()
        Info("Read file <%s> OK! All <%d> peaks." %(self.peakfile, count))


    def ScoreCalc(self, distance):
        #calculate each gene's regulatory potential sg = ...
        gene = open(self.genome,'r')
        self.geneInfo = []
        for line in gene:
            if line.startswith('#') or not line.strip():
                continue
            else:
                line = line.strip()
                line = line.split('\t')
                if line[1] in chroms:
                    info = [line[1], line[0], str(line[3]), str(line[4]), line[2], line[5]]
                    self.geneInfo.append(info)
                else :
                    continue

        count = 0
        for igene in self.geneInfo:
            if igene[4]=="+" :
                gTSS = int(igene[2])
            else :
                gTSS = int(igene[3])
            try:
                peaks = self.peaklist[igene[0]]
            except KeyError:
                peaks = []

            peaksInDistance = []
            for t in peaks :
                if gTSS >= t[1] and gTSS <= t[2] :
                    peaksInDistance.append(0*1.0/distance)
                elif abs(t[1]-gTSS) < distance or abs(t[2]-gTSS) < distance :
                    peaksInDistance.append(min(abs(t[1]-gTSS),abs(t[2]-gTSS))*1.0/distance)
            #peaksInDistance = [min(abs(t[1]-gTSS),abs(t[2]-gTSS))*1.0/distance for t in peaks if abs(t[1]-gTSS) < distance or abs(t[2]-gTSS) < distance ]
            peaksInDistance.sort()
            igene.append(Sg(peaksInDistance))
            count += 1

        Info('Process <%d> genes'%count)
        self.geneInfo.sort(key=lambda x:x[-1], reverse=True)

    def Output2File(self, name,distance):
        outf = open("{}_{}bp_gene2score.txt".format(name,distance),"w")#peaks score and rank file
        outf.write('#chrom\ttxStart\ttxEnd\tgene\tscore\tstrand\tsymbol\trank\n')
        r=1
        for line in self.geneInfo:
            if str('%.3f'%line[6]) == '0.000':
                #if one gene's score is zero, this gene will not be used to rank
                outf.write('%s\t%d\t%d\t%s\t%.3f\t%s\t%s\t%s\n'%(line[0], int(line[2]), int(line[3]), line[1], line[6], line[4], line[5], 'NA'))
            else :
                outf.write('%s\t%d\t%d\t%s\t%.3f\t%s\t%s\t%d\n'%(line[0], int(line[2]), int(line[3]), line[1], line[6], line[4], line[5], r))
            r += 1
        outf.close()
        Info("Finished! Gene2Score result output to <{}_{}bp_gene2score.txt>".format(name,distance))

    def OutputPeak2Gene(self, name,distance):
        outf = open("{}_{}bp_peak2gene.txt".format(name,distance), "w")
        outf.write('#chrom\tpStart\tpEnd\tgene\tsymbol\tdistance\tscore\n')
        for igene in self.geneInfo:
            if igene[4]=="+" :
                gTSS = int(igene[2])
            else :
                gTSS = int(igene[3])
            try:
                peaks = self.peaklist[igene[0]]
            except KeyError:
                peaks = []

            for t in peaks:
                if gTSS >= t[1] and gTSS <= t[2]:
                    peak2gene_distance = 0
                    score = math.exp(-0.5)
                    outf.write('%s\t%d\t%d\t%s\t%s\t%d\t%.3f\n'%(t[0], int(t[1]), int(t[2]),igene[1],igene[5], int(peak2gene_distance),float(score)))
                elif abs(t[1]-gTSS) < distance or abs(t[2]-gTSS) < distance :
                    if abs(t[1]-gTSS) <= abs(t[2]-gTSS):
                        peak2gene_distance = t[1]-gTSS
                    else :
                        peak2gene_distance = t[2]-gTSS
                    normalized_distance = abs(peak2gene_distance)*1.0/distance
                    score = math.exp(-0.5-4.0*normalized_distance)
                    outf.write('%s\t%d\t%d\t%s\t%s\t%d\t%.3f\n'%(t[0], int(t[1]), int(t[2]),igene[1],igene[5], int(peak2gene_distance),float(score)))
        outf.close()
        Info("Finished! Peak2Gene result output to <{}_{}bp_peak2gene.txt>".format(name,distance))


def main():
    start = time.time()
    opts=prepare_optparser()
    g = PScore(opts)
    g.readfile()
    g.ScoreCalc(opts.distance)
    g.Output2File(opts.name,opts.distance)
    g.OutputPeak2Gene(opts.name,opts.distance)
    end = time.time()
    total = end - start
    hour = int(total/3600)
    minite = int(total - hour*3600)/60
    second = int(total - hour*3600 - minite*60)
    print 'total time: %s:%s:%s '%(hour, minite, second)

if __name__ == "__main__":
    main()
    # python drosophila_BETA.py -p ~/ChinmoST/output/chinmo_cut_tag/MACS2/WT_D57_rep1_filtered_peaks.narrowPeak -n WT_D57_rep1 -g ~/DB/dm6/annotation/BETA_input.txt
