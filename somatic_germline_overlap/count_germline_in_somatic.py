#!/bin/python
#03 February 2016 - Kuan-Lin Huang @ WashU - 
 
import sys
import getopt
import gzip

class autovivification(dict):
    '''Implementation of perl's autovivification feature.'''
    def __init__( self , *args , **kwargs ):
        super( autovivification , self ).__init__( *args , **kwargs )
        self.itemlist = super( autovivification , self ).keys()
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def main():
    def usage():
        print """
    count_germline_in_somatic.py : why do I exist?

    USAGE: count_germline_in_somatic.py [-h] <filename>
     -h    print this message
     <filename>    input file
        """

    genomic_pos_count = autovivification()
    genomic_exact_change_count = autovivification()
    peptide_pos_count = autovivification()
    peptide_exact_change_count = autovivification()
    peptide_pos_start = autovivification()

    #use getopt to get inputs
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h') #:after option meaning required arguments
    except getopt.GetoptError:
        print "count_germline_in_somatic.py  <germline file> <somatic file>"

    for opt, arg in opts: #store the input options
        if opt == '-h': # h means user needs help
            usage(); sys.exit()

   

    try:
        germlineFH= "../../TCGA_data/germline/PCA_pathVar_integrated_filtered_adjusted.tsv"
        germlineF = open(germlineFH,"r")
    except IOError:
        print("Germline file does not exist!")

    header = germlineF.readline().strip()
    # initiate dicts for germline variants
    for line in germlineF:
        line=line.strip()
        F = line.split("\t")
        i = [1,8,4,9,10]
        varList = [F[j] for j in i]
        var = "_".join(varList)
        alt = F[11]
        genomic_pos_count[var] = 0
        genomic_exact_change_count[var][alt] = 0

        k= [1,8,115]
        pepVarList = [F[j] for j in k]
        pepVar = "_".join(pepVarList)
        pepVarPos = pepVar[:-1]
        start = F[4]

        peptide_pos_count[pepVarPos] = 0
        peptide_pos_start[pepVarPos] = start
        peptide_exact_change_count[pepVar] = 0

    germlineF.close()

    #open input file
    try:
        somaticFH= "../../TCGA_data/somatic/mc3.v0.2.8.PUBLIC.maf.gz"
        somaticF = gzip.open(somaticFH,"r")
    except IOError:
        print("Somatic file does not exist!")

    header = somaticF.readline().strip()
    # accumulates
    for line in somaticF:
        line=line.strip()
        F = line.split("\t")

        # case for same genomic variant
        i = [0,4,5,6,10]
        varList = [F[j] for j in i]
        var = "_".join(varList)
        alt = F[12]

        if var in genomic_pos_count:
            genomic_pos_count[var] += 1

        if var in genomic_exact_change_count:
            if alt in genomic_exact_change_count[var]:
                genomic_exact_change_count[var][alt] += 1

        k= [0,4,36]
        pepVarList = [F[j] for j in k]
        pepVar = "_".join(pepVarList)
        pepVarPos = pepVar[:-1]
        start = F[4]
        pos = int(F[5])

        if pepVarPos in peptide_pos_count:
            start = int(peptide_pos_start[pepVarPos])
            if pos > (start - 2) and pos < (start + 2):
                peptide_pos_count[pepVarPos] +=1

        if pepVar in peptide_exact_change_count:
            start = int(peptide_pos_start[pepVarPos])
            if pos > (start - 2) and pos < (start + 2):
                peptide_exact_change_count[pepVar] +=1

    somaticF.close()

    outFH = "../../TCGA_data/germline/PCA_pathVar_integrated_filtered_adjusted_wSomaticCounts.tsv"
    sys.stdout=open(outFH,"w")

    
    try:
        germlineFH= "../../TCGA_data/germline/PCA_pathVar_integrated_filtered_adjusted.tsv"
        germlineF = open(germlineFH,"r")
    except IOError:
        print("Germline file does not exist!")

    header = germlineF.readline().strip()
    headerF = header.split("\t")
    print "\t".join(headerF) + "\tSomaticMutation_genomic_location_count\tSomaticMutation_genomic_change_count\tSomaticMutation_peptide_location_count\tSomaticMutation_peptide_change_count" 

    for line in germlineF:
        line=line.strip()
        F = line.split("\t")
        i = [1,8,4,9,10]
        varList = [F[j] for j in i]
        var = "_".join(varList)
        alt = F[11]

        k= [1,8,115]
        pepVarList = [F[j] for j in k]
        pepVar = "_".join(pepVarList)
        pepVarPos = pepVar[:-1]
        start = F[4]


        GenomicPosCount = genomic_pos_count[var]
        GenomicChangeCount = genomic_exact_change_count[var][alt]
        PeptidePosCount = peptide_pos_count[pepVarPos]
        PeptideChangeCount = peptide_exact_change_count[pepVar]

        
        print "\t".join(F) + "\t"  + str(GenomicPosCount)+ "\t"  + str(GenomicChangeCount)+ "\t" + str(PeptidePosCount) + "\t" + str(PeptideChangeCount)

    germlineF.close()
    sys.stdout.close()

if __name__ == "__main__":
    main()
