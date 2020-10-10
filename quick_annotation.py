#!/usr/bin/env python
#test commit

import os
import sys
from subprocess import call 
from collections import defaultdict
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
from Bio.Graphics import GenomeDiagram
from BCBio import GFF    
from reportlab.lib.colors import HexColor, red, black, yellow
from config import metagenemarkloc, blastdbloc, pfamloc
###
### Global Variables
###


pks    = ['PF02801', 'PF00109', 'PF08659', 'PF08541', 'PF08545']
nrps   = ['PF08415', 'PF00668', 'PF00975']
pptase = ['PF01648']
acp    = ['']

colors = {  'AirForceBlue'    :  HexColor('#5d8aa8'),    
            'AirForceLight'   :  HexColor('#aac1d1'),    
            'BleuDeFrance'    :  HexColor('#318ce7'),
            'Canary'          :  HexColor('#ffff99'),
            'Fern'            :  HexColor('#71bc78'),
            'LightCadOrange'  :  HexColor('#f3b075'),
            'Int_Orange'      :  HexColor('#ff4f00'),
            'ScreaninGreen'   :  HexColor('#76ff7a'),
            'Wisteria'        :  HexColor('#c9a0dc'),
            'VividAuburn'     :  HexColor('#9f1d35'),
            'LightBrandeis'   :  HexColor('#4194ff'),
            'Red'             :  HexColor('#ff281a'),
            'LightRed'        :  HexColor('#ff4d41') }


def main(fastafile):
    """
    This is the central function that will take a multi-fasta file and then do a few things:
        1. Find ORFS with MetageneMark
        2. Blast the ORFs
        2. Scann HMM on the ORFs
        3. Annotate the Genes
        4. Draw the Contigs
        5. Write out the Blast Tables
        6. Combine Everything into an RST-structured document
        7. use pandoc to make a single html file or a dz-slide model of everything.
    """
    fasta  = "{0}".format(os.path.splitext(fastafile)[0])    
    
    
    run_metagenemark(fastafile, fasta=fasta,length=5000, met_fold   = metagenemarkloc)
    #run_metagenemark(fastafile, fasta=fasta,length=5000, met_fold   = "/home/zcharlop/software/MetaGeneMark_linux64")
    blast(fastafile, fasta=fasta, blastdb = blastdbloc)
    #blast(fastafile, fasta=fasta, blastdb = "/home/gesmlab/NCBI_DATABASES/FOR_BLAST/NR_FOR_BLAST/NR_20130723/nr")
    #blast(fastafile, fasta=fasta, blastdb = "/home/zcharlop/BLAST/NR_20110917/nr")
    hmmscan(fastafile, fasta=fasta, hmmdir=pfamloc )
    #hmmscan(fastafile, fasta=fasta, hmmdir="/home/gesmlab/NCBI_DATABASES/PFAM_A/")
    #get_taxonomy(fastafile, fasta=fasta)
    
    write_gbk(fastafile, fasta=fasta)
    print "annotating..."
    annotate_gbk(fastafile, fasta=fasta)
    #PPTase_GIs_to_file(fastafile)
    #GIs_to_file(fastafile)
    #GI_to_Taxonomy(fastafile)
    #Taxonomy_to_Taxname(fastafile)
    #write_PPTase_positive_hits(fastafile)
    write_linearmap(fastafile)
    #write_circular_maps(fastafile)

###
### Annotation Functions
###

def run_metagenemark(fastafile, fasta, met_fold, length=5000):
    print "running metagenemark..."
    if not os.path.exists("{0}".format(fasta)): os.mkdir("{0}".format(fasta))
    
    #take only sequence >10kb
    length_10kb = []
    longrecords = "{0}/{0}_long.fasta".format(fasta)
    for record in SeqIO.parse(fastafile, "fasta"):
        if len(record.seq) >= length:
            record.description = (record.description).split()[0]
            length_10kb.append(record)
    SeqIO.write(length_10kb, longrecords, "fasta")
    metagene   = "{0}/gmhmmp -a -d  -f G -m {0}/MetaGeneMark_v1.mod -o {1}/{1}.gff {2}".format(met_fold, fasta, longrecords)
    nucleotide = "{0}/nt_from_gff.pl < {1}/{1}.gff > {1}/{1}_nucleotides.fasta".format(met_fold, fasta)
    amino_acid = "{0}/aa_from_gff.pl < {1}/{1}.gff > {1}/{1}_proteins.fasta".format(met_fold, fasta)
    call(metagene, shell=True)
    call(nucleotide, shell=True)    
    call(amino_acid,  shell=True)
    
def blast(fastafile, fasta, blastdb):
    """
    blast a file a file against a database
    """    
    if os.path.exists("{0}/{0}_blast.xml".format(os.path.splitext(fastafile)[0])):
         print "blast has already been run....."
         return
    
    print "running blast"
    print "Blasting....."        
    proteins = "{0}/{0}_proteins.fasta".format(fasta)
    blast    = "blastp -db {0} -query {1} -outfmt 5 -out {2}/{2}_blast.xml -num_threads 20 -max_target_seqs 15".format(blastdb, proteins, fasta)
    print blast
    call(blast,  shell=True)
    print "Blast Complete!"
    print "parsing blast outputs....."
    
def hmmscan(fastafile, fasta, hmmdir):
    """
    hmmscan a file against known HMM profiles.
    """
    
    if os.path.exists("{0}/{0}.hmmtblout".format(os.path.splitext(fastafile)[0])):
        print "hmmscanning has already been run....."
        return  
    print "running hmm"
    proteins = "{0}/{1}_proteins.fasta".format(fasta,fasta)
    hmmscan  = "hmmscan -o {2}/{2}.hmmout --tblout {2}/{2}.hmmtblout  --domtblout {2}/{2}.hmmdomtblout -E .0001 {0}/Pfam-A.hmm {1} ".format(hmmdir,proteins, fasta)
    #print hmmscan
    call(hmmscan,  shell=True)
    print "HMMSCAN Complete!"

def get_taxonomy(genbank_file, blast_results_xml):
    """
    this function will take a genbank file and the XML blast results that sue the CDS from that field as a query.
    It will return a dicitonary of GI numbers
    """
    #GIs = {}
    
    #CDS names = [record.feature.id for record.feature.id] 
    pass

def write_gbk(fastafile, fasta):
    print "parsing blastoutput"
    longrecords = "{0}/{0}_long.fasta".format(fasta)
    gff_file = "{0}/{0}.gff".format(fasta)
    out_file = "{0}/{0}.gbk".format(fasta) 
    fasta_input = SeqIO.to_dict(SeqIO.parse(longrecords, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(gff_iter, out_file, "genbank")

def hmmscan_to_list(fastafile, fasta):
    print "making list of hmmscan output"
    hmmscan  =  open("{0}/{0}.hmmtblout".format(fasta), "r")
    hmm_dict = defaultdict(list)
    hmm_list = []
    for number, line in enumerate(hmmscan):
        if number > 2 and len(line.split())>15:
            #e-val = line.split()[4]
            #gene  = line.split()[2]
            #pfam  = line.split()[1]
            if float(line.split()[4]) < 10e-15:           
                hmm_dict[line.split()[2]].append((line.split()[1]).split(".")[0])
                hmm_list.append((line.split()[1]).split(".")[0])
    #print hmm_dict.values()
    #hmm_hits = set(hmm_list)
    #for x in hmm_hits: print x
    return hmm_dict
   
def generate_colors():
    
    AirForceBlue    =  HexColor('#5d8aa8')
    AirForceLight   =  HexColor('#aac1d1')
    BleuDeFrance    =  HexColor('#318ce7')
    Canary          =  HexColor('#ffff99')
    Fern            =  HexColor('#71bc78')
    LightCadOrange  =  HexColor('#f3b075')
    Int_Orange      =  HexColor('#ff4f00')
    ScreaninGreen   =  HexColor('#76ff7a')
    Wisteria        =  HexColor('#c9a0dc')
    VividAuburn     =  HexColor('#9f1d35')
    LightBrandeis   =  HexColor('#4194ff')
    Red             =  HexColor('#ff281a')
    LightRed        =  HexColor('#ff4d41')

    #important to put this list in reverse order of assigning.
    transferase     = clustercolor( yellow,        ['PF13489', 'PF01209', 'PF07596', 'PF03737', 'PF12847', 'PF08242', 'PF02574', 'PF13847', 'PF08241', 'PF00891' , 'PF13302', 'PF00698', 'PF13469', 'PF02515', 'PF04820', 'PF01039'], "Transferase" )
    sugars          = clustercolor( Fern,          ['PF00908','PF02585','PF02746','PF13579', 'PF00535', 'PF00534', 'PF02397' ,'PF06964', 'PF02719',  'PF01943', 'PF02563', 'PF01270'], "Sugars")
    aa_biosynthesis = clustercolor( AirForceLight, ['PF07587', 'PF01220', 'PF01488', 'PF00857', 'PF05378', 'PF01968', 'PF02423' ,'PF03301', 'PF01425', 'PF00848' , 'PF00733'], "Amino Acid Biosynthesis")
    red_ox          = clustercolor( VividAuburn,   ['PF00420', 'PF01512', 'PF01257', 'PF01370', 'PF00146', 'PF00499','PF00361', 'PF00662', 'PF00106' ,'PF01613' ,'PF00141' ,'PF02578'], "Redox Proteins")
    membrane        = clustercolor( Wisteria,      ['PF07516', 'PF07517', 'PF01043', 'PF00902','PF00324', 'PF13520', 'PF01312','PF00664',  'PF01061', 'PF00005','PF03741', 'PF00691', 'PF07158', 'PF07947','PF02133', 'PF00528', 'PF01580', 'PF01514', 'PF07715', 'PF02472', 'PF01311', 'PF01313', 'PF00771', 'PF00999', 'PF00939', 'PF01618','PF00006'], "Membrane")
    pptase          = clustercolor( AirForceBlue,  ['PF01648'], "PPTase")
    transcription   = clustercolor( LightCadOrange, ['PF03704', 'PF13453', 'PF01614', 'PF03466', 'PF02482', 'PF02909', 'PF06613', 'PF00072', 'PF00196', 'PF08535','PF00480'], "Transcription")
    fas             = clustercolor( ScreaninGreen, ['PF08541', 'PF08545'], "fas")
    pks             = clustercolor( red  ,         ['PF02801', 'PF00109', 'PF08659'], "pks")
    nrps            = clustercolor( LightRed ,     ['PF08415', 'PF00668', 'PF03621' ,'PF14399', 'PF00550', 'PF00975', 'PF00501'], "nrps")
    annotations     = [ transferase, pptase, transcription, aa_biosynthesis, red_ox, sugars, membrane, fas, pks,nrps]
    return annotations

def annotate_gbk(fastafile, fasta):
    #setup files and dicts
    print "annotating genbank file"
    gbkfile  = "{0}/{0}.gbk".format(fasta)
    proteins = "{0}/{0}_proteins.fasta".format(fasta)
    protein_dict = SeqIO.to_dict(SeqIO.parse(proteins, "fasta", generic_dna))
    blastxml = "{0}/{0}_blast.xml".format(fasta)
    blastdict = make_blast_dictionary(blastxml)
    hmm_dict = hmmscan_to_list(fastafile, fasta) 
    output_handle = open("{0}/{0}_mod.gbk".format(fasta), "w")
    sequences = []
    
     #looping
    for record in SeqIO.parse(open(gbkfile, "r"), "gb"):
        for feature in record.features :
            if feature.type=="CDS" :
                feature.qualifiers['domains'] = (hmm_dict["gene_id_{0}".format(feature.qualifiers['gene_id'][0])])
                if feature.qualifiers['gene_id'][0] in blastdict.keys():
                    #print "Blastdict of Gene ID  {0}".format(blastdict[feature.qualifiers['gene_id'][0] ])
                    feature.qualifiers['product'] = blastdict[feature.qualifiers['gene_id'][0] ]
                    #print "Protein Dict of Gene_ID:  {0}".format(protein_dict[feature.qualifiers['gene_id'][0] ] )
                    feature.qualifiers['translation'] = (protein_dict["gene_id_{0}".format(feature.qualifiers['gene_id'][0])]).seq
            else: break
        sequences.append(record)
    SeqIO.write(sequences, output_handle, "genbank")
    output_handle.close()

def write_PPTase_positive_hits(fastafile):
    """
    loops through the gbk file and looks for proteins that have the pptase_pfam signature
    """
    def intersect(a, b):
        """ returns the number of common members in two lists """
        return len(list(set(a) & set(b)))
    print "writing out PPTase_positive_GBK"
    fasta    = "{0}".format(os.path.splitext(fastafile)[0])
    gbkfile  = "{0}/{0}_mod.gbk".format(fasta)
    output_handle  = open("{0}_PPTases_proteins.faa".format(fasta), "w")
    PPTase_clones = []
    pptase_pfam = ['PF01648']
    output_handle = open("{0}/{0}_PPTases.gbk".format(fasta), "w")    
    for record in SeqIO.parse(gbkfile, "gb"):
        #print len(record.seq)
        if 10000 < len(record.seq) < 55000: #set the size cutoff here
            for feature in record.features:
                    if feature.type == "CDS":  #change features and attributes here.
                        if 'domains' in feature.qualifiers.keys():
                            if intersect( feature.qualifiers['domains'], pptase_pfam ) > 0:
                                #print "PPTase!"
                                PPTase_clones.append(record)
                                if 'product' in feature.qualifiers.keys():
                                    output_handle.write(">{0}_{1} \n{2} \n".format( record.name, str(feature.qualifiers['product'][0])[2:2],feature.qualifiers['translation'][0]))
                                
                                
    PPTases = set(PPTase_clones)
    #for record in PPTases:
    #    record.name = "{1} {0}".format(record.name.split("contig")[1], fasta.split("_")[0])
    SeqIO.write(PPTases, output_handle, "genbank")

def PPTase_GIs_to_file(fastafile):
    """ 
    this will make alist of PPTase positive clones from the annotated file, make a list of proteins in from these clones.
    It will make a dicitonary of the the top blast hits and then will query the PPTase proteins against that list to get
    the GI nubmers and it will write the GIs to a file.
    """
    print "getting PPTase positive genes and writing GI's of top blast hits....."
    fasta        = "{0}".format(os.path.splitext(fastafile)[0])
    gbkfile      = "{0}/{0}_mod.gbk".format(fasta)
    blastxmlfile = "{0}/{0}_blast.xml".format(fasta)
    output_handle = open("{0}/Gi_list_{0}.txt".format(fasta), "w")
   # output_handle_pptase_GIs = open("{0}/{0}_PPTASE_genes.faa".format(fasta),"w")
    
    blastdict         = {}
    #blastdict_pptases = {}
    
    def intersect(a, b):
        """ returns the number of common members in two lists """
        return len(list(set(a) & set(b)))
    
    #make list of proteins in cosmids with PPTases
    PPTase_clones = []
    #pptases       = []
    pptase_genes  = []
    pptase_pfam   = ['PF01648']
    
    #loop to get a list of genes from PPTase Positive clones
    for record in SeqIO.parse(gbkfile, "gb"):
        if 10000 < len(record.seq) < 55000: #set the size cutoff here
            for feature in record.features:
                    if feature.type == "CDS":  
                        if 'domains' in feature.qualifiers.keys():
                            if intersect( feature.qualifiers['domains'], pptase_pfam ) > 0:
                                          
                                PPTase_clones.append(record)
    PPTase_clones = set(PPTase_clones)
    for record in PPTase_clones:
        for feature in record.features:
                if feature.type == "CDS":
                    pptase_genes.append(feature.qualifiers['gene_id'][0])
    #print pptase_genes
    
    #Loop to make a dictionary of Blast results with GIs.
    for record in NCBIXML.parse(open(blastxmlfile)):
        for number, align in enumerate(record.alignments):
            for x, hsp in enumerate(align.hsps):
                if x == 2: break
                #print "Number: {0}".format(x)
                if hsp.expect < 10e-25:
                    blastdict[str(record.query).split("_")[2]] =  str((align.title).split(">")[0]).split("|")[1] # this returns the GI
    #print blastdict.keys()
    #cross reference the gene list with the blast dictionary
    for pptase_gene in pptase_genes:
        if pptase_gene in blastdict.keys():
            output_handle.write("{0}\n".format(blastdict[pptase_gene]))
    output_handle.close()

def GIs_to_file(fastafile):
    """ 
    clone of PPTase_GIs only will write all ORFs for 50 cosmids 
    """
    print "getting PPTase positive genes and writing GI's of top blast hits....."
    fasta        = "{0}".format(os.path.splitext(fastafile)[0])
    gbkfile      = "{0}/{0}_mod.gbk".format(fasta)
    blastxmlfile = "{0}/{0}_blast.xml".format(fasta)
    output_handle = open("{0}/Gi_list_{0}.txt".format(fasta), "w")
    blastdict = {}
    
    def intersect(a, b):
        """ returns the number of common members in two lists """
        return len(list(set(a) & set(b)))
    
    #make list of proteins in cosmids with PPTases
    PPTase_clones = []
    pptase_genes  = []
    #pptase_pfam   = ['PF01648']
    
    #loop to get a list of genes from PPTase Positive clones
    x = 1
    for record in SeqIO.parse(gbkfile, "gb"):
        x  += 1
        if x == 50: break
        if 10000 < len(record.seq) < 55000: #set the size cutoff here
            for feature in record.features:
                    if feature.type == "CDS":  
                        if 'domains' in feature.qualifiers.keys():
                            PPTase_clones.append(record)
    PPTase_clones = set(PPTase_clones)
    for record in PPTase_clones:
        for feature in record.features:
                if feature.type == "CDS":
                    pptase_genes.append(feature.qualifiers['gene_id'][0])
    #print pptase_genes
    
    #Loop to make a dictionary of Blast results with GIs.
    for record in NCBIXML.parse(open(blastxmlfile)):
        #temp=[]     
        for number, align in enumerate(record.alignments):
            for x, hsp in enumerate(align.hsps):
                if x == 2: break
                #print "Number: {0}".format(x)
                if hsp.expect < 10e-25:
                    blastdict[str(record.query).split("_")[2]] =  str((align.title).split(">")[0]).split("|")[1] # this returns the GI
    #print blastdict.keys()
    #cross reference the gene list with the blast dictionary
    for pptase_gene in pptase_genes:
        if pptase_gene in blastdict.keys():
            output_handle.write("{0}\n".format(blastdict[pptase_gene]))
    output_handle.close()


def GI_to_Taxonomy(fastafile):
    print "generating Taxonomy information from GI's of top blast hits....."
    fasta         = "{0}".format(os.path.splitext(fastafile)[0])
    gi_file       = open("{0}/Gi_list_{0}.txt".format(fasta))
    gi_list       = [line.rstrip() for line in gi_file]
    gi_file.close()
    
    #note this is the protein file
    taxonomy_file = open("/home/gesmlab/NCBI_DATABASES/TAXONOMY/gi_taxid_prot.dmp", "r")
    output_handle = open("{0}/{0}_PPTaseclone_taxonomylist.txt".format(fasta), "w")
    
    
    #print "GI LIst = {0}".format(gi_list)
    #print gi_list
    
    for line in taxonomy_file: #gi is the first column, tax ID the second
        if line.split()[0] in gi_list: output_handle.write("{0}\n".format(line.split()[1]))
        #print line.split()[0]
        #print line.split()[1]
    output_handle.close()
    taxonomy_file.close()

def Taxlist_To_Newick(fastafile):
    print "generating Newick tree"
    fasta         = "{0}".format(os.path.splitext(fastafile)[0])
    taxdir     = "/home/zcharlop/PPTASE_Abundance/jhcepas-ncbi_taxonomy-7e29f37"
    currentdir = os.getcwd()
    command    = "python2.6 ncbi_query.py -tf  {1}/{0}/{0} -x > {1}/{0}/{0}.pre_nw".format(fasta, currentdir, taxdir)
    #python2.6 ncbi_query.py -tf ALL_PPTASE_Taxonomy.txt -x > ALL_PPTASE.pre_nw




def Taxonomy_to_Taxname(fastafile):
    print "generating Taxonomy Names from TaxIDs....."
    fasta         = "{0}".format(os.path.splitext(fastafile)[0])
    taxidfile       = open("{0}/{0}_PPTaseclone_taxonomylist.txt".format(fasta))
    taxidlist       = [line.rstrip() for line in taxidfile]
    taxidfile.close()
    
    #note this is the protein file
    taxonomy_file = open("/home/gesmlab/NCBI_DATABASES/TAXONOMY/names.dmp", "r")
    output_handle = open("{0}/{0}_PPTaseclone_taxonomy_names.txt".format(fasta), "w")
    
    for line in taxonomy_file: #taxid is the first column, name the second, type the third
        if 'scientific' in line.split("|")[3]:
            #print line.split("|")[0]  
            if (line.split("|")[0]).strip() in taxidlist:
                output_handle.write("{0}\n".format(line.split("|")[1]))
                print line.split("|")[0], line.split("|")[1], line.split("|")[3]

        #print line.split()[0]
        #print line.split()[1]
    output_handle.close()
    taxonomy_file.close()

def get_GI_dict(blastxml):
    """
    retun a dict of gene:GI for top hsp of a blast query
    """
    blastdict = {}
    for record in NCBIXML.parse(open(blastxml)):
        for number, align in enumerate(record.alignments):
            for x, hsp in enumerate(align.hsps):
                if x == 2: break
                #print "Number: {0}".format(x)
                #if hsp.expect < 10e-25:
                blastdict[record.query] =  str(align).split("|")[1]   
    #print blastdict
    return blastdict

def get_gene_dict(gbk):
    """
    parse a record and return a set of genes
    """
    genedict ={}
    for record in SeqIO.parse(open(gbk,'r'),'genbank'):
        #print record.id
        cds = [feature.qualifiers['gene_id'][0] for feature in record.features if feature.type=='CDS']
        genedict[record.id] = cds
    return genedict

def record_GI_dict(gene_dict,blastdict):
    """
    take two dicitonaries - one of the Record:Genes
    and one with Gene:GI, and return Record:GI
    """
    GI_dict = defaultdict(list)
    for record,genes in gene_dict.itervalues():
        if len(genes) >0: 
            for gene in genes:
                GI_dict[record].append(blastdict[gene])
                
    return GI_dict
    
def make_blast_dictionary(blastxml):
    """ 
    this will take a blast output xml file and return the top hit for each
    query in the form of a dicitonary.
    blastxml is the file.
    """
    blastdict = {}
    for record in NCBIXML.parse(open(blastxml)):
        for number, align in enumerate(record.alignments):
            for x, hsp in enumerate(align.hsps):
                if x == 2: break
                #print "Number: {0}".format(x)
                #if hsp.expect < 10e-25:
                blastdict[str(record.query).split("_")[2]] =  str((align.title).split(">")[0]).split("|")[4]   
    #print blastdict
    return blastdict

###
### Display-Related Functions
###

def make_blast_dictionary_for_output(blastxml):
    """ 
    this will take a blast output xml file and return the top hit for each
    query in the form of a dicitonary.
    blastxml is the file.
    """
    blastdict = {}
    for record in NCBIXML.parse(open(blastxml)):
        #print record.query
        temp=[]     
        for number, align in enumerate(record.alignments):
            temp2=[]
            temp2.append(str((align.title).split("|")[3]))
            temp2.append((str((align.title).split(">")[0]).split("|")[4]).split("[")[0][:30])
            temp.append(temp2)
            if number == 10: break
        blastdict[str(record.query).split("_")[2]] = temp
            #print align
            #print (align.title).split("|")[3]
            #print (str((align.title).split(">")[0]).split("|")[4]).split("[")[0]
            #print (str((align.title).split(">")[0]).split("|")[4]).split("[")[0][:30]   
            #for x, hsp in enumerate(align.hsps):
             #   if x == 10: break
                #print "Number: {0}".format(x)
                #if hsp.expect < 10e-25:
                #print hsp
                #blastdict[str(record.query).split("_")[2]] =  str((align.title).split(">")[0]).split("|")[4]   
    #print blastdict
    return blastdict

def write_linearmap(fastafile):
    
    def intersect(a, b):
        """ returns the number of common members in two lists """
        return len(list(set(a) & set(b)))
    
    print "write linear map"
    fasta    = "{0}".format(os.path.splitext(fastafile)[0])
    gbkfile  = "{0}/{0}_mod.gbk".format(fasta)
    picdir   = "{0}/pics".format(fasta)
    montage_list = []
    gb_list = []
    #get colors and inverted colors for assignment
    colors = generate_colors()
    inverted_colors = {}
    for genetype in colors:
        for prop in genetype.properties:
            inverted_colors[prop] = genetype.color
    for record in SeqIO.parse(gbkfile, "gb"):
        gd_diagram = GenomeDiagram.Diagram(record.name, tracklines=False, track_size=.5) #Define genome Diagram
        gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features", scale=0) #define track using the file name
        gd_feature_set = gd_track_for_features.new_set() #define featureset         
        for feature in record.features:
            if feature.type == "CDS":  #change features and attributes here.
                #print feature.qualifiers
                if not 'domains' in feature.qualifiers.keys():
                    gd_feature_set.add_feature(feature, border=black, sigil="BIGARROW",color="white", label=True,arrowshaft_height=1,
                                                label_size = 0, label_angle=0)
                elif intersect( feature.qualifiers['domains'],inverted_colors.keys() ) == 0:
                    gd_feature_set.add_feature(feature, border=black, sigil="BIGARROW",color="white", label=True,arrowshaft_height=1,
                                               label_size = 0, label_angle=0)
                else:
                    for genetype in colors:
                        if intersect( feature.qualifiers['domains'],genetype.properties ) > 0:
                            gd_feature_set.add_feature(feature, color = genetype.color, border=black, sigil="BIGARROW", label=True,arrowshaft_height=1,
                                                       label_size = 0, label_angle=0)
        gd_diagram.draw(format="linear", orientation="landscape", pagesize=(len(record.seq)/20,100), x=0.01, y=0.01, fragments=1, start=0, end=len(record))
        if not os.path.exists(picdir): os.mkdir(picdir)
        imgname = "{0}/{1}.png".format(picdir,record.name)
        gd_diagram.write(imgname , "PNG")
        imgname_pdf = "{0}/{1}.pdf".format(picdir,record.name)
        gd_diagram.write(imgname_pdf , "PDF")
        #append_value = " -label {0}  {1}imgname ".format(record.name,imgname)
        montage_list.append( " -append" " " + imgname + " ")
        gb_list.append(imgname)
        #note that the purpose of tHe convert command is to reduce whitespace. this will need to be changes if pagesize changes
        convert =  " convert " + imgname + " -shave 0x25 "  + imgname
        call(convert.split())
    ## make the montage
    print "Contigs who images are written: {0}".format(len(montage_list))
    #command =  "convert {0} {1}/{2}.png".format(''.join(montage_list), picdir, fasta)
    #call(command.split())



def write_circular_maps(gbk_file):
    def intersect(a, b):
        """ returns the number of common members in two lists """
        return len(list(set(a) & set(b)))
    
    print "write circular map"
    gbk    = "{0}".format(os.path.splitext(gbk_file)[0])
    gbkfile  = "{0}/{0}_mod.gbk".format(gbk)
    if not os.path.exists("{0}".format(gbk)): os.mkdir("{0}".format(gbk)) 
    picdir   = "{0}/circularpics".format(gbk)
    montage_list = []
    montage_notitle = []
    gb_list = []
    #get colors and inverted colors for assignment
    colors = generate_colors()
    inverted_colors = {}
    for genetype in colors:
        for prop in genetype.properties:
            inverted_colors[prop] = genetype.color
    for record in SeqIO.parse(open(gbkfile), "gb"):
        gd_diagram = GenomeDiagram.Diagram(record.name, tracklines=False, track_size=.2) #Define genome Diagram
        gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features", scale=0) #define track using the file name
        gd_feature_set = gd_track_for_features.new_set() #define featureset         
        for feature in record.features:
            if feature.type == "CDS":  #change features and attributes here.
                #print feature.qualifiers
                if not 'domains' in feature.qualifiers.keys():
                    gd_feature_set.add_feature(feature, border=black, sigil="BIGARROW",color="white", label=True,
                    shaft_height=1,
                                                label_size = 0, label_angle=0, strand=None)
                elif intersect( feature.qualifiers['domains'], pptase ) > 0:
                    gd_feature_set.add_feature(feature, border=black, sigil="BIGARROW",color= yellow, label=True,arrowshaft_height=1,
                                               label_size = 0, label_angle=0, strand=None)
                elif intersect( feature.qualifiers['domains'], nrps ) > 0:
                    gd_feature_set.add_feature(feature, border=black, sigil="BIGARROW",color= red, label=True,arrowshaft_height=1,
                                               label_size = 0, label_angle=0, strand=None)
                elif intersect( feature.qualifiers['domains'], pks ) > 0:
                    gd_feature_set.add_feature(feature, border=black, sigil="BIGARROW",color= red, label=True,arrowshaft_height=1,
                                               label_size = 0, label_angle=0, strand=None)
                else:
                    gd_feature_set.add_feature(feature, border=black, sigil="BIGARROW",color="white", label=True,arrowshaft_height=1,
                                                label_size = 0, label_angle=0, strand=None)

        gd_diagram.draw(format="circular", orientation="landscape", pagesize=(400,400), x=0.01, y=0.01, fragments=1, start=0, end=len(record))
        if not os.path.exists(picdir): os.mkdir(picdir)
        imgname = "{0}/{1}.png".format(picdir,record.name)
        print "Image Name = {0}".format(imgname)
        gd_diagram.write(imgname , "PNG")
        imgname_pdf = "{0}/{1}.pdf".format(picdir,record.name)
        gd_diagram.write(imgname_pdf , "PDF")
        montage_list.append( " -label " + "{0}_{1}kb".format(record.name, str(len(record.seq))[:2]) + " " + imgname + " ")
        montage_notitle.append(" " + imgname + " ")
        gb_list.append(imgname)
        #note that the purpose of tHe convert command is to reduce whitespace. this will need to be changes if pagesize changes
        convert =  " convert {0} -gravity center -extent 250x250 {0} ".format(imgname)
        call(convert.split())
        
    ## make the montage
    print "Circular contigs who images are written: {0}".format(len(montage_list))
    #command =  "montage {0} {1}/{2}.png".format(''.join(montage_list), picdir, gbk)
    command =  "montage -font Myriad-Pro-Regular -title {2} {0} {1}/{2}.png".format(''.join(montage_list), picdir, gbk)
    call(command.split())
    command =  "montage -font Myriad-Pro-Regular -title {2} {0} {1}/{2}_notitles.png".format(''.join(montage_notitle), picdir, gbk)
    call(command.split())
###
### Classes
###

def make_html(fastafile):
    def write_jinja_template():
        fasta  = "{0}".format(os.path.splitext(fastafile)[0])
        f = open("{0}/{0}.rst".format(fasta),"w")
        f.write(
            """
            
            """
        )
        f.close()
            
    
    fasta  = "{0}".format(os.path.splitext(fastafile)[0])
    f = open("{0}/{0}.rst".format(fasta),"w")
    f.write(
        """
        
        """
    )
    f.close()
    
#def write_documentation(fastafile):
#    test = 'a'   

###
### Classes
###
class clustercolor(object):
        def __init__(self, color, properties, name):
            self.color = color
            self.properties = properties
            self.name = name
    
if __name__ == "__main__":
    main(sys.argv[1])
