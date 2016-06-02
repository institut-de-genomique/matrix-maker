#!/usr/bin/env python

################################################################################
##                      
##                             Rscript
##                            
##----------------------- Project: MMaker.py -------------------------------
##
##
## Purpose: This is a RScript that intends to replace the whole CodAnalyze 
##          written in 2007 which job was to build gene models that are 
##          subsequently used by our own software AMIGene. Among improvements, 
##          developements concern the full automation of the process and in 
##          particular, the classification step which in principle requires 
#           human expertise.
##
## Author: Stephane Cruveiller
## Bug reports: scruveil@genoscope.cns.fr
## First version: 0.1.0 (Nov 3, 2014)
## Current version: 1.1.0 (2 June, 2015)
##
## Dependencies:
## Python (=> 2.7.5)
## MatrixMaker.R (version => 0.7.0)
################################################################################ 

## ChangeLog
## 03/10/2014: Initial version
## 15/04/2015: - Bug fix when MatrixMaker.R is launch but not complete (the wrapper now
##               stops and throws 1 as exit code.

from __future__ import division

__author__ = 'Stephane Cruveiller'
__copyright__ = 'Copyright 2015, Stephane Cruveiller, CEA'
__license__ = 'CeCILL'
__version__ = '1.0.0'
__email__ = 'scruveil@genoscope.cns.fr'
__status__ = 'Flying Saucer'


import os
import re
import math
import csv
import sys
import errno
import pymysql
import argparse
import shlex
import time
from collections import defaultdict
from subprocess import call, Popen, PIPE
from posix import getpid

##---------------------------------------------- FUNCTION: version_string --------
##
## Purpose: This function should return a tuple containing infos about the 
##          database to connect which in principle should be set via environment
##          variables.
##
## Parameters: Nothing
##             
## Returns: a tuple containing connexion infos
##------------------------------------------------------------------------------
def version_string ():
    return __version__ +" (" + __status__ + ")"

##---------------------------------------------- FUNCTION: GetEnvVar --------
##
## Purpose: This function should return a tuple containing infos about the 
##          database to connect which in principle should be set via environment
##          variables.
##
## Parameters: Nothing
##             
## Returns: a tuple containing connexion infos
##------------------------------------------------------------------------------
def GetEnvVar():
    # returns Environment which is in principle set when using modules 
    return(os.environ.get('MYAGCUSER'),os.environ.get('MYAGCPASS'),os.environ.get('MYAGCHOST'),os.environ.get('MYAGCDB'))


##---------------------------------------------- FUNCTION: CodFreq -------------
##
## Purpose: This function initialize a tuple which will subsequently store 
##          codon counts and write out to a file .
##
## Parameters: strg (STR) A CDS string to analyze
##             fileout (HF) A handle of file indicating where to output the 
##             counting results
##             
## Returns: Nothing
##------------------------------------------------------------------------------
def CodFreq(strg,fileout):
    # initialize variables
    codon_counts = {'AAA' : 0,'AAC' : 0,'AAG' : 0,'AAT' : 0,'ACA' : 0,'ACC' : 0,'ACG' : 0,'ACT' : 0,'AGA' : 0,'AGC' : 0,'AGG' : 0,'AGT' : 0,'ATA' : 0,'ATC' : 0,'ATG' : 0,'ATT' : 0,'CAA' : 0,'CAC' : 0,'CAG' : 0,'CAT' : 0,'CCA' : 0,'CCC' : 0,'CCG' : 0,'CCT' : 0,'CGA' : 0,'CGC' : 0,'CGG' : 0,'CGT' : 0,'CTA' : 0,'CTC' : 0,'CTG' : 0,'CTT' : 0,'GAA' : 0,'GAC' : 0,'GAG' : 0,'GAT' : 0,'GCA' : 0,'GCC' : 0,'GCG' : 0,'GCT' : 0,'GGA' : 0,'GGC' : 0,'GGG' : 0,'GGT' : 0,'GTA' : 0,'GTC' : 0,'GTG' : 0,'GTT' : 0,'TAA' : 0,'TAC' : 0,'TAG' : 0,'TAT' : 0,'TCA' : 0,'TCC' : 0,'TCG' : 0,'TCT' : 0,'TGA' : 0,'TGC' : 0,'TGG' : 0,'TGT' : 0,'TTA' : 0,'TTC' : 0,'TTG' : 0,'TTT' :0}
    # Counting codons over the string sequence
    for i in range(0, len(strg) - len(strg)%3, 3):
        codon = strg[i:i+3]
        codon_counts[codon] = codon_counts.get(codon, 0) + 1
    codons = codon_counts.keys()
    codons.sort()
    for codon in codons[:-1]:
        fileout.write(str(codon_counts[codon]) + "\t")
    fileout.write(str(codon_counts['TTT']))


##---------------------------------------------- FUNCTION: Codon_Count --------
##
## Purpose: This function reads the input file containing CDSs in fasta format 
##          and outputs the codon counts in a temporary file that will be passed
##          to the R-script routine.
##
## Parameters: fnafile (STR) A string corresponding to the name of input file 
##             (i.e. CDSs in fasta format)
##             cdufile (STR) A string corresponding to the name of cdu file
##             (i.e. raw count of codon for each CDS)
##             
## Returns: a tuple containing connexion infos
##------------------------------------------------------------------------------   
def Codon_Count(fnafile,cdufile):
    # Initialize variables
    reg=re.compile("^>")
    seqbuffer=""
    countline=0

    # Print CDUfile header
    cdufile.write('Label\tAAA    AAC    AAG    AAT    ACA    ACC    ACG    ACT    AGA    AGC    AGG    AGT    ATA    ATC    ATG    ATT    CAA    CAC    CAG    CAT    CCA    CCC    CCG    CCT    CGA    CGC    CGG    CGT    CTA    CTC    CTG    CTT    GAA    GAC    GAG    GAT    GCA    GCC    GCG    GCT    GGA    GGC    GGG    GGT    GTA    GTC    GTG    GTT    TAA    TAC    TAG    TAT    TCA    TCC    TCG    TCT    TGA    TGC    TGG    TGT    TTA    TTC    TTG    TTT\n')
    
    # Now we read input file... 
    with open(fnafile,"rb") as processed_file:
        strg=processed_file.readlines()
        for line in strg:
            if reg.match(line):
                if countline == 0:
                    cdufile.write(line[1:-1]+"\t")
                    countline+=1
                else:
                    CodFreq(seqbuffer,cdufile)
                    seqbuffer=""
                    # Moving to the next header
                    cdufile.write("\n" + line[1:-1] + "\t")
                    countline+=1
            else:
                seqbuffer = seqbuffer + line[:-1]    
        CodFreq(seqbuffer,cdufile)
        cdufile.write("\n")
        print "[INFO] Number of CDSs processed: " + str(countline)
        
##---------------------------------- FUNCTION: create_Seqsatcks_fromSid --------
##
## Purpose: This function should create a file containing all CDSs for a given 
##            S_id via queries on a mysql database.
##
## Parameters: Sid (INT) -> An existing S_id from the database to query.
##               DBinfos (TUP) -> Connexion infos to access to the database (inherited
##             from environment variables  
##             
##
## Returns: Nothing (i.e. files are directly written to disk...)
##------------------------------------------------------------------------------        
def create_Seqsatcks_fromSid(Sid,DBinfos,FNAinput):
    
    # Preparing connexion to DB
    db = pymysql.connect(host=DBinfos[2],user=DBinfos[0],passwd=DBinfos[1],db=DBinfos[3])
    
    # Preparing various queries...
    # Getting all infos about Organism
    get_OrgInfos="select S_id, O_name, O_dirname,O_strain, S_name, S_length, R_genetic_code from Sequence INNER JOIN Replicon using (R_id) INNER JOIN Organism using(O_id) where S_id=" + str(Sid) + ";"
    # Getting current GOs for S_id
    get_GOs="select GO_id,GO_Begin, GO_end, GO_frame from Genomic_Object where S_id=" + str(Sid) + """ and GO_update="current" and GO_status!="artefact" and GO_type="CDS" ORDER BY GO_begin;"""
    
    # Create a crusor object to be able to query database
    cur=db.cursor()
    # ... and execute query!!!
    cur.execute(get_OrgInfos)
    # in principle we should get 0 or 1 row for the current query. So we try
    # to catch the first one if it exists!!!
    unirow=cur.fetchone()
    cur.close()
    if unirow:
        print unirow
        cur=db.cursor()
        cur.execute(get_GOs)
        gos=cur.fetchall()
        cur.close()
        try:
            with open(FNAinput,"wb") as outfile:
                for (goid, gobeg,goend,goframe) in gos:
                    print "Extracting Object " + str(goid) + ":" + str(gobeg) + "<-->" + str(goend) + "(" + goframe + ")"
                    # Extract sequences from database
                    get_GOs_Sequence="CALL NUC_SID_FASTA(" + str(Sid) + "," + str(gobeg) + "," + str(goend) + ",'" + goframe + "'," + str(goid)+ ");"
                    cur=db.cursor()
                    cur.execute(get_GOs_Sequence)
                    seq=cur.fetchall()
                    cur.close()
                    #print seq[0][0]
                    outfile.write(seq[0][0] +"\n")
        except IOError:
          print("[FATAL] Sorry, could not open file for reading/writing")  
    else:
        print "No Sid :" + str(Sid) + "in database " + DBinfos[3] + "!"

##---------------------------------- FUNCTION: prokov_learn_ISOK --------
##
## Purpose: This function should test the presence of prokov_learn used 
##          to build matrices
##
## Parameters: None  

## Returns: Nothing but exits if fails
##------------------------------------------------------------------------------  
def prokov_learn_ISOK():
    '''Tests whether prokov_learn can be successfully called.'''

    try:
        print "[INFO] Testing whether prokov_learn is installed and can be launched...\n",
        sys.stdout.flush()
        
        prokov_status=Popen(['prokov_learn','-h'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        prokov_status_out, prokov_status_err= prokov_status.communicate(b"input data that is passed to subprocess' stdin")
        retcode=prokov_status.returncode
        print "[INFO] Fine! Prokov_learn is here and runnable."
    except:
        print ("[FATAL] Could not launch prokov_learn. Have you installed it? Is it in your path?\n[FATAL] Exiting!\n")
        sys.exit(1)

##---------------------------------- FUNCTION: R_base_ISOK --------
##
## Purpose: This function should test the presence of R_base used 
##          to perform analysis and is at least at version 3.X.X
##
## Parameters: None  

## Returns: Nothing but exits if fails
##------------------------------------------------------------------- 
def R_base_ISOK():
    '''Tests whether R can be successfully called.'''

    try:
        print "[INFO] Testing whether R (V 3.x.x) is installed and can be launched ...\n",
        sys.stdout.flush()
        
        p = Popen(['R', '--version'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        rc = p.returncode
        Rversion=output.split()[2]
        print "[INFO] R version found : " + Rversion
        Rmajor=Rversion.rsplit('.')[0]
        ## Debug
        ## print "[INFO] R major version : " + Rmajor
        ## Beware Here: Rmajor is a str 
        ## Have to cast to int!!!
        if int(Rmajor) >= 3 :
            print "[INFO] Fine! R (v 3.x.x) is here and runnable."
        else :
            print "[FATAL] This software requires R version 3.0.0 and onward!\n[FATAL] Exiting!\n"
            sys.exit(1)
    except :
        print ("[FATAL] Could not find required R-base.\n[FATAL] Exiting!\n")
        sys.exit(1)        

##---------------------------------- FUNCTION: read_fasta --------
##
## Purpose: This function should read a fasta input file
##
## Parameters: fasta file to read  
##
## Returns: a Dict containing all sequences
##------------------------------------------------------------------- 
def read_fasta(fnafile):
    list=defaultdict(str)
    name = ''
    totlength=0
    for line in fnafile:
        #if your line starts with a > then it is the name of the following sequence
        if line.startswith('>'):
            name = line[1:-1]
            continue #this means skips to the next line
        #This code is only executed if it is a sequence of bases and not a name.
        list[name]+=line.strip()
        totlength+=len(list[name])
    return list,totlength

##---------------------------------- FUNCTION: print_dict --------
##
## Purpose: This function should print to stdout the content of a 
##          dict
##
## Parameters: seqs (DICT) - the dictionnary of read sequences
##
## Scope: DEBUG 
##
## Returns: Nothing
##------------------------------------------------------------------- 
def print_dict(seqs):
    for x in seqs:
        print (x) + "," + seqs[x]

##---------------------------------- FUNCTION: choose_K -------------
##
## Purpose: This function should return the appropriate K parameter 
##          (the order of the Markov Model) to be used by prokov-learn
##
## Parameters: seqs (DICT) - the dictionnary of read sequences
##
##
## Returns: k (INT) - the order of the Markov Model 
##-------------------------------------------------------------------
def choose_K(cls_sizes,cls_file) :
    ## Initialize KMM
    k=0
    
    if cls_sizes[cls_file] < 25600:
        print "[INFO] Too few sequences to build high order HMM! Defaulting to k=3."
        k=3
    elif cls_sizes[cls_file] >= 25600 and cls_sizes[cls_file] < 102400 :
        k=4
    elif cls_sizes[cls_file] >= 102400 and cls_sizes[cls_file] <= 409600 :
        k=5
    else :
        k=6
    print "[INFO] Order of the Markov Model : " + str(k) + " for file " + cls_file + "."
    ## Return value
    return(k)

##--------------------------------------------------- MAIN ROUTINE -------------
##  All jobs are done here!!!
##
def main():
    #Defines various variables
    metrics=re.compile("rscu|fmax")
    distances=re.compile("euclidean|correlation")
    pid=0  ## Process id used to prefix various outfiles
    clscds=0 ## number of CDSs for a given class
    clslgth=0 ## cumulated length for a given class
    clssizes={} ## a dict that will contain stats sizes of classes
    KMM=0 ## Order of the Markov Model
    
    
    #Defines a parser to parse cmdline arguments
    parser=argparse.ArgumentParser()
    
    # Define a mutually exclusive group (i.e. either Sid or a fna file)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s','--S_id', action='store',dest='Sid',help='S_id to process')
    group.add_argument('-i','--infile',action='store',dest='infile',help='CDS file to process')
    parser.add_argument('-n','--normalize',action='store',dest='normalized',help='Compute normalized values for codon usage. Allowed values rscu|fmax|none.')
    parser.add_argument('-e','--exclude_codons',action='store',dest='excluded_codons',help='Codons to discard from analysis (comma separated values; specified codons here will be appended to the default list (i.e. TGC,TGT,X,TAA,TAG,TGA)')
    parser.add_argument('-d','--distance',action='store',dest='dist',help='Metric to use during the data dimensions reduction. Allowed values euclidean|correlation (default).')
    # Parsing command line
    args=parser.parse_args()
    
    
    ## get general information about Process such as pid
    ## to individualize files avoiding confusions for processes run
    ## in the same working folder
    pid=getpid()
    start_time = time.time()
    ## Display init infos...
    print "[INFO] MatrixMaker analysis : " + time.strftime("%d/%m/%Y [%H:%M:%S]")
    print "[INFO] MatrixMaker version  : " + version_string()
    
    ## STEP1: 
    ## We got either a Sid or a fna file to analyze
    ## So we try to generate a Codon Usage file (in this version)
    ## TODO: See in a further verison if tuple is usable to avoid
    ## possible IO errors
    
    # STEP 1.a: let's see if we got a S_id....
    # which should be an integer!!!
    if args.Sid:
        print "[INFO] S_id to analyze: " + str(args.Sid)
        FNAfname= str(pid) + "_" + str(args.Sid) + ".fna"
        CDUfname= str(pid) + "_" + str(args.Sid) + ".cdu"
        ## TODO Remove comments of the GetEnvVar function
        #First collect environment infos
        infosdb=GetEnvVar()
        print infosdb
        #Extracting CDSs sequences from DB
        create_Seqsatcks_fromSid(args.Sid,infosdb,FNAfname)
        #Preparing output CDU file
        CDUhof=open (CDUfname,"wb") 
        Codon_Count(FNAfname,CDUhof)
        CDUhof.close()
    
    # STEP 1.b: let's see if we got a a filename....
    # for a fasta fna file!!!
    else:
        try:
            if args.infile and os.stat(args.infile)[6]!=0:
                print "[INFO] Processing file: " + args.infile
                FNAfname=args.infile
                CDUfname=str(pid) + "_" + args.infile.rsplit('.',1)[0]+".cdu"
                #Preparing output CDU file
                with open (CDUfname,"wb") as CDUhof:
                    Codon_Count(args.infile,CDUhof)
        except OSError:
            print "[FATAL] The file " + args.infile + " does not exist or is not readable!\n[FATAL] Aborting!"
            return 1
    ## At this Step, we should have a file containing raw counts
    ## of codons written on disk.
    
    ## STEP 2:
    ## Now we check whether args to be passed to MatrixMaker.R are 
    ## correctly formatted
    if args.dist and not distances.match(args.dist):
        print "[FATAL] The distance specified does not exist or is not implemented!\n[FATAL] Allowed values are: euclidean, correlation.\n"
        sys.exit(1)
    if args.normalized and not metrics.match(args.normalized):
        print "[FATAL] The normalization method specified does not exist or is not implemented!\n[FATAL] Allowed values are: rscu, fmax, none.\n"
        sys.exit(1)
    
    
    # If we reach that point this means that everything is Ok to         
    # Checking first if CDU table file exists
    try:        
        if os.stat(CDUfname)[6]!=0:
            # os.stat throws OSError when fails!!!
            #Now we check parameters for MatrixMaker.R
            if args.normalized is None:
                print "[WARN] No normalizing method specified! Defaulting to rscu..."
                args.normalized="rscu"
            if args.dist is None:
                print "[WARN] No distance specified! Defaulting to correlation distance..."
                args.dist="correlation"
            if args.excluded_codons is None:
                print "[WARN] No codon to exclude specified! Defaulting to  TGC,TGT,X,TAA,TAG,TGA"
                discarded_codons = "TGC,TGT,X,TAA,TAG,TGA"
            else:
                discarded_codons = "TGC,TGT,X,TAA,TAG,TGA," + args.excluded_codons
                print "[INFO] List of discarded codons during analysis: " + discarded_codons
            
            ## Testing if installed R is compliant!!!!
            R_base_ISOK()
            ## If yes then do computation...wooohh!!!!
            try:
                mycmd="MatrixMaker.R -c " + CDUfname + " -n " + args.normalized + " -d " + args.dist + " -e " + discarded_codons + " -p " + str(pid)
                ## Initial version here
                #call(shlex.split(mycmd))
                ## But the problem here resides here in the low level function
                ## call which does not give the possibility to handle error of
                ## subprocesses
            
                ## After complete rewriting
                runMMaker=Popen(shlex.split(mycmd),stderr=PIPE)
                runMMaker_stdout,runMMaker_stderr=runMMaker.communicate()
                # Wait until subprocess finishes otherwise we could 
                # get wrong return code!!!
                runMMaker.wait()
                retcodeMM=runMMaker.returncode
                ## If exit code from MatrixMaker.R is not 0 (i.e. successfully completed)
                ## Then raise exception and return exit code 61=errno.ENODATA=No data available!!!
                if retcodeMM !=0 :
                    raise ValueError("[FATAL] Something went wrong with MatrixMaker.R which exited with code: ", retcodeMM,"[FATAL] Possible cause : ",runMMaker_stderr.split('\n', 1)[0],"[FATAL] Aborting MMaker!")
            except ValueError as err:
                print err.args[0] + str(err.args[1]) + "\n" + err.args[2] + err.args[3] + "\n" + err.args[4]
                sys.exit(errno.ENODATA)
                
            
    except OSError:
        print "[FATAL] The codon usage table does not exist!\nAborting...\n"
        return 1
    
    ## In principle we should have somewhere files used to build the matrices
    try:
        ## First check once for all if prokov-learn is available!!!
        prokov_learn_ISOK()
        ## If we reach here then proceed...
        if os.stat(FNAfname)[6]!=0:
            with open (FNAfname,"rb") as FNAhof:
                    (Seqs,lgth)=read_fasta(FNAhof)
                    for file in os.listdir("."):
                        try:
                            if file.startswith(str(pid) + "_" + "Class_") and file.endswith(".lst"):
                                fnaoutfile=file.rsplit('.',1)[0]+".fna"
                                clscds=0
                                clslgth=0
                                KMM=0
                                with open(file,"rb") as Class_file, open(fnaoutfile,"wb") as Class_fna:
                                        for ids in Class_file:
                                            lookup=re.compile(ids.strip())
                                            for key in Seqs:
                                                if lookup.match(key):
                                                    clscds +=1
                                                    clslgth += len(Seqs[key])
                                                    Class_fna.write(">" + key + "\n")
                                                    seqsplit=[Seqs[key][i:i+60] for i in range(0, len(Seqs[key]), 60)]
                                                    Class_fna.write('\n'.join(seqsplit))
                                                    Class_fna.write("\n")
                                print "[INFO] "+ str(clscds) + " CDSs in " + file.rsplit('.',1)[0] + " for a cumulated length of " + str(clslgth) +" bp."
                                clssizes[file.rsplit('.',1)[0]] = clslgth
                                KMM=choose_K(clssizes,file.rsplit('.',1)[0])
                                try:
                                    prkcmd="prokov_learn -k " + str(KMM) + " -I -o " + file.rsplit('.',1)[0] + ".bin " + file.rsplit('.',1)[0] + ".fna"
                                    runPKLearn=Popen(shlex.split(prkcmd),stderr=PIPE)
                                    runPKLearn_stdout, runPKLearn_stderr=runPKLearn.communicate()
                                    
                                    ## Avoid simple call subprocess to enable error catching!!!
                                    ## call(shlex.split(prkcmd))
                                    
                                    # Wait until subprocess finishes otherwise we could 
                                    # get wrong return code!!!
                                    runPKLearn.wait()
                                    retcodePK=runPKLearn.returncode
                                    if retcodePK !=0 :
                                        raise ValueError("[FATAL] Something went wrong with Prokov_Learn which exited with code: ", retcodePK,"[FATAL] Possible cause : ",runPKLearn_stderr.split('\n', 1)[0],"[FATAL] Aborting MMaker!")
                                except ValueError as err:
                                    print err.args[0] + str(err.args[1]) + "\n" + err.args[2] + err.args[3] + "\n" + err.args[4]
                                    sys.exit(1)               
                        except OSError:
                            print "[FATAL] The class file could not be read!\nAborting...\n"
                            return 1
            ## In principle reaching this point ensures that the analysis
            ## has been sucessfully done... So greetings are displayed
            elapsed=time.time() - start_time
            print "[INFO] MMaker analysis successfully done...\n[INFO] Elapsed time : " + str(elapsed) + " s.\n[INFO] Bye!"
            return 0
            
    except OSError:
        print "[FATAL] The fna file has not been created!\nAborting...\n"
        return 1
        
main()
