#! usr/bin/env python

import argparse
import sys
import numpy as np
import pandas as pd
#from collections import defaultdict #https://stackoverflow.com/questions/5900578/how-does-collections-defaultdict-work

# NB. USARE I GIUSTI OGGETTI IN BASE ALLE ESIGENZE
# NB. USARE GLI OGGETTI DI PYTHON 
#topics
#https://www.biostars.org/p/300281/
#https://www.biostars.org/p/299866/#300021

def readVCF(vcf):
    """
    Docs:
    This function take as input a vcf file [e.g. vcf_file=open(vcf_file_path, "r")]
    and return a dictionary of lists in which the vcf is parsed (CHROM. POS. ID. REF. ALT. QUAL. FILTER. INFO. FORMAT. NORMAL. TUMOR.)
    The INFO field is further parsed with other functions
    
    N.B VEP and SnpEff annotation are still unparsed for transcripts and pipes.
    """
    filevcf = parsereadVCF(vcf)
    #filevcf ={ "CHROM" :[] , "POS" :[], "ID" :[], "REF" :[], "ALT":[], "QUAL":[],"FILTER":[], "INFO":[], "FORMAT":[], "NORMAL":[], "TUMOR":[] }
    
    vcf.seek(0)
    line=vcf.readline()
    while line != "":
        if line[0] == '#':
            pass
        else:
            vcf_fields=line.split()
            for i,key in enumerate(filevcf.keys()):
                filevcf[key].append(vcf_fields[i])
            #filevcf["CHROM"].append(vcf_fields[0])
            #filevcf["POS"].append(vcf_fields[1])
            #filevcf["ID"].append(vcf_fields[2])
            #filevcf["REF"].append(vcf_fields[3])
            #filevcf["ALT"].append(vcf_fields[4])
            #filevcf["QUAL"].append(vcf_fields[5])
            #filevcf["FILTER"].append(vcf_fields[6])
            #filevcf["INFO"].append(vcf_fields[7])
            #filevcf["FORMAT"].append(vcf_fields[8])
            #filevcf["NORMAL"].append(vcf_fields[9])
            #filevcf["TUMOR"].append( vcf_fields[10])
        line=vcf.readline()
    
    #############INFO###################
    #1
    vcf.seek(0)
    header=readHeader(vcf)
    infoDict=info2dictionary(header)
    #2
    vcf.seek(0)
    infoVal=readINFO(vcf)
    info= add2infodictionary(infoVal, infoDict)
    #https://stackoverflow.com/questions/40923429/delete-first-item-in-each-key-in-a-dictionary
    info= {k: info[k][1:] for k in info}
    filevcf.update(info)
    
    ###########FORMAT
    #1
    formatDict=format2dictionary(header)
    #print(formatDict)
    #2
    vcf.seek(0)
    formatVal = readFORMAT(vcf)
    format = add2infodictionary(formatVal, formatDict)
    format= {k: format[k][1:] for k in format}
     ##########
     #QUESTO LO POSSO FARE SOLO QUA PERCHE'non FUNZIONA CON ADD2infoDICT
    format = {k + '_normal': v for k, v in format.items()} 
     ########
    filevcf.update(format)
    
    return filevcf


def parsereadVCF(vcf):
    """
    Docs:
    ...
    """
    header=[]
    line=vcf.readline()
    while line != "":
        if line[0] == '#' and line[1] != '#':
            header.append(line.rstrip("\n"))
        line=vcf.readline()
    header[0] = header[0][1:]
    header = header[0].split("\t")
    data = {k: [] for k in header}
    return data


def readHeader(vcf):
    """
    Docs:
    This function take as input a vcf file [e.g. vcf_file=open(vcf_file_path, "r")]
    and return only the header of the vcf file
    """
    header=[]
    line=vcf.readline()
    while line != "":
        if line.startswith('#'):
            header.append(line.rstrip("\n"))
        line=vcf.readline()
    return header


def info2dictionary(header):
    """
    Docs:
    This function take as input the output of readHeader(vcf)
    and a dictionary as {"INFO_FIELD" : "Description=....."}
    """
    infoDict={}
    for line in header:
        if line.startswith("##INFO="):
                info=line
                info=info.split(",")
                infoID =info[0][11:]
                infoDesc = info[3][:-2]
                infoDict[infoID]= [infoDesc]

    return infoDict


def format2dictionary(header):
    """
    Docs:
    This function take as input the output of readHeader(vcf)
    and a dictionary as {"INFO_FIELD" : "Description=....."}
    """
    infoDict={}
    for line in header:
        if line.startswith("##FORMAT="):
                info=line
                info=info.split(",")
                infoID =info[0][13:]
                infoDesc = info[3][:-2]
                infoDict[infoID]= [infoDesc]

    return infoDict


def readINFO(vcf):
    """
    Docs:
    This function take as input a vcf file [e.g. vcf_file=open(vcf_file_path, "r")]
    and return a list of dictionaries for each vcf row.
    """
    filevcf =[]
    line=vcf.readline()
    while line != "":
        if line[0] == '#':
            pass
        else:
            vcf_fields=line.split()
            INFO=vcf_fields[7]
            INFO=INFO.split(";")
            valSplit = [val.split("=") for val in INFO]
            for val in valSplit: 
                if len(val) < 2:
                    val.append("TRUE")
            vcfDict={}
            for val in valSplit:
                vcfDict[val[0]] = val[1]
            filevcf.append(vcfDict)
        line=vcf.readline()
    return filevcf


def readFORMAT(vcf):
    """
    Docs:
    This function take as input a vcf file [e.g. vcf_file=open(vcf_file_path, "r")]
    and return a list of dictionaries for each vcf row.
    """
    filevcf =[]
    line=vcf.readline()
    while line != "":
        if line[0] == '#':
            pass
        else:
            vcf_fields=line.split()
            nc = len(vcf_fields) - 9

            d = {}
            if nc == 1:
                #WORK IN PROGRESSSSSSSS
                print("Devo fare una tpl di 1 val")
                
            elif nc == 2:
                FORMAT=vcf_fields[8]
                FORMAT=FORMAT.split(":")
                sample1 = vcf_fields[9]
                sample1 = sample1.split(":")
                sample2 = vcf_fields[9]
                sample2 = sample2.split(":")
                d=zip(FORMAT ,sample1)
                d = dict(d)
                #d = {k: (sample1,sample2) for k in FORMAT}
                filevcf.append(d)
        line=vcf.readline()
        
    
    return filevcf


def add2infodictionary(infoVal,infoDict):
    """
    Docs:
    This add values from readINFO(vcf) to info2dictionary(header)
    """
    for kd in infoDict.keys():
        #print(kd)
        for line in infoVal:
            #print(line)
            if kd in line:
                infoDict[kd].append(line[kd])
            else:
                infoDict[kd].append("NA")
    return infoDict             



def transcipts2ListVepSnpEff(vcfDict, annotator_to_split):
    """
    Docs:
    This parse the CSQ and ANN fields, adding a new key (VepCSQSplit/snpeffANNSplit) to the dictionary,
    where different transcripts are parsed to a list.
    
    Options:
    annotator_to_split = "snpeffANNSplit" [parse only ANN field]
                         "VepCSQSplit" [parse only CSQ field]
                         "snpeffANNSplit,VepCSQSplit" [parse both ANN and CSQ fields]
    """
    annotator_dict = {"snpeffANNSplit" : "ANN",  "VepCSQSplit" : "CSQ" , "snpeffANNSplit,VepCSQSplit" : "ANN,CSQ"}
    annotator = annotator_dict[annotator_to_split]
    if annotator == "ANN,CSQ":
        annotator = annotator.split(",")
        annotatorsnpeff = annotator[0]
        annotatorvep = annotator[1]
        annotator_to_split = annotator_to_split.split(",")
        annotator_to_splitsnpeff = annotator_to_split[0]
        annotator_to_splitvep =annotator_to_split[1]
        vcfDict[annotator_to_splitvep] = []
        vcfDict[annotator_to_splitsnpeff] = []       
        for key,v in vcfDict.items():
            ######################################################
            if key == annotatorvep:
                tranSplit= [val.split(",") for val in v] 
                for transplitted in tranSplit:
                    vcfDict[annotator_to_splitvep].append(transplitted)
            ######################################################
            if key == annotatorsnpeff:
                tranSplit= [val.split(",") for val in v] 
                for transplitted in tranSplit:
                    vcfDict[annotator_to_splitsnpeff].append(transplitted) 
    else:
        vcfDict[annotator_to_split] = []   
        for key,v in vcfDict.items():
            if key == annotator:
                tranSplit= [val.split(",") for val in v] 
                for transplitted in tranSplit:
                    vcfDict[annotator_to_split].append(transplitted)
    return vcfDict




def splitransciptsVepSnpEff(vcfDict, annotator_to_split):
    """
    Docs:
    This function take as input the output of transcipts2ListVepSnpEff(vcfDict)
    and separate the transcripts in different rows.
    It returns a pandas.DataFrame
    
    Options:
    annotator_to_split = "snpeffANNSplit" [parse only ANN field]
                         "VepCSQSplit" [parse only CSQ field]
                         "snpeffANNSplit,VepCSQSplit" [parse both ANN and CSQ fields]
    """
    df =  pd.DataFrame.from_dict(vcfDict)
    #Split Column containing lists into different rows in pandas
    #https://stackoverflow.com/questions/42012152/unstack-a-pandas-column-containing-lists-into-multiple-rows
    #Split a Pandas column of lists into multiple columns
    #https://stackoverflow.com/questions/35491274/split-a-pandas-column-of-lists-into-multiple-columns
    if annotator_to_split == 'VepCSQSplit':
        lst_col = 'VepCSQSplit'
        df= pd.DataFrame({col:np.repeat(df[col].values, df[lst_col].str.len()) for col in df.columns.difference([lst_col])}).assign(**{lst_col:np.concatenate(df[lst_col].values)})[df.columns.tolist()]
    elif annotator_to_split == 'snpeffANNSplit':
        lst_col = 'snpeffANNSplit'
        df = pd.DataFrame({col:np.repeat(df[col].values, df[lst_col].str.len()) for col in df.columns.difference([lst_col])}).assign(**{lst_col:np.concatenate(df[lst_col].values)})[df.columns.tolist()]    
    elif annotator_to_split == "snpeffANNSplit,VepCSQSplit":
        lst_col = 'VepCSQSplit'
        dfvep= pd.DataFrame({col:np.repeat(df[col].values, df[lst_col].str.len()) for col in df.columns.difference([lst_col])}).assign(**{lst_col:np.concatenate(df[lst_col].values)})[df.columns.tolist()]
        lst_col = 'snpeffANNSplit'
        df = pd.DataFrame({col:np.repeat(dfvep[col].values, dfvep[lst_col].str.len()) for col in dfvep.columns.difference([lst_col])}).assign(**{lst_col:np.concatenate(dfvep[lst_col].values)})[dfvep.columns.tolist()]    
    return df


def pipe2Col(vcf, vcfTransplitted, annotator_to_split):
    """
    Docs:
    This function take as input the output of splitransciptsVepSnpEff(vcfDict)
    and separate the pipe into different columns.
    It returns a pandas.DataFrame
    
    Options:
    annotator_to_split = "snpeffANNSplit" [parse only ANN field]
                         "VepCSQSplit" [parse only CSQ field]
                         "snpeffANNSplit,VepCSQSplit" [parse both ANN and CSQ fields]
    """
    annotator_dict = {"snpeffANNSplit" : "ANN",  "VepCSQSplit" : "CSQ" , "snpeffANNSplit,VepCSQSplit" : "ANN,CSQ"}
    annotatorInfo = annotator_dict[annotator_to_split]
    df = vcfTransplitted.copy()
    vcf.seek(0)
    header= readHeader(vcf)
    info = info2dictionary(header)
    if annotatorInfo == "ANN,CSQ":
        annotator = annotatorInfo.split(",")
        annotatorsnpeff = annotator[0]
        annotatorvep = annotator[1]
        annotator_to_split = annotator_to_split.split(",")
        annotator_to_splitsnpeff = annotator_to_split[0]
        annotator_to_splitvep =annotator_to_split[1]
        ########################SNPEFF###################
        annotatorsnpeff = info[annotatorsnpeff]
        annotatorSplit= str(annotatorsnpeff).split(":")
        annotatorANN= annotatorSplit[1][:-1]
        annotatorANNsplit = annotatorANN.split("|")
        annotatorANNsplit = ["SNPEFF_" + s for s in annotatorANNsplit]
        df[annotatorANNsplit] = df[annotator_to_splitsnpeff].str.split("|",expand=True,)
        
        ####################VEP####################
        annotatorvep = info[annotatorvep]
        annotatorSplit= str(annotatorvep).split(":")
        annotatorCSQ= annotatorSplit[1][:-1]
        annotatorCSQsplit = annotatorCSQ.split("|")
        annotatorCSQsplit = ["VEP_" + s for s in annotatorCSQsplit]
        df[annotatorCSQsplit] = df[annotator_to_splitvep].str.split("|",expand=True,)
    else:
        annotator = info[annotatorInfo]
        annotatorSplit= str(annotator).split(":")
        annotatorCSQ= annotatorSplit[1][:-1]
        annotatorCSQsplit = annotatorCSQ.split("|")
        #################VEP/SNPEFF###################
        
        df[annotatorCSQsplit] = df[annotator_to_split].str.split("|",expand=True,)
    return df



if __name__=="__main__":
    parser=argparse.ArgumentParser(description="vcfParser")
    parser.add_argument('-v','--vcf',action='store',type=str,help="Path of the VCF file", required=True, default=None)
    parser.add_argument('-o','--output',action='store',type=str,help="Path of the output file", required=True, default=None)
    parser.add_argument('-t','--splitAnnotator', action='store', type=str, help="Split snpEff/VEP transcripts and internal pipe separator by rows and columns ('s' for snpeff, 'v' for vep, 's,v' for both snpeff and vep)", required=False , default="None")
    args=parser.parse_args()

    vcf_file_path = args.vcf
    output_file_path = args.output
    annotator= args.splitAnnotator
    annotator_dict = { "s" : "snpeffANNSplit", "v" : "VepCSQSplit", "s,v" : "snpeffANNSplit,VepCSQSplit" , "None" : "skip"}
    
    try:
        annotator_to_splitTranscripts = annotator_dict[annotator]
    except:
        sys.exit("The annotator '-t {}' option doesn't exist".format(annotator))
    try:
        vcf_file=open(vcf_file_path, "r" )
    except:
        sys.exit("The file {} doesn't exist".format(vcf_file_path))
        

    vcfFile = readVCF(vcf_file)
    #print(annotator_to_splitTranscripts)
    ###################################BUG################################################
    #vcfFileTransplisted = transcipts2ListVepSnpEff(vcfFile)
    #print(annotator_to_splitTranscripts)
    ######################################################################################
    
    if annotator_to_splitTranscripts != "skip" :
        vcfFileTransplisted = transcipts2ListVepSnpEff(vcfFile, annotator_to_splitTranscripts)
        vcfFileTransplitted = splitransciptsVepSnpEff(vcfFileTransplisted, annotator_to_splitTranscripts)
        #devo aggiungere annotator_to_splitTranscripts
        vcfFileTransplittedPipe = pipe2Col(vcf_file, vcfFileTransplitted, annotator_to_splitTranscripts)
        vcfFileTransplittedPipe.replace(".", "NA", inplace=True)
        vcfFileTransplittedPipe.replace("", "NA", inplace=True)
        vcfFileTransplittedPipe.replace("NA", np.nan, inplace=True)
        
        print("File CSV succesfully generated")
        print("General Description of the data parsed:")
        print("********************************************Summary*************************************")
        #print(vcfFileTransplitted)
        vcfFileTransplittedPipe.info(verbose=True)
        print("****************************************************************************************")
        print("****************************************************************************************")
        #test.to_csv("test.csv")
        print("Information about the missing data:")
        print("*********************************************NA%****************************************")
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        #print(vcfFileTransplittedPipe.isnull().sum() * 100 / len(vcfFileTransplittedPipe))
        percent_missing = vcfFileTransplittedPipe.isnull().sum() * 100 / len(vcfFileTransplittedPipe)
        missing_value_df = pd.DataFrame({'column_name': vcfFileTransplittedPipe.columns,'percent_missing': percent_missing})
        missing_value_df.sort_values('percent_missing', inplace=True)
        print(missing_value_df.to_string(index=False))
        #vcfFileTransplittedPipe.dropna(axis=1, how='all', inplace=True)
        print("********************************************yussab**************************************")
        
        vcfFileTransplittedPipe.to_csv(output_file_path)
    else:
        vcfFile = pd.DataFrame.from_dict(vcfFile)
        vcfFile.replace(".", "NA", inplace=True)
        vcfFile.replace("", "NA", inplace=True)
        vcfFile.replace("NA", np.nan, inplace=True)
        
        print("File CSV succesfully generated")
        print("General Description of the data parsed:")
        print("****************************************************************************************")
        #print(vcfFile)
        vcfFile.info(verbose=True)
        print("****************************************************************************************")
        print("****************************************************************************************")
      
        print("Information about the missing data:")
        print("*********************************************NA%****************************************")
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        #print(vcfFileTransplittedPipe.isnull().sum() * 100 / len(vcfFileTransplittedPipe))
        percent_missing = vcfFile.isnull().sum() * 100 / len(vcfFile)
        missing_value_df = pd.DataFrame({'column_name': vcfFile.columns,'percent_missing': percent_missing})
        missing_value_df.sort_values('percent_missing', inplace=True)
        print(missing_value_df.to_string(index=False))
        #vcfFileTransplittedPipe.dropna(axis=1, how='all', inplace=True)
        print("********************************************yussab**************************************")
        vcfFile.to_csv(output_file_path)
        



