#!/usr/bin/env python
#========================================================================
# BlueSTARR Version 0.1
#
# Adapted from DeepSTARR by Bill Majoros (bmajoros@alumni.duke.edu)
#
#========================================================================
import gc
import gzip
import time
import math
import tensorflow as tf
import keras
import keras.layers as kl
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from keras.layers import Dropout, Reshape, Dense, Activation, Flatten
from keras.layers import BatchNormalization, InputLayer, Input, LSTM, GRU, Bidirectional, Add, Concatenate, LayerNormalization, MultiHeadAttention
import keras_nlp
from keras_nlp.layers import SinePositionEncoding
from keras import models
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.callbacks import EarlyStopping, History, ModelCheckpoint
import keras.backend as backend
from keras.backend import int_shape
import pandas as pd
import numpy as np
import ProgramName
import sys
import IOHelper
import SequenceHelper
import random
from scipy import stats
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr
from NeuralConfig import NeuralConfig
from Rex import Rex
rex=Rex()


#========================================================================
#                                GLOBALS
#========================================================================
config=None
NUM_DNA=None # number of DNA replicates
NUM_RNA=None # number of RNA replicates
#RANDOM_SEED=1234
ALPHA={"A":0,"C":1,"G":2,"T":3}
BATCH_SIZE=1

#=========================================================================
#                                main()
#=========================================================================
def main(infile,modelFilestem1,modelFilestem2,output_dir):
    #startTime=time.time()

    # Load model-DMSO
    model1=None
    with open(modelFilestem1+'.json', "r") as json_file:
        model_json=json_file.read()
        model_json = keras.saving.serialize_keras_object(model_json)
        model1 = tf.keras.models.model_from_json(model_json)
        model1.load_weights(modelFilestem1+'.h5')
 
    # Load model-Dex
    model2=None
    with open(modelFilestem2+'.json', "r") as json_file:
        model_json=json_file.read()
        model_json = keras.saving.serialize_keras_object(model_json)
        model2 = tf.keras.models.model_from_json(model_json)
        model2.load_weights(modelFilestem2+'.h5')

    # Load data
    IN=open(infile,"rt")
    with open(output_dir,"w") as file:
        for line in IN:
            fields=line.rstrip().split()
            # if(len(fields)<6): continue
            ID=fields[0]; ref=fields[1]
            
            if(not rex.find("ref=(.)",ref)):
                raise Exception("Can't parse ref: "+ref)
            ref=rex[1]
            #alleles=[fields[2],fields[4],fields[6],fields[8]]
            #seqs=[fields[3],fields[5],fields[7],fields[9]]
            alleles=[]; seqs=[]
            i=2
            while(i<len(fields)):
                alleles.append(fields[i])
                seqs.append(fields[i+1])
                seq = fields[i+1]
                i+=2
            Y1=[]
            Y2=[]
            for seq in seqs:
                X=oneHot(seq)
                X=X.reshape((1,X.shape[0],X.shape[1]))
                pred1=model1.predict(X,batch_size=1,verbose=0)
                Y1.append(pred1[0][0][0])

                pred2=model2.predict(X,batch_size=1,verbose=0)
                Y2.append(pred2[0][0][0])
                del X
            recs1=getScores(ref,alleles,Y1)
            recs2=getScores(ref,alleles,Y2)
            line=[ID]
            # for rec in recs: line.extend([str(x) for x in rec])
            for rec in recs1: line.append(str(rec))
            for rec in recs2: line.append(str(rec))
            line.append(seq)
            # print("\t".join(line))
            file.write("\t".join(line)+'\n')
            del recs1; del recs2; del fields; del line; del Y1; del Y2; del seqs; del alleles
            del ref; del ID
            gc.collect()

    # Report elapsed time
    #endTime=time.time()
    #seconds=endTime-startTime
    #minutes=seconds/60
    #print("Elapsed time:",round(minutes,2),"minutes")

#========================================================================
#                               FUNCTIONS
#========================================================================
def oneHot(seq):
    L=len(seq)
    X=np.zeros((L,4))
    for i in range(L):
        c=seq[i]
        cat=ALPHA.get(c,-1)
        if(cat>=0): X[i,cat]=1
    return X

def findRef(ref,alleles):
    n=len(alleles)
    for i in range(n):
        if(alleles[i]==ref): return i
    raise Exception("Can't find ref allele")

def getScores(ref,alleles,scores):
    r=findRef(ref,alleles)
    refScore=scores[r]
    n=len(alleles)
    recs=[]
    # modified to return reference score
    recs.append(refScore)
    # for i in range(n):
    #     if(i==r): continue
    #     log2FC=scores[i]-refScore
    #     log2FC=round(log2FC,2)
    #     rec=[alleles[i],log2FC]
    #     recs.append(refScore)
    return recs

#=========================================================================
#                         Command Line Interface
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <model-filestem1> <model-filestem2> <data> <output-dir>\n")
(modelFilestem1,modelFilestem2,infile,output_dir)=sys.argv[1:]
main(infile,modelFilestem1,modelFilestem2,output_dir)


