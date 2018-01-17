# -*- coding: utf-8 -*-
'''
# Python script created by jesse greenslade Jan2018

Reads CAABA/MECCA OUTPUT
'''

### Modules ###

# plotting module, and something to prevent using displays(can save output but not display it)
import numpy as np
from datetime import datetime, timedelta
import pandas

###############
### GLOBALS ###
###############

__VERBOSE__=False




###############
### METHODS ###
###############

def read_cm_ascii(folder):
    '''
        Read the caaba_
    '''

def read_csv(filename, delimiter=',', hasheader=True):
    '''
        read a csv into a structure
        headerline is nth line read as the names for the columns
    '''
    print("Reading %s"%filename)
    #data=np.genfromtxt(filename, delimiter=delimiter, names=hasheader)

    ret={}
    with open(filename) as csvfile:
        reader=csv.DictReader(csvfile)

        for i,row in enumerate(reader):
            #print(i, row)
            # headerline is titles:
            for k in row.keys():
                if i == 0:
                    ret[k]=[]
                ret[k].append(row[k])

    return ret

