#Reads SPARSE output from decomp code
import numpy as np
from scipy import sparse
from sys import stdout

class site:
    elem = None
    crdsD = [None, None, None]
    crdsC = [None, None, None]

    def __init__(self, envLine, A):
        self.elem = None
        self.crdsD = [None, None, None]
        self.crdsC = [None, None, None]

        split = envLine.split()
        self.elem = int(split[1])
        self.crdsD[0] = float(split[6][1:-1])
        self.crdsD[1] = float(split[7][:-1])
        self.crdsD[2] = float(split[8][:-1])
        for i in range(0, 3): ##dot product not messed up: A is stored as
            self.crdsC[i] = 0.##A = [ [ax, ay, az], [bz, by, bz], [cx, cy, cz] ]
            for j in range(0, 3):
                self.crdsC[i] += A[j][i]*self.crdsD[j]

def SiteDistSqd(a, b):
    return  ((b.crdsC[0] - a.crdsC[0])*(b.crdsC[0] - a.crdsC[0]) + \
            (b.crdsC[1] - a.crdsC[1])*(b.crdsC[1] - a.crdsC[1]) + \
            (b.crdsC[2] - a.crdsC[2])*(b.crdsC[2] - a.crdsC[2]) + \
            (float(b.elem) - float(a.elem))*(float(b.elem) - float(a.elem)))

def SiteDistSqd_NoElem(a, b):
    return  ((b.crdsC[0] - a.crdsC[0])*(b.crdsC[0] - a.crdsC[0]) + \
            (b.crdsC[1] - a.crdsC[1])*(b.crdsC[1] - a.crdsC[1]) + \
            (b.crdsC[2] - a.crdsC[2])*(b.crdsC[2] - a.crdsC[2]))

class env:
    nSites = None
    sites = None
    fingerprint = None
    nrg = None

    def __init__(self, envLines, A):
        self.nSites = None
        self.sites = [] ##yes, this is necessary
        self.nrg = None

        ##Set sites
        self.nSites = len(envLines)
        for line in envLines:
            self.sites.append(site(envLine=line, A=A))

    def SetEnergy(self, nrg):
        self.nrg = nrg







"""
        ##Set distance array
        size = (self.nSites)*(self.nSites - 1)//2 ##// is int division
        self.distArr = [None for i in range(0, size)]
        loc = 0
        for i in range(0, self.nSites - 1):
            for j in range(i + 1, self.nSites):
                self.distArr[loc] = SiteDist(self.sites[i], self.sites[j])
                loc += 1
        self.distArr.sort()
"""

#Reads a *.JOB file and returns the important info as a dict:
#{A: [[], [], []], nSites: , nElements: , cutoffRadius: ,
#nRows: , nCols: , format: (s or d for sparse or dense)}
def ReadJobFile(infileLoc, encoding="utf-8"):
    def LatLineToLatFloat(s):
        ret = [None, None, None]
        split = s.split()
        ret[0] = float(split[6][1:])
        ret[1] = float(split[7])
        ret[2] = float(split[8][:-1])
        return ret

    ret = dict()
    with open(infileLoc, 'r', encoding=encoding) as infile:
        ##Lattice
        ret["A"] = [[None, None, None], [None, None, None], [None, None, None]]
        s = infile.readline()
        ret["A"][0] = LatLineToLatFloat(s)
        s = infile.readline()
        ret["A"][1] = LatLineToLatFloat(s)
        s = infile.readline()
        ret["A"][2] = LatLineToLatFloat(s)
        ##Number of sites, elements, cutoff radius
        ret["nSites"] = int(infile.readline().split()[3])
        ret["nElements"] = int(infile.readline().split()[3])
        ret["cutoffRadius"] = float(infile.readline().split()[2])
        infile.readline() ##empty space
        ##Matrix dimensions, matrix format
        infile.readline() ##runtime
        line = infile.readline().split()
        ret["nRows"] = int(line[9])
        ret["nCols"] = int(line[11])
        ret["format"] = infile.readline().split()[2][0]
        infile.close()

    return ret

#Reads a *.ENV file and returns a list of chemical enviornments in the order that they were read
def ReadEnvFile(infileLoc, A, encoding="utf-8"):
    ret = []
    with open(infileLoc, 'r', encoding=encoding) as infile:
        currentLines = []
        for line in infile:
            if(line[0] != 'E'): ##Every important line begins with the word "Element"
                ret.append(env(envLines=currentLines, A=A))
                currentLines = []
                continue
            currentLines.append(line)
        ##Don't miss the last enviornment
        ret.append(env(envLines=currentLines, A=A))
        infile.close()

    return ret

#Reads a *.MAT file and returns a matrix corresponding to decompositions.
#Returns a scipy COOrdinate matrix.
##note that if you've run the cell decomp program on windows, you'll probably need to use utf-16 encoding
##for some insane reason
#This function reads the sparse-matrix formatted output file (no explicit zeros)
def ReadMatFile(infileLoc, nRows=None, nCols=None,
               printProgress=True, printEvery=1000, encoding="utf-8"):
    if(printProgress):
        print("ReadInfile(): Start")
    #Read in decomposition data, fill matrix
    if (printProgress):
        print("Filling (sparse) matirx entries... ")

    rs, cs, nonzeroVals = [], [], []
    with open(infileLoc, 'r', encoding=encoding) as infile:
        for i, line in enumerate(infile):
            if(i%printEvery == 0):
                stdout.write("\r")
                stdout.write(' ' + f'{i/nRows:.2f}')
                stdout.flush()

            line_ = line.split()
            for j in range(1, len(line_), 2):
                rs.append(int(line_[0]))
                cs.append(int(line_[j]))
                nonzeroVals.append(int(line_[j+1]))

    if(printProgress):
        print("\nConstructing COOrdinate matrix ... ")
    V = sparse.coo_matrix((nonzeroVals, (rs, cs)), shape=(nRows, nCols))
    if(printProgress):
        print("ReadInfile(): End\n")
    return V
