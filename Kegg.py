#!/usr/bin/env python
# encoding: utf-8
"""
Kegg.py

Making human's metabolic (bi-partite) network, which consists of compounds and enzymes. 
Enzymes are connected only to compounds, and vice versa. 

Created by Yong-Yeol Ahn on 2008-09-20.
"""
import sys, os, urllib, time
from KeggEntry import *
from MetabolicNetwork import *

        
class KEGG:
    rawDataFiles = {'compound': 'ftp://ftp.genome.jp/pub/kegg/ligand/compound/compound',
                    'enzyme': 'ftp://ftp.genome.jp/pub/kegg/ligand/enzyme/enzyme'}
                    
    def __init__(self, forcedUpdate=False, organismCode="HSA"):
        """Initialization. HSA is the organism code for human."""
        self.update(forced=forcedUpdate)
        self.net = MetabolicNetwork()
        self.organismCode = organismCode
    
    def update(self, forced=False):
        """update data files. The default directory is data/."""
        if forced:
            self.retrieveFiles(self.rawDataFiles.keys())
        else:
            dataFiles = os.listdir('data/')
            toBeUpdated = filter(lambda x: x not in dataFiles, self.rawDataFiles.keys())
            self.retrieveFiles(toBeUpdated)
    
    def retrieveFiles(self, files):
        """retrieve files from KEGG database."""
        if len(files) == 0:
            print "There is no file to be updated."
        else:
            print "updating following files:", ','.join(files)
            for f in files:
                print "downloading %s..." % (f)
                urllib.urlretrieve(self.rawDataFiles[f], 'data/' + f)

    def writeNet(self):
        self.net.writeNet()
    
    def constructNetwork(self):
        """docstring for constructNetwork"""
        self.parseEnzyme(organismCode=self.organismCode)
        self.parseCompound()
        self.net.connectNodes()
    
    def parseEnzyme(self, filename="data/enzyme", organismCode="HSA"):
        print "parsing the enzyme list..."
        entries = open(filename).read().split('///')
        for aEntry in entries:
            enzyme = Enzyme()
            organismFlag = 0
            for line in aEntry.split('\n'):
                temp = line[:12].strip()
                context = temp if temp != '' else context
                data = line[12:].strip()
                if context == "ENTRY":
                    enzyme.entry = data.split()[1]
                elif context == "NAME":
                    enzyme.names.append(data.strip(';'))
                elif context == "CLASS":
                    enzyme.classes.append(data.strip(';'))
                elif context in ["SUBSTRATE", "PRODUCT", "COFACTOR"]:
                    try:
                        enzyme.compounds.append(data.split(':')[1].strip('];'))
                    except IndexError:
                        continue
                elif context == "PATHWAY":
                    try:
                        enzyme.pathways.append(data.split(':')[1].split('  '))
                    except IndexError:
                        enzyme.pathways.append(('', data.strip()))
                elif context == "GENES":
                    if organismCode in data: organismFlag = 1
                else: 
                    continue
            if not organismFlag:
                continue
            else:
                self.net.addNode(enzyme)
                                
    def parseCompound(self, filename="data/compound"):
        """docstring for parseCompound"""
        print "parsing the compound list..."
        entries = open(filename).read().split('///')
        for aEntry in entries:
            compound = Compound()
            for line in aEntry.split('\n'):
                temp = line[:12].strip()
                context = temp if temp != '' else context
                data = line[12:].strip()
                if context == "ENTRY":
                    compound.entry = data.split()[0]
                elif context == "NAME":
                    compound.names.append(data.strip(';'))
                elif context == "REACTION":
                    [compound.reactions.append(r) for r in data.split()]
                elif context == "PATHWAY":
                    try:
                        compound.pathways.append(data.split(':')[1].split('  '))
                    except:
                        compound.pathways.append(('', data.strip()))
                elif context == "ENZYME":
                    [compound.enzymes.append(e) for e in data.split()]
                else: 
                    continue
            self.net.addNode(compound)
            
def main():
	kegg = KEGG()
	kegg.constructNetwork()
	kegg.writeNet()
    
if __name__ == '__main__':
	main()

