#!/usr/bin/env python
# encoding: utf-8
"""
Kegg.py

Making human's metabolic (bi-partite) network, which consists of compounds and enzymes. 
Enzymes are connected only to compounds, and vice versa. 

Created by Yong-Yeol Ahn on 2008-09-20.
"""
import sys, os, urllib, time
from collections import defaultdict

class KeggEntry:
    """Parent class for compound & enzyme"""
    def __init__(self):
        self.entry = ''
        self.names = []
        self.pathways = []
    
    def __str__(self):
        return '\t'.join([self.entry, self.names[0]])
        
class Compound(KeggEntry):
    """Compound class"""
    def __init__(self):
        KeggEntry.__init__(self)
        self.reactions = []
        self.enzymes = []
        
class Enzyme(KeggEntry):
    """docstring for Enzyme"""
    def __init__(self):
        KeggEntry.__init__(self)
        self.classes = []
        self.compounds = []
        
class MetabolicNetwork:
    """Metabolic network of human. The network is bi-partite. i.e. the enzymes are 
       only connected to the compounds and vice versa. 
    """
    def __init__(self):
        self.net = {}
        self.net['Enzyme'] = defaultdict(set)
        self.net['Compound'] = defaultdict(set)
        self.enzymeDict = {}
        self.compoundDict = {}
        
    def addNode(self, anObject):
        if anObject.__class__.__name__ == 'Enzyme':
            self.enzymeDict[anObject.entry] = anObject
        elif anObject.__class__.__name__ == 'Compound':
            self.compoundDict[anObject.entry] = anObject
        
    def addLink(self, enz, comp):
        e_entry = enz if type(enz) == type('') else enz.entry
        c_entry = comp if type(comp) == type('') else comp.entry
        
        if self.enzymeDict.has_key(e_entry) and self.compoundDict.has_key(c_entry):
            self.net['Enzyme'][e_entry].add(c_entry)
            self.net['Compound'][c_entry].add(e_entry)
    
    def writeNet(self, filename="hs_bi.net"):
        """write the network into given file."""
        f_entry = open(filename, 'w')
        f_name = open(filename[:-4] + '_name.net', 'w')
        for e_entry in self.enzymeDict.keys():
            for c_entry in self.net['Enzyme'][e_entry]:
                f_entry.write('\t'.join([e_entry, c_entry]) + '\n')
                try:
                    f_name.write('\t'.join([ self.enzymeDict[e_entry].names[0], self.compoundDict[c_entry].names[0]]) + '\n')
                except KeyError:
                    print c_entry
        
    def connectNodes(self):
        print "Number of enzymes (Homo sapiens):", len(self.enzymeDict)
        print "Number of compound (All):", len(self.compoundDict)
        
        validCompoundFlag = {}
        for e_entry, e in self.enzymeDict.iteritems():
            for c_entry in e.compounds:
                self.addLink(e_entry, c_entry)
                #valide compound c
                validCompoundFlag[c_entry] = True

        # if the compound doesn't have any connection to the enzymes of Homo sapiens, 
        # prune the compound.    
        for c_entry,c in self.compoundDict.iteritems():
            validCompoundFlag[c_entry] = validCompoundFlag[c_entry] if validCompoundFlag.has_key(c_entry) else False
            for e_entry in c.enzymes:
                if self.enzymeDict.has_key(e_entry):
                    self.addLink(e_entry, c_entry)
                    validCompoundFlag[c_entry] = True
                else:
                    continue
        
        garbages = filter(lambda x: not validCompoundFlag[x], self.compoundDict.keys())
        for c_entry in garbages:
            del(self.compoundDict[c_entry])
            
        print "Number of compound after pruning (Homo sapiens):", len(self.compoundDict)
    
    def projectTo(self, whom="Compound"):
        """One-mode projection to compounds network or enzymes network. 
        """
        pass
        
class KEGG:
    rawDataFiles = {'compound': 'ftp://ftp.genome.jp/pub/kegg/ligand/compound/compound',
                    'enzyme': 'ftp://ftp.genome.jp/pub/kegg/ligand/enzyme/enzyme'}
                    
    def __init__(self, forcedUpdate=False):
        self.update(forced=forcedUpdate)
        self.net = MetabolicNetwork()
    
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
        self.parseEnzyme()
        self.parseCompound()
        self.net.connectNodes()
    
    def parseEnzyme(self, filename="data/enzyme"):
        entries = open(filename).read().split('///')
        for aEntry in entries:
            enzyme = Enzyme()
            humanFlag = 0
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
                    if 'HSA:' in data: humanFlag = 1
                else: 
                    continue
            if not humanFlag:
                continue
            else:
                self.net.addNode(enzyme)
                                
    def parseCompound(self, filename="data/compound"):
        """docstring for parseCompound"""
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

