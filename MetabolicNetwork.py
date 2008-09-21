#!/usr/bin/env python
# encoding: utf-8
"""
MetabolicNetwork.py

Created by Yong Yeol Ahn on 2008-09-21.
"""

import sys, os
from collections import defaultdict

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
    
    def writeNet(self, filename="data/hs_bi_metabolic.net"):
        """write the network into given file."""
        f_entry = open(filename, 'w')
        f_name = open(filename[:-4] + '_withname.net', 'w')
        for e_entry in self.enzymeDict.keys():
            for c_entry in self.net['Enzyme'][e_entry]:
                f_entry.write('\t'.join([e_entry, c_entry]) + '\n')
                try:
                    f_name.write('\t'.join([ self.enzymeDict[e_entry].names[0], self.compoundDict[c_entry].names[0]]) + '\n')
                except KeyError:
                    print c_entry
        
    """def writeInfo(self):
        f_enz = open('data/enzymes.info','w')
        f_cmp = open('data/compound.info','w')
        for e_entry, e in self.enzymeDict.iteritems():
            f_enz.write("" % ())
            
        for c_entry, c in self.compoundDict.iteritems():
            pass
    """
        
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
            validCompoundFlag[c_entry] = True if validCompoundFlag.has_key(c_entry) else False
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


def main():
    pass


if __name__ == '__main__':
    main()

