#!/usr/bin/env python
# encoding: utf-8
"""
KeggEntry.py

Created by Yong Yeol Ahn on 2008-09-21.
"""

import sys, os

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


def main():
    pass

if __name__ == '__main__':
    main()

