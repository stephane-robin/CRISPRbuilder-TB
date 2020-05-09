import pathlib
from shutil import move
import xmltodict # interprets XML files like JSON files
import subprocess # allows to connect to input/output/error pipes of processes
from os import mkdir, chdir, remove, listdir, rename, getcwd, system, stat
from pickle import load, dump
from xlrd import open_workbook
from openpyxl import load_workbook
import csv
import shutil
from shutil import copyfile
from Bio import pairwise2, SeqIO
from Bio.pairwise2 import format_alignment
from Bio import Entrez # provides code to access NCBI over the Web
from datetime import datetime
import argparse
import csv
from collections import namedtuple
from pathlib import PurePath
from pathlib import Path
import pathlib
import tempfile
import time

item = 'SRR8368689'
dico_afr = {}
dico_afr[item] = {}
rep = 'REP/sequences/'+item+'/'
print("   - On blaste les spoligos")
dico_afr[item]['spoligo'] = ''
dico_afr[item]['spoligo_new'] = ''
completed = subprocess.run("blastn -num_threads 12 -query data/spoligo_old.fasta -evalue 1e-6 -task blastn -db "+rep+item+" -outfmt '10 qseqid sseqid sstart send qlen length score evalue' -out /tmp/"+item+"_old.blast", shell = True)
assert completed.returncode == 0
completed = subprocess.run("blastn -num_threads 12 -query data/spoligo_new.fasta -evalue 1e-6 -task blastn -db "+rep+item+" -outfmt '10 qseqid sseqid sstart send qlen length score evalue' -out /tmp/"+item+"_new.blast", shell = True)
assert completed.returncode == 0
print("   - On écrit les spoligos obtenus dans le fichier csv")
for pos, spol in enumerate(['old', 'new']):
    with open('/tmp/'+item+'_'+spol+'.blast') as f:
        matches = f.read()
        nb = open('data/spoligo_'+spol+'.fasta').read().count('>')
        for k in range(1,nb+1):
            if matches.count('espaceur'+spol.capitalize()+str(k)+',')>=5:
                dico_afr[item]['spoligo'+['','_new'][pos]] += '\u25A0'
            else:
                dico_afr[item]['spoligo'+['','_new'][pos]] += '\u25A1'
    dico_afr[item]['spoligo'+['','_new'][pos]+'_nb'] = [matches.count('espaceur'+spol.capitalize()+str(k)+',') for k in range(1,nb+1)]

print("     " + dico_afr[item]['spoligo'])
print("     " + str(dico_afr[item]['spoligo_nb']))
print("     " + dico_afr[item]['spoligo_new'])
print("     " + str(dico_afr[item]['spoligo_new_nb']))