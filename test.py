import xmltodict  # interprets XML files like JSON files
import subprocess  # allows to connect to input/output/error pipes of processes
from os import mkdir, chdir, remove, listdir, rename, getcwd, system, \
    stat  # system
from pickle import load, dump
from xlrd import open_workbook
from openpyxl import load_workbook
import csv
import shutil  # TESTER UTILITE D'IMPORTER TOUT LE MODULE
from shutil import copyfile
from Bio import pairwise2, SeqIO
from Bio.pairwise2 import format_alignment
from Bio import Entrez  # provides code to access NCBI over the Web
from datetime import datetime


def to_Brynildsrud():
    """
    This function extracts data from 'data/Brynildsrud_Dataset_S1.xls'
    representing the Brynidsrud lineage and put them in a dictionary called
    Brynildsrud.

    Returns:
        Brynildsrud (dict):

    Note:


    """
    wwb = open_workbook('data/Brynildsrud_Dataset_S1.xls')
    wws = wwb.sheet_by_index(0)

    Brynildsrud = {}
    source, author, study, location, date = '', '', '', '', ''

    for row in range(1, wws.nrows):
        srr = wws.cell_value(row, 4)

        if len(srr) > 1:
            Brynildsrud[srr] = {}

            source = wws.cell_value(row, 5).replace(',', '')

            if source == 'This study':
                author = 'Brynildsrud et al.'

            source = source.replace('This study', 'Global expansion of '
                                'Mycobacterium tuberculosis lineage 4 shaped by '
                                'colonial migration and local adaptation')

            if len(wws.cell_value(row, 3)) > 1:
                study = wws.cell_value(row, 3)

            if len(wws.cell_value(row, 6)) > 1:
                location = wws.cell_value(row, 6)

            if wws.cell_value(row, 0).count('_') == 2:
                dat = wws.cell_value(row, 0).split('_')[-1]
                une_date = True
                for w in dat:
                    if w not in '0123456789':
                        une_date = False
                if une_date:
                    date = dat

    Brynildsrud[srr] = {
        'Source': source,
        'Author': author,
        'study accession number': study,
        'location': location,
        'date': date
    }
    return Brynildsrud

print(to_Brynildsrud())
