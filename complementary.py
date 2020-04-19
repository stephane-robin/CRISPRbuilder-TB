# attention, désinstallation de subprocess.run 0.0.8 qui empêche le
# fonctionnement


"""
fastq : from INSDC archive NCBI
fastq-dump --fasta --split-files SRR8368696
for item in manilla:
    system('mkdir'+item)

REP = '/home/christophe/Documents/codes/66.MTBC_old_DNA/'
REP = '/home/christophe/Documents/codes/66.MTBC_Kinshasa/'
REP = '/run/media/christophe/76005607-4d83-44fd-9550-4027f19e822b/'
REP = '/run/media/christophe/b02a1067-f287-4665-a4c2-9129853c0b26/'
REP = '/media/christophe/MTBC/'
"""


# ===============
# DATABASE OF SRA
# ===============

# Mycobacterium tuberculosis variant caprae ERR1462634
# Mycobacterium tuberculosis variant microti SRR3647357
# H37Rv ERR305600
Nouvel_article = ['ERR3335723', 'ERR3335724', 'ERR3335725', 'ERR3335726', 'ERR3335727', 'ERR3335728', 'ERR3335729', 'ERR3335730', 'ERR3335731', 'ERR3335732', 'ERR3335733', 'ERR3335734', 'ERR3335735', 'ERR3335736', 'ERR3335737', 'ERR3335738', 'ERR3335739', 'ERR3335740', 'ERR3335741', 'ERR3335742', 'ERR3335743', 'ERR3335744', 'ERR3335745', 'ERR3335746', 'ERR3335747', 'ERR3335748', 'ERR3335749', 'ERR3335750', 'ERR3335751', 'ERR3335752', 'ERR3335753', 'ERR3335754', 'ERR3335755', 'ERR3335756', 'ERR3335757', 'ERR3335758', 'ERR3335759', 'ERR3335760', 'ERR3335761', 'ERR3335762', 'ERR3335763', 'ERR3335764', 'ERR3335765', 'ERR3335766', 'ERR3335767', 'ERR3335768', 'ERR3335769', 'ERR3335770']
marinum  = ['SRR8368696', 'SRR8368678', 'SRR8368689']
Momies = ['ERR650569', 'ERR650970', 'ERR650971', 'ERR650973', 'ERR650974', 'ERR650975', 'ERR650976', 'ERR650977', 'ERR650978', 'ERR650979', 'ERR650980', 'ERR650981', 'ERR650982', 'ERR650983', 'ERR650984', 'ERR650985', 'ERR650986', 'ERR650987', 'ERR650988', 'ERR650990', 'ERR650991', 'ERR650992', 'ERR650993', 'ERR650994', 'ERR650995', 'ERR650996', 'ERR650997', 'ERR650998', 'ERR650999', 'ERR651000', 'ERR651001', 'ERR651002', 'ERR651003', 'ERR651005', 'ERR651006', 'ERR651007', 'ERR651008', 'ERR651009', 'ERR651010', 'ERR650989', 'ERR651004', 'ERR650972']
Manilla = ['DRR099684', 'DRR099686', 'DRR099689', 'DRR099692', 'DRR099683', 'DRR099685', 'DRR099687', 'DRR099688', 'DRR099690', 'DRR099691', 'DRR099693', 'DRR099694']
#for item in Manilla:
#    mkdir('REP/sequences/'+item+'/')

# we browse the list of directories and files in 'REP/sequences' called item
# and define a path called 'rep' for each of them. Each 'item' represents a SRA.
#for item in listdir(REP + 'sequences/'):
    #rep = REP + 'sequences/' + item + '/'

# TODO check the relevance of this function. It's not being used
def entrez_to_dico(dico, loc='', Strain=''):
    """
    This function takes data form dico and creates a dictionary called
    dico_afr.

    Args
        dico (dict): dictionary of experimental MTBC references
        loc (str): location of an experimental MTBC reference
        Strain (str): strain of an experimental MTBC reference

    Returns:
        dico_afr (dict): the structure of dico_afr is the following
        {
            a_SRR:
                {
                    'accession': a_value,
                    'location': a_value,
                    'date': a_value,
                    'SRA': a_value,
                    'center': a_value,
                    'strain': a_value,
                    'taxid': a_value,
                    'name': a_value,
                    'spoligo': ,
                    'spoligo_new': ,
                    'SIT': ,
                    'lineage_Coll': ,
                    'IS6110': ,
                    'desc': a_value,
                    'title': a_value
                },
            ...
        }

    Note:
        - if the dico 'platform' is 'ILLUMINA' and the dico 'LIBRAIRY_LAYOUT' is
          'PAIRED', then we assign 'SAMPLE_ATTRIBUTE' to a temporary list called
          'attributes', which we'll use to assign values to dico_afr. We
          browse the list 'attributes' to define dico0['location'],
          dico0['date'], dico0['sra'], dico0['center'], dico0['strain'],
        - if dico['SAMPLE_ATTRIBUTE'] is empty, then dico_afr stays empty,
        - if the parameters loc and Strain are defined when the function is
          called, then we assign loc to location and Strain to strain, instead
          of  their values from 'attributes',
        - we rename dico[...]['RUN'] into SRR, which is either a list of
          dictionaries or a dictionary. We transform SRR into a list of values
          called srr, separating the cases when SRR is dictionary or a list of
          dictionaries. Those srr correspond to former SRR['@accession'],
           which are actually real SRR numbers (out of the context of the code).
        - for each srr in SRR, we assign the different values found above to the
          dictionary dico_afr[srr]

    """
    dico_afr = {}
    location, date, sra, center, strain = '', '', '', '', ''

    if 'ILLUMINA' in dico['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE'][
        'EXPERIMENT']['PLATFORM'] and 'PAIRED' in dico[
        'EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['EXPERIMENT'][
        'DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_LAYOUT']:

        try:
            attributes = dico['EXPERIMENT_PACKAGE_SET'][
                'EXPERIMENT_PACKAGE']['SAMPLE']['SAMPLE_ATTRIBUTES'][
                'SAMPLE_ATTRIBUTE']
        except:
            return {}

        for k in attributes:
            if k['TAG'] == 'geographic location (country and/or sea)':
                location = k['VALUE']
            elif k['TAG'] == 'collection date':
                date = k['VALUE']
            elif k['TAG'] == 'SRA accession':
                sra = k['VALUE']
            elif k['TAG'] == 'INSDC center name':
                center = k['VALUE']
            elif k['TAG'] == 'Strain':
                strain = k['VALUE']

        if loc != '':
            location = loc
        if Strain != '':
            strain = Strain

        SRR = dico['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE'][
            'RUN_SET']['RUN']
        if not isinstance(SRR, list):
            SRR = [SRR['@accession']]
        else:
            SRR = [u['@accession'] for u in SRR]

        for srr in SRR:
            if srr not in dico_afr:
                dico_afr[srr] = {
                    'accession': srr,
                    'location': location,
                    'date': date,
                    'SRA': sra,
                    'center': center,
                    'strain': strain,
                    'taxid': dico['EXPERIMENT_PACKAGE_SET'][
                        'EXPERIMENT_PACKAGE']['SAMPLE'][
                        'SAMPLE_NAME'].setdefault('TAXON_ID', ''),
                    'name': dico['EXPERIMENT_PACKAGE_SET'][
                        'EXPERIMENT_PACKAGE']['SAMPLE'][
                        'SAMPLE_NAME'].setdefault('SCIENTIFIC_NAME', ''),
                    'SIT': '',
                    'spoligo': '',
                    'spoligo_new': '',
                    'lineage_Coll': '',
                    # 'IS_mapper': '',
                    'IS6110': ''
                }

                """
                # Ajout des MIRU
                for mir in ['MIRU02', 'Mtub04', 'ETRC', 'MIRU04', 'MIRU40', 'MIRU10', 'MIRU16', 'Mtub21', 'MIRU20', 'QUB11b', 
                            'ETRA', 'Mtub29', 'Mtub30', 'ETRB', 'MIRU23', 'MIRU24', 'MIRU26', 'MIRU27', 'Mtub34', 'MIRU31', 
                            'Mtub39', 'QUB26', 'QUB4156', 'MIRU39']:
                    dico_afr[srr][mir] = ''
                """
                '''
                dico_afr[srr]['title'] = dico['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['STUDY']['DESCRIPTOR']['STUDY_TITLE']
                desc = '' 
                for k in ['STUDY_ABSTRACT', 'STUDY_DESCRIPTION']:
                    if k in dico['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['STUDY']['DESCRIPTOR']:
                        txt = dico['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['STUDY']['DESCRIPTOR'][k]
                        if txt != desc:
                            desc = txt
                dico_afr[srr]['desc'] = desc
                '''

    return dico_afr


"""
A RAJOUTER


A_rajouter = []
#[('SRP131742','')]
if A_rajouter != []:
    for k in A_rajouter:
        print('-',k[0])
        handle = Entrez.esearch(db="sra", retmax=10000, term=k[0])
        X = Entrez.read(handle)['IdList']
# On parcourt les SRA trouvés, et on enrichit le dictionnaire des spoligos déjà obtenus.

for ids in X:
    print(ids)
    ret = Entrez.efetch(db="sra", id=ids, retmode="xml")
    dico = xmltodict.parse(ret.read())
    dico_afr = dict(dico_afr, **entrez_to_dico(dico, loc='', Strain=k[1]))
# Taxonomy browser : https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
handle = Entrez.esearch(db="sra", retmax=10000, term="txid78331[Organism:exp]")
X = Entrez.read(handle)['IdList']
# On parcourt les SRA trouvés, et on enrichit le dictionnaire des spoligos déjà obtenus.
for ids in X:
    print(ids)
    ret = Entrez.efetch(db="sra", id=ids, retmode="xml")
    dico = xmltodict.parse(ret.read())
    dico_afr = dict(dico_afr, **entrez_to_dico(dico, loc='', Strain=k[1]))

# On rajoute les SRA du fichier TB_Roychowdhury_s2
wb = open_workbook('data/TB_Roychowdhury_s2.xls')
ws = wb.sheet_by_index(1)
# ou by_names
for row in range(1, ws.nrows):
    acc = ws.cell_value(row, 0).split('_')[0]
    print(acc)
    if acc not in dico_afr:
        print(' -> A rajouter')
        handle = Entrez.esearch(db="sra", retmax=10000, term=acc)
        X = Entrez.read(handle)['IdList']
        # On parcourt les SRA trouvés, et on enrichit le dictionnaire des spoligos déjà obtenus.
for ids in X:
    ret = Entrez.efetch(db="sra", id=ids, retmode="xml")
    dico = xmltodict.parse(ret.read())
    dico_afr = dict(dico_afr, **entrez_to_dico(dico))
    if fichier in listdir('data/'):
        with open(REP+'data/'+fichier, 'rb') as f:
            dico_afr = load(f)
            # On rajoute les SRA du fichier NIHMS512109
wb = load_workbook(filename='data/NIHMS512109-supplement-3.xlsx', read_only=True)
ws = wb['Hoja1']
for row in ws.iter_rows(row_offset=5):
# Est-ce du paired-end ?
    if row[2].value == 'PE':
        print(row[0].value)
        handle = Entrez.esearch(db="sra", retmax=10000, term="(txid1773[Organism])and "+row[0].value)
        X = Entrez.read(handle)['IdList']
        # On parcourt les SRA trouvés, et on enrichit le dictionnaire des spoligos déjà obtenus.
for ids in X:
    ret = Entrez.efetch(db="sra", id=ids, retmode="xml")
    dico = xmltodict.parse(ret.read())
    dico_afr = dict(dico_afr, **entrez_to_dico(dico, loc=row[7].value, Strain=row[0].value)) else:
    dico_afr = {}
    # On récupère la liste des SRA (Sequence Read Archive) Mycobacterium tuberculosis variant africanum :
handle = Entrez.esearch(db="sra", retmax=10000, term="txid33894[Organism]")
X = Entrez.read(handle)['IdList']
# On parcourt les SRA trouvés, et on enrichit le dictionnaire des spoligos déjà obtenus.
for ids in X: ret = Entrez.efetch(db="sra", id=ids, retmode="xml")
dico = xmltodict.parse(ret.read())
dico_afr = dict(dico_afr, **entrez_to_dico(dico)) with open(REP+'data/'+fichier, 'wb') as f:
    dump(dico_afr, f)

FIN A RAJOUTER
"""
