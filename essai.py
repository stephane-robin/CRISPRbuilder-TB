from pathlib import PurePath

def h_37rv():
    p_nc = str(PurePath('data', 'NC_000962.3.txt'))
    with open(p_nc, 'r') as file:
        p_nc_reader = file.read()
    return p_nc_reader

H = h_37rv()