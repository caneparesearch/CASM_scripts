# this code is to show energies for adjusting casm fitting
# Usage: run "python energy {0.5}" where 0.5 is the composition
import numpy as np
import os, sys
import pandas as pd
from casm.project import Project, Selection, write_eci
pd.set_option('display.max_rows', 1000)  # or 1000

proj = Project()
sel = Selection(proj, 'casm_learn_input', all=True)
comp = 'comp(a)'
Ef = 'formation_energy'
clEf = 'clex(formation_energy)'
hull_dist = 'hull_dist(casm_learn_input,comp)'
cl_hull_dist = 'clex_hull_dist(casm_learn_input,comp)'
sel.query([comp, Ef, clEf, hull_dist, cl_hull_dist])
data = sel.data.sort_values([comp])
data.set_index("name", inplace=True)

def read_file(fname): #read data from file and return as a pandas dataframe (header contains a hashtag)
    with open(fname) as f:
        header_line = f.readline()
    header = header_line.strip().split()[1:]
    df = pd.read_csv(fname, delim_whitespace=True, comment='#',names=header,index_col = 0)
    return df

print('Na_x Pb_1-x   (x = 0 ... 1)')
data['Clex-DFT'] = data['clex(formation_energy)'].sub(data['formation_energy'])
#data['comp_n(Na)'] = data['comp_n(Na)'].apply(lambda x: x/1.0) # convert number of Na to concentration x
data['Clex-DFT'] = data['Clex-DFT'].apply(lambda x: x*1000) # convert number of unit to meV/f.u.
data = data[[comp, Ef, clEf, 'selected','weight', 'Clex-DFT']]
data.columns = ['x','E_DFT','E_clex','selected','weight', 'E_clex-E_DFT (meV)']
data = data[abs(data['x'] - float(eval(sys.argv[1]))) < 1e-8]
selected_data = data[data["selected"] == True]
max_dft = selected_data.loc[selected_data['E_DFT'].idxmax()]
min_dft = selected_data.loc[selected_data['E_DFT'].idxmin()]
max_cx = selected_data.loc[selected_data['E_clex'].idxmax()]
min_cx = selected_data.loc[selected_data['E_clex'].idxmin()]
#sort data
data = data.sort_values(by=['E_DFT'],ascending=False)
print(data)
print ("min E_DFT\t", str(min_dft['E_DFT']),'\tat\t',min_dft.name,'\tweight\t',min_dft['weight'])
print ("min E_CX\t", str(min_cx['E_clex']),'\tat\t',min_cx.name,'\tweight\t',min_cx['weight'])
print ("max E_DFT\t", str(max_dft['E_DFT']),'\tat\t',max_dft.name,'\tweight\t',max_dft['weight'])
print ("max E_CX\t", str(max_cx['E_clex']),'\tat\t',max_cx.name,'\tweight\t',max_cx['weight'])
