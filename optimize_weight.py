from skopt import gp_minimize
from casm.project import Project, Selection, write_eci
import pandas as pd
import os, re, json
import numpy as np

def read_file(fname="casm_learn_input"): #read data from file and return as a pandas dataframe (header contains a hashtag)
    with open(fname) as f:
        header_line = f.readline()
    header = header_line.strip().split()[1:]
    df = pd.read_csv(fname, delim_whitespace=True, comment='#',names=header,index_col = 0)
    return df
    
def write_to_data(df,fname):
    df.to_csv(fname,sep=' ')
    with open(fname,'r') as f:
        lines=f.readlines()
        lines[0]='# '+lines[0]
    with open(fname,'w') as f:
        f.write(''.join(lines))
        
def read_log(fname): # extract number of ECI, CV and RMS error from fit_log.txt
    with open(fname,'r') as f:
        while True:        
            line = f.readline()
            if "Pickling" in line:
                break
            elif "ranged_rms:" in line:
                line = f.readline()
                rms = float(line.split()[-1])
                break
    return rms
    
def do_fit(): # call fit function externally for fitting (need optimization to use internal functions)
    os.system('rm -rf ../.casm/tmp')
    #os.system('rm -rf ../cluster_expansions/clex.formation_energy/calctype.default/ref.default/bset.default/eci.__tmp')
    os.system('rm fit_*')
    os.system("casm-learn -s fit.json > fit_log.txt")
    os.system("casm-learn -s fit.json --select 0 -q")
    os.system("casm-learn -s fit.json --checkhull > fit_log.txt")
    os.system('casm-learn -s fit.json --hall --indiv 0 --format json > fit-eci.json')
    os.system("cp ../cluster_expansions/clex.formation_energy/calctype.default/ref.default/bset.default/eci.default/eci.json .")
    cv = read_log('fit_log.txt')
    return cv

proj = Project()
sel = Selection(proj, 'casm_learn_input', all=False)
comp = 'comp(a)'
Ef = 'formation_energy'
clEf = 'clex(formation_energy)'
hull_dist = 'hull_dist(casm_learn_input,comp)'
cl_hull_dist = 'clex_hull_dist(casm_learn_input,comp)'
config = 'name'
sel.query([comp, Ef, clEf, hull_dist, cl_hull_dist])
df = sel.data
df = df[df[comp] == 0.5] # change composition here
#df = df[df["weight"] != 20]
#df = df[df[config] != "SCEL8_4_1_2_0_2_0/14"]
#df = df[(abs(df[comp] - 1/5) < 1e-8) | (df[comp] == 0.4) | (abs(df[comp] - 0.6) < 1e-8)]
states = ["SCEL1_1_1_1_0_0_0/0","SCEL8_4_1_2_0_2_3/0","SCEL8_8_1_1_0_4_4/2","SCEL6_6_1_1_0_2_3/1","SCEL8_4_1_2_0_2_0/14","SCEL3_3_1_1_0_2_0/1","SCEL4_4_1_1_0_2_2/2","SCEL1_1_1_1_0_0_0/1","SCEL6_6_1_1_0_2_3/0","SCEL7_7_1_1_0_5_5/15"]
df = df[~df[config].isin(states)]
tune_df = df[config]
print(len(tune_df))

def f(params):
    ## tune weight
    df = read_file('casm_learn_input')
    df.loc[tune_df, "weight"] = params
    
    write_to_data(df,'casm_learn_input')
    rms = do_fit()
    
    proj = Project()
    sel = Selection(proj, 'casm_learn_input', all=False)
    sel.query([config, comp, Ef, clEf])
    data = sel.data
    #data = data[data[config].isin(states)]
    data = data.set_index(config)

    data['Clex-DFT'] = data[clEf].sub(data[Ef])
    data['Clex-DFT'] = data['Clex-DFT'].apply(lambda x: x*1000) # convert number of unit to meV/f.u.
 
    return abs(data.loc['SCEL8_4_1_2_0_2_0/14', 'Clex-DFT']) #sum(data['Clex-DFT'])

to_try = (0.1,1,10,20) # weights to try
params = [to_try for i in range(len(tune_df))]
res = gp_minimize(f, params, acq_func="EI", n_calls=500, n_random_starts=10, n_jobs=-1, verbose=True) 

## refit best weights
df = read_file('casm_learn_input')
best_weights = res.x
final_cv = f(best_weights)
print(final_cv)
