#!/usr/bin/env python
"""
This code contains a set of tools for manually optimizing weights for cluster expansion fitting
Author: Zeyu Deng
Modified by: Damien Lee
"""
import numpy as np 
import pandas as pd
import sys,os,re,json
from casm.project import Project, Selection, write_eci

def read_file(fname): #read data from file and return as a pandas dataframe (header contains a hashtag)
    with open(fname) as f:
        header_line = f.readline()
    header = header_line.strip().split()[1:]
    df = pd.read_csv(fname, delim_whitespace=True, comment='#',names=header,index_col = 0)
    return df

def read_log(fname): # extract number of ECI, CV and RMS error from fit_log.txt
    pattern = r'^\s*\d*:\s*\d*...\s*\d{3}\s*\d.\d*\s*\d.\d*\s*\d.\d*\s*'
    regexp = re.compile(pattern)
    with open(fname,'r') as f:
        for line in f:
            match = regexp.search(line)
            if match:
                nselect=int(line.split()[2])
                cv = float(line.split()[3])
                rms = float(line.split()[4])
                break
    return nselect,cv,rms

def add_weight_field():
    df = read_file('casm_learn_input')
    df['weight'] = np.ones((len(df),1))
    df['selected'] = np.ones((len(df),1))
    write_to_data(df,'casm_learn_input')

def tune_weight(fname,structure_id,weight): # modify casm_learn_input for weight
    df = read_file('casm_learn_input')
    print('Tune ',structure_id,' to ',str(weight),'...')
    if df.at[structure_id,'selected'] == 0.0:
        df.at[structure_id,'selected'] = 1
    if weight == 0.0:
        df.at[structure_id,'selected'] = 0
    df.at[structure_id,'weight']=float(weight)
    write_to_data(df,'casm_learn_input')
    write_log_tune(structure_id,weight)
    
def do_fit(): # call fit function externally for fitting (need optimization to use internal functions)
    os.system('rm -rf ../.casm/tmp')
    #os.system('rm -rf ../cluster_expansions/clex.formation_energy/calctype.default/ref.default/bset.default/eci.__tmp')
    os.system('rm fit_*')
    os.system("casm-learn -s fit.json > fit_log.txt")
    os.system("casm-learn -s fit.json --select 0")
    os.system("casm-learn -s fit.json --checkhull > fit_log.txt")
    os.system('casm-learn -s fit.json --hall --indiv 0 --format json > fit-eci.json')
    #os.system("casm query -k  'comp(a)' 'formation_energy' 'clex(formation_energy)' 'hull_dist(ALL,comp)'  'clex_hull_dist(ALL,comp)' 'comp_n(Na)'  -c casm_learn_input   -o data.dat")
    #os.system("cp data.dat fit_fit.dat")
    os.system("cp ../cluster_expansions/clex.formation_energy/calctype.default/ref.default/bset.default/eci.default/eci.json .")
    nselect_new,cv_new,rms_new = read_log('fit_log.txt')
    write_log(nselect_new,cv_new,rms_new)
    os.system('tail weight_log.txt')

def test_weight(fname,name): # do fit at different weight and choose the best one for minimizing rms error
    weight_list = [0.1,1,4,7,10]
    data = {}
    for weight in weight_list:
        print(' Weight = ',weight)
        tune_weight(fname,name,weight)
        do_fit()
        nselect,cv,rms = read_log('fit_log.txt')
        data[weight]=float(rms)
    best_weight = min(data,key=lambda k:data[k]) # choose the best weight that minimize rms
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('Best weight for',name,'is',best_weight)
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    tune_weight(fname,name,best_weight)
    if weight != best_weight:
        do_fit()

def write_to_data(df,fname):
    df.to_csv(fname,sep=' ')
    with open(fname,'r') as f:
        lines=f.readlines()
        lines[0]='# '+lines[0]
    with open(fname,'w') as f:
        f.write(''.join(lines))

def write_log(nselect_new,cv_new,rms_new):
    with open('weight_log.txt','a') as f:
        new_log_line = '\t'.join((('Fit:',str(nselect_new),str(cv_new),str(rms_new),'\n')))
        f.write(new_log_line)

def write_log_tune(name,weight):
    with open('weight_log.txt','a') as f:
        new_log_line = '\t'.join((('Tune:',name,str(weight),'\n')))
        f.write(new_log_line)

def test_alpha(): # test different alpha for optimizing CV score (need leaveoneout CV estimator, otherwise CV will fluctuate)
    fname = 'fit.json'
    for alpha in [1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2]:
        with open(fname,'r') as f:
            data = json.load(f)
        data['estimator']['kwargs']['alpha']=alpha
        with open(fname,'w') as f:
            f.write(json.dumps(data,indent=4))
        do_fit()
        print('alpha = ',alpha,", nselect,cv,rms = ",read_log('fit_log.txt'))
    get_diff()
    
def deselect_all(forig='casm_learn_input',select=0): # deselect all data points with 0 weight
    df_orig = read_file(forig)
    df_to_be_selected = df_orig.loc[df_orig['weight']==0.0]
    df_orig.loc[df_to_be_selected.index,'selected']=[select]*len(df_to_be_selected)
    write_to_data(df_orig,'casm_learn_input')

def transfer_weight(f_old="casm_learn_input_backup", f_new="casm_learn_input"): # transfer existing weights from old file to new file
    df = read_file(f_new)
    old_df = read_file(f_old)
    for i, row in df.iterrows():
        if i in old_df.index:
            df.loc[i, "weight"] = old_df.loc[i, "weight"]
    write_to_data(df,"casm_learn_input")

"""
Depth-first Search
Ref. Phys. Rev. B 88, 094108 (2013)

Get an initial ECI set using least-square fitting, and get LOOCV
Each candidate ECI is toggled on or off, refit ECI and re-calculate LOOCV
If new LOOCV is lower than the previous one, start the same process from this one and repeat
Should use "DirectSelection" in "feature_selection" section in fit.json
"""

def get_eci_select(fname='fit-eci.json'):
    fname = 'fit-eci.json'
    pattern = r'^\s*"selected":\s*"\d*",'
    regexp = re.compile(pattern)
    with open(fname,'r') as f:
        for line in f:
            match = regexp.search(line)
            if match:
                eci_bitstring = line.strip().split()[-1][1:-2]
                break
    eci_select_list = list(np.nonzero([int(eci) for eci in list(eci_bitstring)])[0])
    return eci_select_list

def depth_first_search():
    # initial lsq fit and get ECI
    # os.system('cp fit_lsq.json fit.json')
    # do_fit()
    nselect_init,cv_init,rms_init = read_log('fit_log.txt')
    eci_select = get_eci_select('fit-eci.json')
    # turn off each eci and recalculate CV score
    for eci in eci_select:
        print('remove ',eci,'...')
        print('now we\'re doing', eci_select)
        # write_to_json(eci_select,'fit.json')
        # do_fit()
        # nselect,cv,rms = read_log('fit_log.txt')
        # if cv < cv_init:
        #     depth_first_search()
        

    # os.system('cp fit_depth_first.json fit.json')


"""
* Number of ECI, CV and RMS error for each steps will be written into weight_log.txt

Examples:

1. Deselect all Na1.25:
select_data('casm_learn_input',1.25,0)
do_fit()  #<- don't forget this function for the results

2. Tune weight for a data point:
tune_weight('casm_learn_input','SCEL2_2_1_1_0_0_0/492',10) #<- this function can be called multiple times for tunning weights together
do_fit()  #<- don't forget this function for the results

3. Test alpha
test_alpha()

4. Test and select the best weight for optimizing RMS for a data point
test_weight('casm_learn_input','SCEL1_1_1_1_0_0_0/115')

5. Deselect all data points with 0 weight in casm_learn_input
deselect_all()
do_fit()  #<- don't forget this function for the results

6. Read data points name in to_do.txt and optimize their weights one by one
optimize_rms_jobs('to_do.txt')

7. Adding new data from casm_learn_input2 to casm_learn_input;
the new casm_learn_learn data structures should be get from
casm query -k "formation_energy corr" -c train -o casm_learn_input2
where train is obtained from
casm select --set is_calculated -o train
"""

#for i in tune_df: #['SCEL16_8_2_1_1_1_3/5']
    #tune_weight('casm_learn_input',i,1) 

tune_weight('casm_learn_input','SCEL3_3_1_1_0_1_2/0',0.1)
#tune_weight('casm_learn_input','SCEL4_4_1_1_0_3_3/0',20)
#tune_weight('casm_learn_input','SCEL4_4_1_1_0_3_3/2',20)
#tune_weight('casm_learn_input','SCEL4_4_1_1_0_3_3/1',20) 
#tune_weight('casm_learn_input','SCEL7_7_1_1_0_5_5/1',20) 
#tune_weight('casm_learn_input','SCEL5_5_1_1_0_3_4/2',20)
#point, weight = sys.argv[1], float(sys.argv[2])
#tune_weight('casm_learn_input',point,weight)
#add_weight_field()
#transfer_weight()
deselect_all()
do_fit()

