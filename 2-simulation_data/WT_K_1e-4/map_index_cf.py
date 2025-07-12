# var1 goes from value1 to value2 with intervals of interval1

import itertools
import pandas as pd
import os
import sys
import numpy as np
import shutil

def main():
    job_id=sys.argv[1]
    var_dict={
    'K': [1e-4], #np.geomspace(1,20,8)*1e-4,
    #'K': [8e-5,], #np.geomspace(1,20,8)*1e-4,
    'int_curvature_CF': [0.5, 1.0, 2.0, 5.0],
    #'R_vert_frac': [0.4],
    #'noise_iterations':[1000, 10000, 100000, 1000000],
    'seed': np.arange(20),#np.arange(1), #np.arange(5),
    #'int_curvature_MD':[0,],#[ -10, -15 ], #here we give negative value because the orientation of the curve is clockwise in the code
    'germ_band_push' : [0.0, 0.2, 0.35, 0.4], 
    'noise_scale' : [1e-5],
    'dt' : [1e-5],
    #'number_of_mitotic_domains' : [4],
    #'K_md_by_K' : [1,],
    #'dt' : [0.000001],
    't_final':[1000.0],
    #'add_MD_bool':[0.0],
    'add_CF_bool':[1.0],
    't_CF':[0.0],
    't_MD':[0.0, 5.0], 
    }

    dirname = ''
    a= list(var_dict.values())
    combinations=list(itertools.product(*a))
    comb_df=pd.DataFrame(combinations, columns=var_dict.keys())
    comb_df = comb_df.sort_values(by=['seed']) #set order of jobs
    comb_df['seed'] = np.arange(len(comb_df))
    input_var=list(comb_df.columns)

    filelist = ["CF_methods.py", "gradDescent.jl", "initialize.py"]

    for i in range(len(comb_df)):

        folder = '_'.join([var+ '_' +str(comb_df.loc[i,var]) for var in input_var])+'/'
        comb_df.loc[i, 'folder_name'] = folder

        os.makedirs(folder, exist_ok=True)
        comb_df[comb_df.index == i].to_csv(folder+'params_to_vary.csv', index=False)

        #for file in filelist:
        #    shutil.copy(file, folder)

    comb_df.to_csv(dirname+'map_index_'+job_id+'.csv', index=False)


if __name__ == '__main__':
    main()
