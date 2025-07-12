
# coding: utf-8

# In[1]:


#import os
#os.chdir('/projects/modeslab/data_version_5/WT_K_1e-4/') #


# In[2]:


from CF_methods import *


# In[3]:


groupby_cols = ["K", "germ_band_push", 'int_curvature_CF']
path_to_dfs = "dfs/"
os.makedirs(path_to_dfs, exist_ok=True)


# In[4]:


def get_fold_analysis(row, t = None, colName_suffix = "", fold_threshold = 0.035):
    
    path = row["path"]
    
    #####
    # compute Select a timepoint
    #####
    
    df_energy = pd.read_csv(path + "summed_prop_per_timepoint.csv")
    df_all_times = pd.read_csv(path + "database.csv")
    times = df_energy["t"].values
    
    if (t is None) and (colName_suffix == "_at_final" ):
        t = max(times)
            
    if (t is None) and (colName_suffix == "_at_maxBendEnergy" ):
        if row["nb_folds_at_final"] == 0: #would be nice to. check here if "nb_folds_at_final" is in columns
            #if nb of folds is zero then we take the last timepoint
            t = max(times)
        else:
            #we take the last t if there are multiple times where max bend energy is achieved
            t = df_energy.loc[df_energy["bend_energy"] == max(df_energy["bend_energy"]),"t"].values[-1] 
        #row["t"+colName_suffix] = t
        
    if (t is None) and (colName_suffix == "_at_avgMaxBendEnergy" ):
        K = row["K"]
        g = row["germ_band_push"]
        t_avgMaxBendEnergy = group_K_g.loc[(group_K_g["K"] == K) & (group_K_g["germ_band_push"] == g),"t_at_maxBendEnergy"].values[0]
        t = find_nearest(array=times, value=t_avgMaxBendEnergy)
        #row["t"+colName_suffix] = t
        
    if (t is None) and (colName_suffix == "_at_maxBendEnergy_afterMD" ):
        t_MD = row["t_MD"]
        if row["nb_folds_at_final"] == 0: #would be nice to. check here if "nb_folds_at_final" is in columns
            #if nb of folds is zero then we take the last timepoint
            t = max(times)
        else:
            #we take the last t if there are multiple times where max bend energy is achieved
            df_energy_temp = df_energy[df_energy["t"]>t_MD]
            if len(df_energy_temp) == 0:
                t = max(times)
            else:
                t = df_energy_temp.loc[df_energy_temp["bend_energy"] == max(df_energy_temp["bend_energy"]),"t"].values[-1] 
        
    if (t is not None) and (colName_suffix != ""):
        row["t"+colName_suffix] = t
        
    #get the timepoint
    df = df_all_times[df_all_times["t"] == t]
    df["dist_to_ellipse"] = df.apply(dist_to_ellipse,axis = 1)
    
    #segmenting folds
    df["fold_mask"] = (df["dist_to_ellipse"] > fold_threshold).astype('float')
    df["fold_mask_left_boundary"] = (df["fold_mask"].diff() > 0).astype('float')
    df["fold_id"] = df["fold_mask"]*(df["fold_mask_left_boundary"].cumsum())
    
    #counting number and size of folds
    nb_folds = len(np.unique(df["fold_id"])) - 1 #non-folds have fold id 0 which is a unique id
    df["fold_width_ref"] = df.groupby("fold_id")["fold_id"].transform(lambda x : np.sum(df.loc[x.index, "dr0"])*df.loc[x.index, "fold_mask"])
    df["fold_width_current"] = df.groupby("fold_id")["fold_id"].transform(lambda x : np.sum(df.loc[x.index, "d_r_mod_alternate"])*df.loc[x.index, "fold_mask"])
    
    #checking if a fold is inside the CF
    CF_ind = df[df["ko"] != 0].index
    if len(CF_ind) > 0:
        max_CF_arclength_frac = max(df.loc[CF_ind,"arclength_frac"])
        min_CF_arclength_frac = min(df.loc[CF_ind,"arclength_frac"])
        df["fold_inside_CF_bool"] = df.groupby("fold_id")["fold_id"].transform(lambda x : (float((max(df.loc[x.index,"arclength_frac"]) > max_CF_arclength_frac) & (min(df.loc[x.index,"arclength_frac"]) < min_CF_arclength_frac)))*df.loc[x.index, "fold_mask"])    
        fold_id_inside_CF = np.sort(np.unique(df["fold_id"]*df["fold_inside_CF_bool"]))[-1] #sort and then take last value to ensure fold_id is caught (if no fold in CF then 0 is returned)
        CF_area_of_influence_ref = df[df["fold_id"] == fold_id_inside_CF]["fold_width_ref"].values[0]
        CF_area_of_influence_current = df[df["fold_id"] == fold_id_inside_CF]["fold_width_current"].values[0]
    else:
        df["fold_inside_CF_bool"] = 0
        CF_area_of_influence_ref = 0
        CF_area_of_influence_current = 0
        
    #we can also check if folds are intersecting with MD or are they at boundary of MD?
    #is the diff of MD_id inside a fold
    #is the maxima of fold inside a MD?
        
    #saving this timepoint
    df["CF_area_of_influence_ref"] = CF_area_of_influence_ref
    df["CF_area_of_influence_current"] = CF_area_of_influence_current
    df.to_csv(path + "balls_w_fold_analysis" + colName_suffix + ".csv", index = False)
    
    row["nb_folds" + colName_suffix] = nb_folds
    row["avg_depth" + colName_suffix] = np.mean(df["dist_to_ellipse"])
    row["max_depth" + colName_suffix] = np.max(df["dist_to_ellipse"])
    row["CF_area_of_influence_ref" + colName_suffix] = CF_area_of_influence_ref
    row["CF_area_of_influence_current" + colName_suffix] = CF_area_of_influence_current
    
    return(row)

###################


# In[5]:


########
# Reading the fold threshold #
########

try:
    f = open("fold_threshold.txt", "r")
    fold_threshold = float(f.read())
    print("reading fold threshold from text file")
except:
    fold_threshold = 0.035
    print("could not file file to read fold threshold")


# In[6]:


#%%time
#preparing a dataframe with one row for each simulation with all parameters
grouped_df = pd.read_csv(glob.glob("map_index_*.csv")[0])
grouped_df = grouped_df.rename(columns={"folder_name": "path"})
sim_ran_bool = grouped_df.apply(lambda row: os.path.exists(row["path"]+"database.csv"),axis = 1)
grouped_df = grouped_df[sim_ran_bool]
print(f"{len(grouped_df)} simulations ran from {len(sim_ran_bool)} simulations")


# In[7]:


#%%time 
print("fold analysis at final timepoint")
#this can take time
grouped_df = grouped_df.apply(get_fold_analysis, t = None, colName_suffix = "_at_final", fold_threshold = fold_threshold, axis = 1)


# In[8]:


#%%time
print("fold analysis at max bend energy after MD")
grouped_df = grouped_df.apply(get_fold_analysis, t = None, colName_suffix = "_at_maxBendEnergy", fold_threshold = fold_threshold, axis = 1)


# In[9]:


#%%time
print("group simulation by K and g")
#make sure every other parameter that you want to keep constant is constant
#we will average over different seed values
group_K_g = grouped_df.groupby(groupby_cols).agg(["mean", "std", "max"]).reset_index()
colnames = [x[0] if ((x[1] == "mean") or (x[0] in groupby_cols)) else x[0]+'_'+x[1] for x in group_K_g.columns]
group_K_g.columns = colnames #removing multi-indexing
#group_K_g.to_csv("dfs/avg_prop_over_K_and_g_combination.csv",index = False)


# In[10]:


#%%time
print("save avg_prop_over_K_and_g_combination.csv")
#saving csvs
group_K_g.to_csv(path_to_dfs + "avg_prop_over_K_and_g_combination.csv",index = False)
grouped_df.to_csv(path_to_dfs + "avg_prop_over_each_sim.csv",index = False)
#group_K_g = pd.read_csv("dfs/avg_prop_over_K_and_g_combination.csv")
#grouped_df = pd.read_csv("dfs/avg_prop_over_each_sim.csv")


# In[26]:


#%%time
print("combine energy of simulations")
df_energy_all = combine_files_from_folders(glob_path = "*/", filename="summed_prop_per_timepoint.csv",verbose = False)
df_energy_all.to_csv(path_to_dfs + "combined_energy.csv", index=False)


# In[27]:


#%%time
print("group energy measurements")
#here we know that t values are all int but might differ in the 5th decimal or so
df_energy_all["t"] = df_energy_all["t"].astype("int")
groupby_cols_temp = ['t'] + groupby_cols
#here we pool discs within a devstage and calculate the mean and std
grouped_energy_df = df_energy_all.groupby(groupby_cols_temp).agg(['mean', 'std']).reset_index()
colnames = [x[0]+'_'+x[1] if x[0] not in groupby_cols_temp else x[0] for x in grouped_energy_df.columns]
grouped_energy_df.columns = colnames #removing multi-indexing
grouped_energy_df.to_csv(path_to_dfs + "grouped_energy_df.csv",index=False)


# In[28]:


#%%time
print("combining all simulations at max Bend Energy")
df_allSims_at_maxBendEnegy = combine_files_from_folders(glob_path = "*/", filename="balls_w_fold_analysis_at_maxBendEnergy.csv", verbose = False)
df_allSims_at_maxBendEnegy.to_csv(path_to_dfs + "allSims_at_maxBendEnegy_aferMD.csv",index = False)


# In[29]:


#%%time
print("getting the average depth along the curve at max Bend Energy")
groupby_cols_temp = ['ID'] + groupby_cols
df_mean_shape_at_maxBendEnegy = df_allSims_at_maxBendEnegy.groupby(groupby_cols_temp).agg(["mean", "std"]).reset_index()
colnames = [x[0] if ((x[1] == "mean") or (x[0] in groupby_cols_temp)) else x[0]+'_'+x[1] for x in df_mean_shape_at_maxBendEnegy.columns]
df_mean_shape_at_maxBendEnegy.columns = colnames #removing multi-indexing
df_mean_shape_at_maxBendEnegy.to_csv(path_to_dfs + "mean_depth_along_curve_at_maxBendEnergy.csv", index = False)


# In[30]:


#%%time
print("combining all simulations at final timepoint")
df_allSims_at_final = combine_files_from_folders(glob_path = "*/", filename="balls_w_fold_analysis_at_final.csv", verbose = False)
df_allSims_at_final.to_csv(path_to_dfs + "allSims_at_final.csv",index = False)


# In[31]:


#%%time
print("getting the average depth along the curve at final timepoint")
groupby_cols_temp = ['ID'] + groupby_cols
df_mean_shape_at_final = df_allSims_at_final.groupby(groupby_cols_temp).agg(["mean", "std"]).reset_index()
colnames = [x[0] if ((x[1] == "mean") or (x[0] in groupby_cols_temp)) else x[0]+'_'+x[1] for x in df_mean_shape_at_final.columns]
df_mean_shape_at_final.columns = colnames #removing multi-indexing
df_mean_shape_at_final.to_csv(path_to_dfs + "mean_depth_along_curve_at_final.csv", index = False)

print("Generated all dfs. Happy plotting!")
