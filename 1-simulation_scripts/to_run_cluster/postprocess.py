#import packages
from CF_methods import *

#set working directory
map_index_dest=sys.argv[1]
task_id=int(sys.argv[2])
map_index=pd.read_csv(map_index_dest)
input_var=list(map_index.columns) #variables for which values are being imported
for i in range(len(input_var)):
    if input_var[i] == 'folder_name':
        dirname = map_index.loc[task_id,str(input_var[i])]
        break
    #var_dict[input_var[i]]=map_index.loc[task_id,str(input_var[i])]

#dirname=sys.argv[1]
os.chdir(dirname)

#get timepoints
files = glob.glob("N_iter_*.csv")

if len(files) > 0:
    print("compiling all timepoints")
    times_str = [re.search('N_iter_(.*).csv',x).group(1) for x in files]
    times = []
    for i in range(len(times_str)):
        if len(times_str[i]) == 0:
            continue
        times.append(int(times_str[i]))
    times = np.sort(times)

    #####################################################
    #combining different csv files into one single file #
    #####################################################
    
    df_common = pd.read_csv("init_balls.csv")
    #df_common = df_common[["ID", "mitotic_domain", "mitotic_domain_id","dr0", "K", "L+1", "L-1","ko","arclength", "arclength_frac"]]
    
    database = pd.DataFrame()
    #loop in timepoints to make combine datafraem
    for t in times:
	    df_temp = pd.read_csv("N_iter_"+str(t)+'.csv')
	    #df_temp = pd.concat([df_common, df_temp], axis = 1)
	    database = pd.concat([database, df_temp])
    #save database
    database.to_csv("database.csv", index = False)
    #loop in timepoints to delete file
    for t in times:
	    os.remove("N_iter_"+str(t)+'.csv')

#####################################################


######################
# Calculating energy #
######################

#reading the file
database = pd.read_csv("database.csv")
#calculation
database["bend_energy"] = database["K"]*((database["curvature"] - database["ko"])**2)*database["dr0"]
database["stretch_energy"] = database["L+1"]*( (database['d_r_mod'] - database["dr0"])**2 )/(database["dr0"])
df_energy = database[["N_iter", "t", "bend_energy", "stretch_energy"]].groupby(["N_iter", "t"]).agg('sum').reset_index()
df_energy = df_energy.sort_values(by = ["N_iter"]).reset_index(drop = True)
df_energy["total_energy"] = df_energy["bend_energy"] + df_energy["stretch_energy"]
df_energy["bend_energy_normalized"] = df_energy["bend_energy"]/df_energy["total_energy"][0]
df_energy["stretch_energy_normalized"] = df_energy["stretch_energy"]/df_energy["total_energy"][0]
df_energy["total_energy_normalized"] = df_energy["total_energy"]/df_energy["total_energy"][0]

#save the files
database.to_csv("database.csv", index = False)
df_energy.to_csv("summed_prop_per_timepoint.csv", index = False)

######################


###############
# Plotting ####
###############

#plot
fig,axs = plt.subplots(2,3, figsize = (18,10))

ax = axs[0,0]
fig,ax = get_subplot(fig,ax,df_energy["t"], df_energy["stretch_energy"],
                 xlabel = "t", ylabel = "stretch"
                )
ax = axs[0,1]
fig,ax = get_subplot(fig,ax,df_energy["t"], df_energy["bend_energy"],
                 xlabel = "t", ylabel = "bend"
                )
ax = axs[0,2]
fig,ax = get_subplot(fig,ax,df_energy["t"], df_energy["total_energy"],
                 xlabel = "t", ylabel = "total"
                )
ax = axs[1,0]
fig,ax = get_subplot(fig,ax,df_energy["t"], df_energy["stretch_energy_normalized"],
                 xlabel = "t", ylabel = "stretch_norm", ylim = (-0.01, 2.01)
                )
ax = axs[1,1]
fig,ax = get_subplot(fig,ax,df_energy["t"], df_energy["bend_energy_normalized"],
                 xlabel = "t", ylabel = "bend_norm",ylim = (-0.01, 0.51)
                )
ax = axs[1,2]
fig,ax = get_subplot(fig,ax,df_energy["t"], df_energy["total_energy_normalized"],
                 xlabel = "t", ylabel = "total_norm",ylim = (-0.01, 2.01)
                )

os.makedirs("plots/", exist_ok = True)
plt.savefig("plots/energy_v_time.pdf", bbox_inches = 'tight')


#####################################################

############
# Plotting #
############

database = pd.read_csv("database.csv")
df = database

max_N_iter = max(df['N_iter'])
min_N_iter = min(df['N_iter'])
N_iters = np.sort(np.unique(df['N_iter']))

fig,ax = plt.subplots(figsize = (13,6))

norm_max = 2
norm_min = 0
norm = plt.Normalize(norm_min, norm_max)
cmap = 'viridis'
lw = 5

skipby = 10


#plot vitelline membrane
thetas = np.linspace(0, np.pi, 1000)
vit_x = np.cos(thetas)
vit_y = 0.4*np.sin(thetas)
ax.plot(1.02*vit_x, 1.02*vit_y, lw = 5, alpha = 1, color = 'gray')

for i in range(len(N_iters)):
    if i%skipby != 0:
        if i != len(N_iters)-1:
            continue
    #get timepoint
    N_iter = N_iters[i]
    df_timepoint = df[df['N_iter'] == N_iter]
    alpha = (N_iter - min_N_iter)/(max_N_iter - min_N_iter)
    #get points
    points_1 = np.array(df_timepoint[['x', 'y']]).reshape(-1, 1, 2)[:-1]
    points_2 = np.array(df_timepoint[['x', 'y']]).reshape(-1, 1, 2)[1:]
    #create a collection of lines
    segments_demo = np.concatenate([points_1, points_2], axis = 1)
    #value by which to color lines
    dydx = df_timepoint["d_r_mod"][:-1]/df_timepoint["dr0"][:-1]
    colors = np.where(df_timepoint["mitotic_domain_id"] == -1, "purple", "green")
    linewidths = np.where(df_timepoint["mitotic_domain_id"] == -1, 2, 4)
    #make line collection
    #lc = LineCollection(segments_demo, array = dydx,cmap=cmap, norm=norm,alpha = alpha,lw = lw,)
    lc = LineCollection(segments_demo, colors =colors,alpha = alpha,linewidths=linewidths)
    line = ax.add_collection(lc)

ax.axis('equal')
ax.set_ylim(-0.1, 0.5), #ax.set_xticks([])
ax.set_xlim(-1.1, 1.1), #ax.set_yticks([])

#cbar = fig.colorbar(line, ax=ax, alpha = 1, 
                    #label = 'current / preferred length',
#                   )
#cbar.solids.set_edgecolor("face")
#cbar.ax.tick_params(labelsize=15) 
#cbar.set_label(label='current / preferred length', size='large',fontsize = 30,
               #pad = 10,weight='bold'
#               )


os.makedirs("plots/", exist_ok = True)
plt.savefig("plots/curve.pdf", bbox_inches = 'tight')





#################
# Fold analysis #
#################

try:
    f = open("../fold_threshold.txt", "r")
    fold_threshold = float(f.read())
    print("reading fold threshold from text file")
except:
    fold_threshold = 0.0227
    print("could not file file to read fold threshold")

#for each timepoint, calculate the number of folds, CF width of influence, 


def dist_to_ellipse(row):
    closest_point = np.array(get_closest_ellipse_pt(1, 0.4, (row.x, row.y)))
    dist = np.linalg.norm(closest_point - np.array([row.x, row.y]))
    return(dist)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def get_fold_analysis(row, t = None, colName_suffix = "", fold_threshold = 0.02266):
    
    #get the timepoint
    df = database[database["t"] == t]
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
    df["fold_inside_CF_bool"] = 0
    CF_area_of_influence_ref = 0
    CF_area_of_influence_current = 0
    CF_depth = 0
    CF_ind = df[df["ko"] != 0].index
    if len(CF_ind) > 0:
        max_CF_arclength_frac = max(df.loc[CF_ind,"arclength_frac"])
        min_CF_arclength_frac = min(df.loc[CF_ind,"arclength_frac"])
        df["fold_inside_CF_bool"] = df.groupby("fold_id")["fold_id"].transform(lambda x : (float((max(df.loc[x.index,"arclength_frac"]) > max_CF_arclength_frac) & (min(df.loc[x.index,"arclength_frac"]) < min_CF_arclength_frac)))*df.loc[x.index, "fold_mask"])    
        fold_id_inside_CF = np.sort(np.unique(df["fold_id"]*df["fold_inside_CF_bool"]))[-1] #sort and then take last value to ensure fold_id is caught (if no fold in CF then 0 is returned)
        if fold_id_inside_CF>0:
            CF_area_of_influence_ref = df[df["fold_id"] == fold_id_inside_CF]["fold_width_ref"].values[0]
            CF_area_of_influence_current = df[df["fold_id"] == fold_id_inside_CF]["fold_width_current"].values[0]
            CF_depth = max(df[df["fold_id"] == fold_id_inside_CF]["dist_to_ellipse"])
        
    row["nb_folds" + colName_suffix] = nb_folds
    row["avg_depth" + colName_suffix] = np.mean(df["dist_to_ellipse"])
    row["max_depth" + colName_suffix] = np.max(df["dist_to_ellipse"])
    row["CF_area_of_influence_ref" + colName_suffix] = CF_area_of_influence_ref
    row["CF_area_of_influence_current" + colName_suffix] = CF_area_of_influence_current
    row["CF_depth" + colName_suffix] = CF_depth
    
    return(row)

df_energy = df_energy.apply(lambda row: get_fold_analysis(row, t = row["t"], colName_suffix = "", fold_threshold = fold_threshold) ,axis = 1)

df_energy.to_csv("summed_prop_per_timepoint.csv", index = False)

#####
# Plotting 
#####

fig,ax = plt.subplots()
ax.plot(df_energy["t"], df_energy["nb_folds"], label = "Nb folds", marker = "o")
ax.plot(df_energy["t"], df_energy["CF_area_of_influence_ref"], label = "CF area of influence", marker = "o")
ax.plot(df_energy["t"], df_energy["CF_depth"], label = "CF depth", marker = "o")
ax.legend(loc = "left top", fontsize = 12)
os.makedirs("plots/", exist_ok = True)
plt.savefig("plots/folds_CF_analysis.pdf", bbox_inches = 'tight')

