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





