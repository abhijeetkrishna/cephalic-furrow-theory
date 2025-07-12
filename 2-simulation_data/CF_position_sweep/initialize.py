from CF_methods import *

#Taking arguments from command line
map_index_dest=sys.argv[1]
task_id=int(sys.argv[2])

#print('task_id : ' + str(task_id))

##############
# Parameters #
##############

#defining a variable dictionary
var_dict={
    'K' : 1e-5,
    'L' : 1.0,
    'R_vit' : 1,
    'R': 1,
    'R_vert_frac': 0.4,
    #'k_vit' : 10,
    't_final' : 25,
    'dt' : 1e-5, #0.0001,
    'seed' : 0,
    #'uniform_growth' : 1,
    #'mitotic_growth' : 1.5,
    'int_curvature_MD' : 0.0,
    'int_curvature_CF' : 20.0,
    'int_curvature_non_MD' : 0.0,
    #'dt_after_noise' : None,
    #'dt_upp_thr' : 10**(-3),
    #'dt_low_thr': 10**(-5),
    #'MD_separation' : 5,
    'vitelline_shape' : 'line', 
    #'MD_growth_const_l0': True, 
    #'adjust_growth_rate': False,
    'germ_band_push' : 0.0,
    'mitotic_domains_bool' : False,
    'noise_iterations' : 1000,
    'noise_scale' : 1e-5,
    'output_intervals' : 600000, #5000,
    'Nballs' : 100,
    'number_of_mitotic_domains' : 0,
    #'K_cf': None,
    'dt_after_noise':None,
    'no_noise_mitotic_domain':False,
    'K_md_by_K' : 1,#0.14,#1,
    'add_CF_bool' : 0.0, #0 - False
    't_CF' : 0.0,
    'add_MD_bool' : 1.0, #1 - True
    't_MD' : 0.0,
    'CF_position' : 0.37,
    'CF_width' : 0.02,
    }

#reading the map_index csv
map_index=pd.read_csv(map_index_dest)
input_var=list(map_index.columns) #variables for which values are being imported
for i in range(len(input_var)):
    if input_var[i] == 'folder_name':
        dirname = map_index.loc[task_id,str(input_var[i])]
        continue
    var_dict[input_var[i]]=map_index.loc[task_id,str(input_var[i])]

#defining variables 
#from dictionary
for key,val in var_dict.items():
    exec(key + '=val')



continued_simulation =  False #this bool tells us whether to start the simulation from start or continue from some point
#if os.path.exists( dirname + 'database.pickle'):
#    print('database already exists - need to continue simulation')
#    # if database.pickle exists
#    continued_simulation = True

if not os.path.exists(dirname):
    os.mkdir(dirname)

#anneal_checkpoint = 2/t_final

anneal_checkpoint = noise_iterations*dt/t_final #this ensures that noise is added only for 10 timesteps

#output_intervals = 1/dt #when total time was 10000
#output_intervals = 1
#noise_scale = R_vit/10000


#dynamic_dt = False


bend_energy = True
stretch_energy = True
noise = True
#debug = False
#vitelline = True
#print('vitelline : '+str(vitelline))
#elastic_vitelline = False
#print('elastic_vitelline : '+str(elastic_vitelline))
#mitotic_domains_bool = True
print('mitotic_domains_bool : ' + str(mitotic_domains_bool))
#mitotic_domains = []

#if mitotic_domains_bool:
#    #generate mitotic domains
#    m1 = np.arange(13, 24)
#    m2 = np.arange(30, 45)
#    #m2 = np.arange(30, 40)
#    m3 = np.arange(55, 70)
#    #m3 = np.arange(45, 60)
#    m4 = np.arange(80, 95)
#    #m4 = np.arange(80, 95)
#    mitotic_domains = np.hstack((m1, m2, m3, m4))


###################
# Initialize ring #
###################
#balls, springs = get_tissue(Nballs, K, L, R=R, R_vert_frac = R_vert_frac,  vitelline_shape = vitelline_shape, close_loop=False, 
#                                mitotic_domains_bool=mitotic_domains_bool, mitotic_domains=mitotic_domains, 
#                               MD_growth_const_l0 = MD_growth_const_lo, 
#                                mitotic_growth = mitotic_growth)

#we don
balls, springs = get_tissue(Nballs, K, L, K_md = K_md_by_K*K, R=R, R_vert_frac = R_vert_frac,  vitelline_shape = vitelline_shape, #close_loop=False,  
                            germ_band_push = germ_band_push,
                            )

if continued_simulation:

    # read the last position of balls
    df_temp = pickle.load(open(glob.glob(dirname + 'database.pickle')[0], "rb"))
    bool_ind = df_temp.N_iter == np.max(df_temp.N_iter)
    df_temp = df_temp[bool_ind].reset_index()

    # update position of balls
    balls[['x', 'y']] = df_temp[['x', 'y']]

    # update springs  - not required

    # rename the database 
    name_found = False
    i = 1
    while not(name_found):

        if os.path.exists(dirname + 'database' + str(i) + '.pickle'):

            i = i + 1

        else:

            #rename the file with name = 'database' + str(i) + '.pickle'
            os.rename(dirname + 'database.pickle', dirname + 'database' + str(i) + '.pickle')

            name_found = True

    print('renamed database.pickle to database' + str(i) + '.pickle')


    noise = False
    dt_after_noise = None
    

#For mitotic domains
#balls["ko"] = int_curvature_non_MD
#balls["ko"] = np.where(balls.ID.isin(mitotic_domains), int_curvature_MD, int_curvature_non_MD)
#For cephalic furrow
#balls["ko"][[ 24, 25, 26, 74, 75, 76]] = int_curvature_CF
#balls["ko_i-1"] = np.array(balls["ko"][(np.array(balls.ID)+len(balls)-1)%len(balls)])
#balls["ko_i+1"] = np.array(balls["ko"][(np.array(balls.ID)+1)%len(balls)])

##########
# export #
##########

balls.to_csv(dirname + 'init_balls.csv', index = False)
balls.to_csv(dirname + 'reference_balls.csv', index = False)
springs.to_csv(dirname + 'init_springs.csv', index = False)

#export simulations
sim_params = pd.DataFrame(columns = ['variable', 'value'])

sim_params.loc[len(sim_params),['variable', 'value']] = ['seed', seed]
sim_params.loc[len(sim_params),['variable', 'value']] = ['t_final', t_final]
sim_params.loc[len(sim_params),['variable', 'value']] = ['t_CF', t_CF]
sim_params.loc[len(sim_params),['variable', 'value']] = ['add_CF_bool', add_CF_bool]
sim_params.loc[len(sim_params),['variable', 'value']] = ['CF_position', CF_position]
sim_params.loc[len(sim_params),['variable', 'value']] = ['CF_width', CF_width]
sim_params.loc[len(sim_params),['variable', 'value']] = ['int_curvature_CF', int_curvature_CF]
sim_params.loc[len(sim_params),['variable', 'value']] = ['t_MD', t_MD]
sim_params.loc[len(sim_params),['variable', 'value']] = ['add_MD_bool', add_MD_bool]
sim_params.loc[len(sim_params),['variable', 'value']] = ['int_curvature_MD', int_curvature_MD]
sim_params.loc[len(sim_params),['variable', 'value']] = ['dt', dt]
sim_params.loc[len(sim_params),['variable', 'value']] = ['vitelline_shape', vitelline_shape]
sim_params.loc[len(sim_params),['variable', 'value']] = ['R', R]
sim_params.loc[len(sim_params),['variable', 'value']] = ['R_vert_frac', R_vert_frac]
sim_params.loc[len(sim_params),['variable', 'value']] = ['L', L]
sim_params.loc[len(sim_params),['variable', 'value']] = ['K', K]
sim_params.loc[len(sim_params),['variable', 'value']] = ['K_md', K_md_by_K*K]
sim_params.loc[len(sim_params),['variable', 'value']] = ['R_vit', R_vit]
sim_params.loc[len(sim_params),['variable', 'value']] = ['bend_energy', bend_energy]
sim_params.loc[len(sim_params),['variable', 'value']] = ['stretch_energy', stretch_energy]
sim_params.loc[len(sim_params),['variable', 'value']] = ['output_intervals', output_intervals]
sim_params.loc[len(sim_params),['variable', 'value']] = ['noise', noise]
sim_params.loc[len(sim_params),['variable', 'value']] = ['noise_scale', noise_scale]
sim_params.loc[len(sim_params),['variable', 'value']] = ['anneal_checkpoint', anneal_checkpoint]
sim_params.loc[len(sim_params),['variable', 'value']] = ['mitotic_domains_bool', mitotic_domains_bool]
sim_params.loc[len(sim_params),['variable', 'value']] = ['output_dir', dirname]
#sim_params.loc[len(sim_params),['variable', 'value']] = ['dt_after_noise', dt_after_noise]
#simulations_param.loc[len(simulations_param),'no_noise_mitotic_domain'] = no_noise_mitotic_domain

sim_params.to_csv(dirname + 'sim_params.csv', index = False)


##############
# simulation #
##############

#np.random.seed(seed)

#balls_timeseries, final_len = energy_minimize(balls,springs,t_final=t_final, dt = dt, 
#                                               vitelline_shape = vitelline_shape, 
#                                               R = R, R_vert_frac = R_vert_frac, #coords_vit = coords_vit,
#                                               L=L, K= K, 
#                                               R_vit = R_vit, #elastic_vitelline=elastic_vitelline,vitelline = vitelline, 
#                                               bend_energy=bend_energy, stretch_energy = stretch_energy,
#                                               #germ_band_pushing = germ_band_pushing, F_germ_band = F_germ_band, push_thresh = push_thresh, 
#                                               output_intervals = output_intervals, 
#                                               noise = noise, noise_scale=noise_scale, anneal_checkpoint = anneal_checkpoint,
#                                               mitotic_domains_bool=mitotic_domains_bool,
#                                               output_dir=dirname,dt_after_noise=dt_after_noise,
#                                               #dt_upp_thr = dt_upp_thr, dt_low_thr = dt_low_thr, #variable_dt_start = variable_dt_start,
#                                               no_noise_mitotic_domain = no_noise_mitotic_domain,
#                                             )

#filehandler = open(dirname+'database.pickle', 'wb') 
#pickle.dump(balls_timeseries, filehandler)






