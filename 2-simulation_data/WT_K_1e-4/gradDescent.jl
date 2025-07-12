#using Pkg
#Pkg.add("DataFrames")
using DataFrames
#Pkg.add("CSV")
using CSV
#Pkg.add("PyCall")
using PyCall
#Pkg.add("Random")
using Random
#Pkg.add("NLsolve")
using NLsolve

function add_CF(balls; int_curvature_CF = 10.0)

    CF_id = balls[0.36 .< balls[:,"arclength_frac"] .< 0.38, "ID"]
    CF_id = CF_id .+ 1 #for Julia since row number starts from 1
    balls[CF_id, "ko"] .= int_curvature_CF
    #save the dataframe before running forward
    return balls
end

function add_MD(balls; int_curvature_MD = 0, K_md = 1)

    #assuming number of mitotic_domains = 4
    md_locations = [(0.1, 0.15), (0.2,0.25), (0.3,0.35), (0.4, 0.45)]

    row_pointer = 0#-1
    in_md = false
    md_pointer = 1 #0
    md_location = md_locations[md_pointer]

    while md_pointer <= length(md_locations)

        row_pointer += 1
        if row_pointer > nrow(balls)
            break
        end

        arclength_frac = balls[row_pointer, "arclength_frac"]

        md_location = md_locations[md_pointer]
        if arclength_frac >= md_location[1]
            if arclength_frac <= md_location[2]
                in_md = true
            else
                #we have come outside the mitotic domain
                in_md = false
                md_pointer += 1
            end
        end

        if in_md
            #particles before and after belong inside the mitotic domain
            balls[row_pointer-1, ["mitotic_domain_id", "mitotic_domain", "K", "ko"]] = [md_pointer - 1, 1.0, K_md, int_curvature_MD] #-1 to keep it consistent with the python code
            balls[row_pointer, ["mitotic_domain_id", "mitotic_domain", "K", "ko"]] = [md_pointer - 1, 1.0, K_md, int_curvature_MD] #-1 to keep it consistent with the python code

            #computing coordinates of new particle
            x_new = (balls[row_pointer-1,"x"] + balls[row_pointer,"x"])/2
            y_new = (balls[row_pointer-1,"y"] + balls[row_pointer,"y"])/2
            #[x_new, y_new] =  [(balls[row_pointer-1,'x'] + balls[row_pointer,'x'])/2, (balls[row_pointer-1,'y'] + balls[row_pointer,'y'])/2]
            arclength_new = balls[row_pointer-1, "arclength"] + sqrt( (balls[row_pointer-1,"x"] - x_new)^2 + (balls[row_pointer-1,"y"] - y_new)^2 )
            arclength_frac_new = arclength_new/balls[nrow(balls), "arclength"]

            #preparing row for new particle
            row = balls[[row_pointer], :] #get a copy of the previous row
            row[1,["x","y","arclength","arclength_frac"]] = [x_new, y_new, arclength_new, arclength_frac_new]
            #add row at row_pointer location
            balls = [view(balls, 1:(row_pointer-1), :); DataFrame(names(balls) .=> Array(row[1,:])); view(balls, row_pointer:nrow(balls), :)]

            #increase row pointer
            row_pointer += 1
        end
    end
        
    balls[:,"ID"] = 0:(nrow(balls)-1)
    balls[!,"ID"] = convert.(Integer,balls[!,"ID"])
    
    balls[1:end-1, "K+1"] = balls[2:end, "K"]
    balls[2:end, "K-1"] = balls[1:end-1, "K"]    
    balls[1:end-1, "ko_i+1"] = balls[2:end, "ko"]
    balls[2:end, "ko_i-1"] = balls[1:end-1, "ko"]

    return balls

end

function energy_minimize(balls; #compulsory and optional arguments must be seperated by semi-colon
                    springs = nothing, seed = 0, 
                    t_final = 100, 
                    t_CF = 0, t_MD = 0, 
                    add_CF_bool = false, add_MD_bool = false, 
                    int_curvature_CF = 0, int_curvature_MD = 0, 
                    dt = 0.01, dynamic_dt = false, 
                    vitelline_shape = "line",
                    R = 1, R_vert_frac = 0, coords_vit = nothing,
                    K = 1, L = 0.1, K_md = 1,
                    k_vit = 1, tau = nothing, tol = 1e-12, dynamic_t_final = true,
                    debug = false, output_intervals = 10,
                    R_vit = 1, vitelline = false, elastic_vitelline = false,
                    bend_energy = true, stretch_energy = true,
                    noise = false, noise_scale = 1, anneal_checkpoint = nothing, no_noise_mitotic_domain = false,
                    growth = false, desired_growth = 1.1, growth_intervals = 0.5, growth_rate = 1.005,
                    mitotic_domains_growth = false, mitotic_domains_bool = false, int_curve_checkpoint = nothing,
                    output_dir = nothing,
                    dt_after_noise = nothing, dt_upp_thr = 10^(-2), dt_low_thr = 10^(-4), variable_dt_start = 0.001,
                    )

    #add_CF and add_MD should be true if you need to add CF or MD after t_CF or t_MD time

    #setting the seed
    rng = MersenneTwister(seed) #random number generator

    first_anneal_bool = true

    if tau  == nothing
        tau = 1/L
    end
    
    R1 = R
    R2 = R_vert_frac
    R_vit = R #this is abitrary - this is a constant value in the dynamics eqation and should not matter
    
    #balls = balls.copy(deep=True)
    #springs =springs.copy(deep=True)
    balls_timeseries = DataFrame()
    
    t = 0.0
    times = []
    N_iter = -1
    
    t_anneal_checkpoint = anneal_checkpoint*t_final
    #for annealing, we will stop the noise and increase dt little by little
    
    #cols = []
    #balls[:,cols] .= 0


    #balls[:, "rowid"] = balls[:, "ID"] .+ 1
    col_to_set_zero = ["x_i-2", "y_i-2", "x_i-1", "y_i-1", "x_i+1", "y_i+1", "x_i+2", "y_i+2", "curvature", "curvature_i-1", "curvature_i+1", "d_r_mod", "d_r_mod_i-1","d_r_mod_i+1","d_r_mod_alternate","d_r_mod_alternate_i-1","d_r_mod_alternate_i+1", "d_curvature_x", "d_curvature_y", "d_curvature_i-1_x","d_curvature_i+1_x","d_curvature_i-1_y","d_curvature_i+1_y","d_energy_x", "d_energy_y", "noise_x", "noise_y"]
    for col in col_to_set_zero
        balls[:,col] .= 0.0 #it is important to write 0.0 to set the type of column as float
    end
    cols_to_write = ["ID", "mitotic_domain", "mitotic_domain_id","dr0", "K", "L+1", "L-1","ko","arclength", "arclength_frac", "t", "N_iter", "x", "y", "curvature", "d_r_mod", "dx_t", "dy_t", "d_r_mod_alternate","noise_x", "noise_y"]

    while t<t_final

        if add_CF_bool && t>=t_CF
            println("Adding CF")
            CSV.write(string(output_dir,"before_adding_CF.csv"), balls)
            balls = add_CF(balls, int_curvature_CF = int_curvature_CF)
            add_CF_bool = false #we don't add CF anymore for rest of simulation
            CSV.write(string(output_dir,"after_adding_CF.csv"), balls)
            CSV.write(string(output_dir,"reference_balls.csv"), balls)
        end

        if add_MD_bool && t>=t_MD
            println("Adding MD")
            CSV.write(string(output_dir,"before_adding_MD.csv"), balls)
            balls = add_MD(balls, int_curvature_MD = int_curvature_MD, K_md = K_md)
            add_MD_bool = false #we don't add MD anymore for rest of simulation
            CSV.write(string(output_dir,"after_adding_MD.csv"), balls)
            CSV.write(string(output_dir,"reference_balls.csv"), balls)
        end

        balls[:,"d_energy_x"] .= 0.0
        balls[:,"d_energy_y"] .= 0.0

        #println(string("time : ", string(t)))

        N_iter = N_iter + 1
        energy = 0
        
        balls[:,"t"] .= t
        balls[:,"N_iter"] .= N_iter
                
        #col_to_set_zero = ["x_i-2", "y_i-2", "x_i-1", "y_i-1", "x_i+1", "y_i+1", "x_i+2", "y_i+2", "curvature", "curvature_i-1", "curvature_i+1", "d_r_mod", "d_r_mod_i-1","d_r_mod_i+1","d_r_mod_alternate","d_r_mod_alternate_i-1","d_r_mod_alternate_i+1", "d_curvature_x", "d_curvature_y", "d_curvature_i-1_x","d_curvature_i+1_x","d_curvature_i-1_y","d_curvature_i+1_y","d_energy_x", "d_energy_y", "noise_x", "noise_y"]
        #balls[!,col_to_set_zero] .= 0
        #for col in col_to_set_zero
        #    balls[:,col] .= 0
        #end
        
        #balls[3:end, "x_i-2"] = balls[balls[3:end,"rowid"] .- 2,"x"]     #np.array(balls["x"][(np.array(balls.ID[2:])-2)])
        #balls[3:end, "y_i-2"] = balls[balls[3:end,"rowid"] .- 2,"y"]     #np.array(balls["y"][(np.array(balls.ID[2:])-2)])
        #balls[2:end, "x_i-1"] = balls[balls[2:end,"rowid"] .- 1,"x"]     #np.array(balls["x"][(np.array(balls.ID[1:])-1)])
        #balls[2:end, "y_i-1"] = balls[balls[2:end,"rowid"] .- 1,"y"]     #np.array(balls["y"][(np.array(balls.ID[1:])-1)])
        #balls[1:end-1, "x_i+1"] = balls[balls[1:end-1,"rowid"] .+ 1,"x"] #np.array(balls["x"][(np.array(balls.ID)[:-1]+1)])
        #balls[1:end-1, "y_i+1"] = balls[balls[1:end-1,"rowid"] .+ 1,"y"] #np.array(balls["y"][(np.array(balls.ID)[:-1]+1)])
        #balls[1:end-2, "x_i+2"] = balls[balls[1:end-2,"rowid"] .+ 2,"x"] #np.array(balls["x"][(np.array(balls.ID)[:-2]+2)])
        #balls[1:end-2, "y_i+2"] = balls[balls[1:end-2,"rowid"] .+ 2,"y"] #np.array(balls["y"][(np.array(balls.ID)[:-2]+2)])

        balls[3:end, "x_i-2"] = balls[1:(end-2),"x"]     #np.array(balls["x"][(np.array(balls.ID[2:])-2)])
        balls[3:end, "y_i-2"] = balls[1:(end-2),"y"]     #np.array(balls["y"][(np.array(balls.ID[2:])-2)])
        balls[2:end, "x_i-1"] = balls[1:(end-1),"x"]     #np.array(balls["x"][(np.array(balls.ID[1:])-1)])
        balls[2:end, "y_i-1"] = balls[1:(end-1),"y"]     #np.array(balls["y"][(np.array(balls.ID[1:])-1)])
        balls[1:end-1, "x_i+1"] = balls[2:end,"x"] #np.array(balls["x"][(np.array(balls.ID)[:-1]+1)])
        balls[1:end-1, "y_i+1"] = balls[2:end,"y"] #np.array(balls["y"][(np.array(balls.ID)[:-1]+1)])
        balls[1:end-2, "x_i+2"] = balls[3:end,"x"] #np.array(balls["x"][(np.array(balls.ID)[:-2]+2)])
        balls[1:end-2, "y_i+2"] = balls[3:end,"y"] #np.array(balls["y"][(np.array(balls.ID)[:-2]+2)])
        
        balls[:,"r"] = (balls[:,"x"].^2 + balls[:,"y"].^2).^0.5 #np.array(np.sqrt((balls["x"])**2+(balls["y"])**2))
        
        # Measure curvature
        
        balls[2:end-1, "curvature"] = ( 4* ((balls[2:end-1,"x_i+1"] - balls[2:end-1,"x_i-1"]).*(balls[2:end-1,"y_i+1"] - 2*balls[2:end-1,"y"] + balls[2:end-1,"y_i-1"]) - (balls[2:end-1,"y_i+1"] - balls[2:end-1,"y_i-1"]).*(balls[2:end-1,"x_i+1"] - 2*balls[2:end-1,"x"] + balls[2:end-1,"x_i-1"]))./(((balls[2:end-1,"x_i+1"] - balls[2:end-1,"x_i-1"]).^2 + (balls[2:end-1,"y_i+1"] - balls[2:end-1,"y_i-1"]).^2).^(3/2)))
        #balls["curvature"][1:-1] += np.array((4* ((balls["x_i+1"] - balls["x_i-1"])*(balls["y_i+1"] - 2*balls["y"] + balls["y_i-1"]) - (balls["y_i+1"] - balls["y_i-1"])*(balls["x_i+1"] - 2*balls["x"] + balls["x_i-1"]))/(((balls["x_i+1"] - balls["x_i-1"])**2 + (balls["y_i+1"] - balls["y_i-1"])**2)**(3/2))))[1:-1]
        balls[2:end,"curvature_i-1"] = balls[1:(end-1),"curvature"]
        #balls["curvature_i-1"][1:] += np.array(balls["curvature"][(np.array(balls.ID[1:])-1)])
        balls[1:end-1,"curvature_i+1"] = balls[2:end,"curvature"]
        #balls["curvature_i+1"][:-1] += np.array(balls["curvature"][(np.array(balls.ID[:-1])+1)])
        balls[1:end-1,"d_r_mod"] = ((balls[1:end-1,"x_i+1"]-balls[1:end-1,"x"]).^2 + (balls[1:end-1,"y_i+1"]-balls[1:end-1,"y"]).^2).^0.5
        #balls["d_r_mod"][:-1] += np.array(((balls["x_i+1"]-balls["x"])**2 + (balls["y_i+1"]-balls["y"])**2)**0.5)[:-1]

        # non boundary nodes
        #balls[["d_r_mod_i-1","d_r_mod_i+1","d_r_mod_alternate","d_r_mod_alternate_i-1","d_r_mod_alternate_i+1"]] = 0
        #idx = ~np.isin(np.arange(len(balls)), [0,1,len(balls)-2, len(balls)-1])
        idx = 3:(nrow(balls)-2)
        
        balls[idx,"d_r_mod_i-1"] = ((balls[idx,"x_i-1"]-balls[idx,"x"]).^2 + (balls[idx,"y_i-1"]-balls[idx,"y"]).^2).^0.5
        balls[idx,"d_r_mod_i+1"] = ((balls[idx,"x_i+2"]-balls[idx,"x_i+1"]).^2 + (balls[idx,"y_i+2"]-balls[idx,"y_i+1"]).^2).^0.5
        balls[idx,"d_r_mod_alternate"] = (((balls[idx,"x_i+1"]-balls[idx,"x_i-1"]).^2 + (balls[idx,"y_i+1"]-balls[idx,"y_i-1"]).^2).^0.5)./2
        balls[idx,"d_r_mod_alternate_i-1"] = (((balls[idx,"x"]-balls[idx,"x_i-2"]).^2 + (balls[idx,"y"]-balls[idx,"y_i-2"]).^2).^0.5)./2
        balls[idx,"d_r_mod_alternate_i+1"] = (((balls[idx,"x_i+2"]-balls[idx,"x"]).^2 + (balls[idx,"y_i+2"]-balls[idx,"y"]).^2).^0.5)./2
        
        #compute curvature gradient
        #balls[["d_curvature_x", "d_curvature_y"]] = 0
        balls[idx,"d_curvature_x"] = (balls[idx,"y_i+1"] - balls[idx,"y_i-1"])./(balls[idx,"d_r_mod_alternate"].^3)
        balls[idx,"d_curvature_y"] = -(balls[idx,"x_i+1"] - balls[idx,"x_i-1"])./(balls[idx,"d_r_mod_alternate"].^3)

        # adding curvature gradient of neighbours
        #balls[["d_curvature_i-1_x","d_curvature_i+1_x","d_curvature_i-1_y","d_curvature_i+1_y"]] = 0
        balls[idx,"d_curvature_i-1_x"] = -(balls[idx,"y_i-1"] - balls[idx,"y_i-2"])./((balls[idx,"d_r_mod_alternate_i-1"]).^3) - 3*(balls[idx,"x"] - balls[idx,"x_i-2"]).*balls[idx,"curvature_i-1"]./((2*balls[idx,"d_r_mod_alternate_i-1"]).^2)
        balls[idx,"d_curvature_i+1_x"] = -(balls[idx,"y_i+2"] - balls[idx,"y_i+1"])./((balls[idx,"d_r_mod_alternate_i+1"]).^3) + 3*(balls[idx,"x_i+2"] - balls[idx,"x"]).*balls[idx,"curvature_i+1"]./((2*balls[idx,"d_r_mod_alternate_i+1"]).^2)
        balls[idx,"d_curvature_i-1_y"] = (balls[idx,"x_i-1"] - balls[idx,"x_i-2"])./((balls[idx,"d_r_mod_alternate_i-1"]).^3) - 3*(balls[idx,"y"] - balls[idx,"y_i-2"]).*balls[idx,"curvature_i-1"]./((2*balls[idx,"d_r_mod_alternate_i-1"]).^2)
        balls[idx,"d_curvature_i+1_y"] = (balls[idx,"x_i+2"] - balls[idx,"x_i+1"])./((balls[idx,"d_r_mod_alternate_i+1"]).^3) + 3*(balls[idx,"y_i+2"] - balls[idx,"y"]).*balls[idx,"curvature_i+1"]./((2*balls[idx,"d_r_mod_alternate_i+1"]).^2)
        
        #compute energy gradient 
        #balls[["d_energy_x", "d_energy_y"]] = 0
        
        if stretch_energy
          
            balls[idx,"d_energy_x"] += balls[idx,"L+1"].*(balls[idx,"d_r_mod"]./balls[idx,"dr0"] .- 1).*(balls[idx,"x"] - balls[idx,"x_i+1"])./(balls[idx,"d_r_mod"]) + 
                                             balls[idx,"L-1"].*(balls[idx,"d_r_mod_i-1"]./balls[idx,"dr0_i-1"] .- 1).*(balls[idx,"x"] - balls[idx,"x_i-1"])./(balls[idx,"d_r_mod_i-1"])
        
            balls[idx,"d_energy_y"] += balls[idx,"L+1"].*(balls[idx,"d_r_mod"]./balls[idx,"dr0"] .- 1).*(balls[idx,"y"] - balls[idx,"y_i+1"])./(balls[idx,"d_r_mod"]) + 
                                             balls[idx,"L-1"].*(balls[idx,"d_r_mod_i-1"]./balls[idx,"dr0_i-1"] .- 1).*(balls[idx,"y"] - balls[idx,"y_i-1"])./(balls[idx,"d_r_mod_i-1"])
        end


        #adding force due to bending energy
        if bend_energy
            
            balls[idx,"d_energy_x"] += 2*balls[idx,"K-1"].*(balls[idx,"curvature_i-1"]-balls[idx,"ko_i-1"]).*balls[idx,"d_curvature_i-1_x"].*balls[idx,"dr0_i-1"] + 
                                            #balls[idx,"K-1"].*((balls[idx,"curvature_i-1"]-balls[idx,"ko_i-1"]).^2).*(balls[idx,"x"] - balls[idx,"x_i-2"])./(4*balls[idx,"d_r_mod_alternate_i-1"]) + 
                                            2*balls[idx,"K"].*(balls[idx,"curvature"] - balls[idx,"ko"]).*balls[idx,"d_curvature_x"].*balls[idx,"dr0"] + 
                                            2*balls[idx,"K+1"].*(balls[idx,"curvature_i+1"] - balls[idx,"ko_i+1"]).*balls[idx,"d_curvature_i+1_x"].*balls[idx,"dr0_i+1"] #+ 
                                            #balls[idx,"K+1"].*((balls[idx,"curvature_i+1"] - balls[idx,"ko_i+1"]).^2).*(balls[idx,"x"] - balls[idx,"x_i+2"])./(4*balls[idx,"d_r_mod_alternate_i+1"])

            balls[idx,"d_energy_y"] += 2*balls[idx,"K-1"].*(balls[idx,"curvature_i-1"]-balls[idx,"ko_i-1"]).*balls[idx,"d_curvature_i-1_y"].*balls[idx,"dr0_i-1"] + 
                                            #balls[idx,"K-1"].*((balls[idx,"curvature_i-1"]-balls[idx,"ko_i-1"]).^2).*(balls[idx,"y"] - balls[idx,"y_i-2"])./(4*balls[idx,"d_r_mod_alternate_i-1"]) + 
                                            2*balls[idx,"K"].*(balls[idx,"curvature"] - balls[idx,"ko"]).*balls[idx,"d_curvature_y"].*balls[idx,"dr0"] + 
                                            2*balls[idx,"K+1"].*(balls[idx,"curvature_i+1"] - balls[idx,"ko_i+1"]).*balls[idx,"d_curvature_i+1_y"].*balls[idx,"dr0_i+1"] #+ 
                                            #balls[idx,"K+1"].*((balls[idx,"curvature_i+1"] - balls[idx,"ko_i+1"]).^2).*(balls[idx,"y"] - balls[idx,"y_i+2"])./(4*balls[idx,"d_r_mod_alternate_i+1"])
        end

        #compute displacement
        balls[:,"dx_t"] = -(balls[:,"d_energy_x"]).*((R_vit)/(tau * L))
        balls[:,"dy_t"] = -(balls[:,"d_energy_y"]).*((R_vit)/(tau * L)) #tau is set so that tau*L is 1

        balls[:,"x_t+1"] = balls[:,"x"] + dt*balls[:,"dx_t"] 
        balls[:,"y_t+1"] = balls[:,"y"] + dt*balls[:,"dy_t"]

        balls[idx,"noise_x"] = noise_scale*randn(rng, Float64, size(idx)) #np.random.normal(loc = 0, scale = noise_scale, size = len(balls))[idx]
        balls[idx,"noise_y"] = noise_scale*randn(rng, Float64, size(idx)) #np.random.normal(loc = 0, scale = noise_scale, size = len(balls))[idx]

        balls[:,"x_t+1"] = balls[:,"x_t+1"] + balls[:,"noise_x"]
        balls[:,"y_t+1"] = balls[:,"y_t+1"] + balls[:,"noise_y"]
            
        #adding vitelline repulsion
        
        if R_vert_frac == 0
            # find points that are outside vitelline
            #balls_outside = balls[balls['y_t+1']>0]
            balls_outside_id = balls[balls[:,"y_t+1"] .>0, "ID"] #assuming that ID matches with the row number
            balls_outside_id = balls_outside_id .+ 1 #for Julia since row number starts from 1

            # intersection of two straight lines
            #for index, row in balls_outside.iterrows():
            for index in balls_outside_id
                
                row = balls[index,:] #here we add 1 to index because row number starts from 1 while ID starts from 0
                P1 = [row["x"], row["y"]]
                P2 = [row["x_t+1"], row["y_t+1"]]
                
                slope = (P2[2] - P1[2])/(P2[1] - P1[1]) #slope
                intercept = P1[2] - slope*P1[1]

                #[x_guess, y_guess] = P1
                
                function f!(F,x)
                    #return [line equation, line equation]
                    F[1] = x[2] - slope*x[1] - intercept
                    F[2] = x[2]
                    #return( (x[1] - slope*x[0] - intercept , x[1]) )
                end
                
                #root = fsolve(func, P1)
                root = nlsolve(f!, P1).zero
                
                balls[index, "x_t+1"] = root[1]
                balls[index, "y_t+1"] = root[2]

            end
            
            
        else
            # find points that are outside vitelline
            balls_outside_id = balls[(balls[:,"x_t+1"].^2)./(R1^2) + (balls[:,"y_t+1"].^2)./(R2^2) .> 1, "ID"]
            balls_outside_id = balls_outside_id .+ 1 #for Julia since row number starts from 1
            
            # intersection of ellipse with line
            #for index, row in balls_outside.iterrows()
            for index in balls_outside_id
                
                row = balls[index,:]
                P1 = [row["x"], row["y"]]
                P2 = [row["x_t+1"], row["y_t+1"]]
                if P1 == P2
                    continue
                end
                
                slope = (P2[2] - P1[2])/(P2[1] - P1[1]) #slope
                intercept = P1[2] - slope*P1[1]

                function f!(F,x)
                    #return [line equation, ellipse equation]
                    #return( (x[1] - slope*x[0] - intercept , (x[0]/R1)**2 + (x[1]/R2)**2 - 1) )
                    F[1] = x[2] - slope*x[1] - intercept 
                    F[2] = (x[1]/R1)^2 + (x[2]/R2)^2 - 1
                end
                
                root = nlsolve(f!, P1).zero

                balls[index, "x_t+1"] = root[1]
                balls[index, "y_t+1"] = root[2]

            end

        end
                
                
        #storing final position of each particle
        balls[:,"x"] = balls[:,"x_t+1"]
        balls[:,"y"] = balls[:,"y_t+1"]
        
        if mod(N_iter,output_intervals) == 0 
            #change output intervals to output more or less number of timepoints in the dataframe
            println(string("timestep : ",string(t)))
            
            #balls_timeseries = pd.concat([balls_timeseries, balls], ignore_index=True)            
            #balls_timeseries = vcat(balls_timeseries, balls[:,cols_to_write]) 
            #save
            #filehandler = open(output_dir+'database.pickle', 'wb') 
            #pickle.dump(balls_timeseries, filehandler)
            CSV.write(string(output_dir,"N_iter_",string(N_iter),".csv"), balls[:,cols_to_write])

            #plot
            #plot_curve(balls, springs, scatter = True, vitelline_shape = vitelline_shape,
            #            R1 = R1, R2 = R2, save_fig=True, mitotic_domains_bool=False, 
            #            filename=output_dir+'niter_'+str(N_iter)+'.png',)

        end
            

        t = t + dt

    end
        
    println(string("exiting with timestep : ",string(t)))

    if mod(N_iter,output_intervals) != 0
        #last timepoint has not been included
        #balls_timeseries = pd.concat([balls_timeseries, balls], ignore_index=True)
        balls_timeseries = vcat(balls_timeseries, balls)
        #plot_curve(balls,springs,scatter = False,save_fig=True, mitotic_domains_bool=False, filename=output_dir+'niter_'+str(N_iter)+'.png')
        #filehandler = open(output_dir+'database.pickle', 'wb') 
        #pickle.dump(balls_timeseries, filehandler)
        CSV.write(string(output_dir,"N_iter_",string(N_iter),".csv"), balls[:,cols_to_write])
        CSV.write(string(output_dir,"final.csv"), balls[:,cols_to_write])
    end
    
    return balls_timeseries

end

#set working dir
map_index_dest=ARGS[1]
task_id= parse(Int,ARGS[2]) + 1 #because indexing starts from 1 in julia

#
#get the file with name of the directory
map_index_df = CSV.read(map_index_dest, DataFrame)
#rename directory
dirname = map_index_df[task_id, "folder_name"]
#set wd
#cd(current_wd)
println("directories set")

#import spring network
balls = CSV.read(string(dirname,"init_balls.csv"), DataFrame)
springs = CSV.read(string(dirname,"init_springs.csv"), DataFrame)
println("reading finished")

#import parameters
sim_params = CSV.read(string(dirname,"sim_params.csv"), DataFrame)

seed = parse(Int, sim_params[sim_params.variable .== "seed", :value][1])
t_final = parse(Float64, sim_params[sim_params.variable .== "t_final", :value][1])
t_CF = parse(Float64, sim_params[sim_params.variable .== "t_CF", :value][1])
add_CF_bool = Bool(parse(Float64, sim_params[sim_params.variable .== "add_CF_bool", :value][1]))
int_curvature_CF = parse(Float64, sim_params[sim_params.variable .== "int_curvature_CF", :value][1])
t_MD = parse(Float64, sim_params[sim_params.variable .== "t_MD", :value][1])
add_MD_bool = Bool(parse(Float64, sim_params[sim_params.variable .== "add_MD_bool", :value][1]))
int_curvature_MD = parse(Float64, sim_params[sim_params.variable .== "int_curvature_MD", :value][1])
dt = parse(Float64, sim_params[sim_params.variable .== "dt", :value][1])
vitelline_shape = sim_params[sim_params.variable .== "vitelline_shape", :value][1]
R = parse(Float64, sim_params[sim_params.variable .== "R", :value][1])
R_vert_frac = parse(Float64, sim_params[sim_params.variable .== "R_vert_frac", :value][1])
L = parse(Float64, sim_params[sim_params.variable .== "L", :value][1])
K = parse(Float64, sim_params[sim_params.variable .== "K", :value][1])
K_md = parse(Float64, sim_params[sim_params.variable .== "K_md", :value][1])
R_vit = parse(Float64, sim_params[sim_params.variable .== "R_vit", :value][1])
# True will throw error due to case of 1st letter , bend_energy = parse(Bool, sim_params[sim_params.variable .== "bend_energy", :value])
# True will throw error due to case of 1st letter , stretch_energy = parse(Bool, sim_params[sim_params.variable .== "stretch_energy", :value])
output_intervals = parse(Float64, sim_params[sim_params.variable .== "output_intervals", :value][1])
# True will throw error due to case of 1st letter , noise = parse(Float64, sim_params[sim_params.variable .== "noise", :value])
noise_scale = parse(Float64, sim_params[sim_params.variable .== "noise_scale", :value][1])
anneal_checkpoint = parse(Float64, sim_params[sim_params.variable .== "anneal_checkpoint", :value][1])
#mitotic_domains_bool = parse(Bool, sim_params[sim_params.variable .== "mitotic_domains_bool", :value])
#dirname = sim_params[sim_params.variable .== "output_dir", :value][1]
#dt_after_noise = sim_params[sim_params.variable .== "dt_after_noise", :value]
bend_energy = true
stretch_energy = true
noise = true
mitotic_domains_bool = true



println("Read parameters")
#seed = 0
Random.seed!(seed)



println("sending to minimize")
#convert this call to julia call
balls_timeseries = energy_minimize(balls,springs = springs, seed = seed, 
                                               t_final=t_final, dt = dt, 
                                               t_CF = t_CF, int_curvature_CF = int_curvature_CF, add_CF_bool = add_CF_bool,
                                               t_MD = t_MD, int_curvature_MD = int_curvature_MD, add_MD_bool = add_MD_bool,
                                               vitelline_shape = vitelline_shape, 
                                               R = R, R_vert_frac = R_vert_frac, #coords_vit = coords_vit,
                                               L=L, K= K, K_md = K_md, 
                                               R_vit = R_vit, #elastic_vitelline=elastic_vitelline,vitelline = vitelline, 
                                               bend_energy=bend_energy, stretch_energy = stretch_energy,
                                               #germ_band_pushing = germ_band_pushing, F_germ_band = F_germ_band, push_thresh = push_thresh, 
                                               output_intervals = output_intervals, 
                                               noise = noise, noise_scale=noise_scale, anneal_checkpoint = anneal_checkpoint,
                                               mitotic_domains_bool=mitotic_domains_bool,
                                               output_dir=dirname,#dt_after_noise=dt_after_noise,
                                               #dt_upp_thr = dt_upp_thr, dt_low_thr = dt_low_thr, #variable_dt_start = variable_dt_start,
                                               #no_noise_mitotic_domain = no_noise_mitotic_domain,
                                             )

println("returned")
