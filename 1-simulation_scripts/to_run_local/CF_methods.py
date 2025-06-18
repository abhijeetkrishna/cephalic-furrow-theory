#In this code we simulate an elastic ring inside a shell
import numpy as np
import pandas as pd
import math


import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib import collections  as mc
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import warnings
import matplotlib.gridspec as gridspec
warnings.filterwarnings('ignore')
import pickle
import os
import sys
import glob
import re

from scipy.special import ellipe
import scipy as sp
import scipy.optimize
from scipy.spatial import distance
from scipy.optimize import fsolve

def get_half_ellipse_xs_ys(R_vert_frac, R, Npoints):
    phi, ellipse_len = angles_in_ellipse(2*Npoints, R_vert_frac, 1)

    phi -= np.pi/2

    xs = 2*R/ellipse_len*np.cos(phi)
    ys = 2*R_vert_frac/ellipse_len*np.sin(phi)

    xs = xs[ys > 0]
    ys = ys[ys > 0]
    
    return(xs, ys)


def squeeze_md_particles(phi, R, R_vert_frac, mitotic_domains, Nballs, MD_growth, ellipse_len):
    
    print(mitotic_domains)
    
    md_inds = []

    if R_vert_frac == 0:
        
        xs = phi
        ys = np.zeros(len(phi))
        
        indices = np.arange(0, len(phi))
        phi_frac = indices*(1/(len(phi)-1))

        xs_del = None
        new_xs = None

        for md in mitotic_domains: 
            #plt.scatter(xs[(xs > md[0]) & (xs < md[1])], ys[(xs > md[0]) & (xs < md[1])], c = 'orange')

            xs_del_temp = xs[(phi_frac > md[0] - 0.5*1/Nballs) & (phi_frac < md[1] + 0.5*1/Nballs)]
            new_xs_temp = np.linspace(min(xs_del_temp), max(xs_del_temp), num = round(MD_growth*(len(xs_del_temp) - 1) + 1))

            if xs_del is not None:
                xs_del = np.concatenate((xs_del, xs_del_temp))
                new_xs = np.concatenate((new_xs, new_xs_temp))
            else: 
                new_xs = new_xs_temp.copy()
                xs_del = xs_del_temp.copy()

        xs = np.setdiff1d(xs, xs_del)
        xs = np.concatenate((xs, new_xs))

        xs.sort()
        
        MD_bool = np.isin(xs, new_xs)
        
        phi = xs

        ys = np.zeros(len(xs))
        
    else: 
        '''
        if R_vert_frac == 1: 
            phi = np.arange(0, np.pi, np.pi/Nballs)
            
            ellipse_len = 2*np.pi*R_vert_frac
            
        else: 
            phi, ellipse_len = angles_in_ellipse(2*Nballs, R_vert_frac, 1)
            phi -= np.pi/2
            phi = phi[(phi <= np.pi) & (phi >= 0)]
        '''
        
            

        phi[::-1].sort()


        indices = np.arange(0, len(phi))
        phi_frac = indices*(1/(len(phi)-1))

        phi_del = None
        new_phi = None

        for md in mitotic_domains: 

            temp_bool = (phi_frac > md[0] - 0.5*1/Nballs ) & (phi_frac < md[1] + 0.5*1/Nballs)

            phi_del_temp = phi[temp_bool]
            
            #round(MD_growth*(len(xs_del_temp) - 1) + 1)
            new_phi_temp = np.linspace(min(phi_del_temp), max(phi_del_temp), num = round(MD_growth*(len(phi_del_temp) - 1) + 1))

            if phi_del is not None:
                phi_del = np.concatenate((phi_del, phi_del_temp))
                new_phi = np.concatenate((new_phi, new_phi_temp))
            else: 
                new_phi = new_phi_temp.copy()
                phi_del = phi_del_temp.copy()



        phi = np.setdiff1d(phi, phi_del)
        phi = np.concatenate((phi, new_phi))

        phi[::-1].sort()
        
        MD_bool = np.isin(phi, new_phi)
        
        xs = 2*R/ellipse_len*np.cos(phi)
        ys = 2*R_vert_frac/ellipse_len*np.sin(phi)
        
        
        
    return(xs, ys, phi, MD_bool)


"""
#Function to initialize ring
def get_tissue(Nballs, K, L, R=1, R_vert_frac = 0, vitelline_shape = 'line', intrinsic_curvature = False,
               mitotic_domains_bool = False, mitotic_domains=[], 
               MD_growth_const_l0 = False, mitotic_growth = 1, 
               germ_band_push = 1,
               K_md = None, L_md = None, close_loop = True):
    
    if vitelline_shape == 'line':
        close_loop = False
    
    if K_md is None:
        K_md = K
    if L_md is None:
        L_md = L
    
    mass = 1

    balls_colnames=['ID', 'x', 'y', 'z', 'ax', 'ay', 'az', 'neighbours', 'spring1', 'spring2','row','column','stack', 'vx', 'vy', 'vz', 'mass']
    balls=pd.DataFrame(0, index=range(0), columns=range(len(balls_colnames)))
    balls.columns=balls_colnames
    balls.index = range(balls.shape[0])
    
    phi = None
    ellipse_len = None
    
    md_bool = None

    if vitelline_shape == 'circle':
        thetas = np.arange(0, 2*np.pi, 2*np.pi/Nballs)
        xs = [R*np.cos(theta)/(2*np.pi) for theta in thetas]
        ys = [R*np.sin(theta)/(2*np.pi) for theta in thetas]

        phi = thetas
        ellipse_len = 2*np.pi*R

    elif vitelline_shape == 'line':
        
        Nballs_squeezed = round(Nballs/(germ_band_push))
        print(Nballs_squeezed)

        if R_vert_frac == 0:
            xs = np.linspace(0,R,Nballs_squeezed) 
            ys = np.zeros(Nballs_squeezed)
            
            xs = xs[:Nballs]
            ys = ys[:Nballs]

            phi = xs
            ellipse_len = R

        elif R_vert_frac == 1:
            thetas = np.arange(0, np.pi, np.pi/Nballs_squeezed)
            thetas[::-1].sort()
            
            thetas = thetas[:Nballs]
            
            xs = [R*np.cos(theta)/(np.pi) for theta in thetas]
            ys = [R*np.sin(theta)/(np.pi) for theta in thetas]

            phi = thetas
            ellipse_len = 2*np.pi*R

        else: 
            phi, ellipse_len = angles_in_ellipse(2*(Nballs_squeezed), R_vert_frac, R)

            #print(ellipse_len)

            phi -= np.pi/2

            phi = phi[(phi <= np.pi) & (phi >= 0)]

            phi[::-1].sort()
            
            print(len(phi))

            phi = phi[:Nballs + 1]
            
            print(len(phi))

            xs = 2*R/ellipse_len*np.cos(phi)
            ys = 2*R_vert_frac/ellipse_len*np.sin(phi)


            
    if MD_growth_const_l0 and mitotic_domains_bool:
        xs, ys, phi, md_bool = squeeze_md_particles(phi, R, R_vert_frac, mitotic_domains, Nballs, mitotic_growth, ellipse_len)

    
    balls["x"] = xs 
    balls["y"] = ys
    balls["z"] = 0
    balls.row = 0
    balls.column = 0
    balls.stack = 0

    k = 1
    ns = 1


    balls[['vx', 'vy', 'vz']] = 0
    balls.mass = mass

    balls.ID = balls.index

    #empty lists for neighbours and springs, will be extended in loop below
    balls.neighbours = [[] for _ in range(len(balls))]
    balls.spring1 = [[] for _ in range(len(balls))]
    balls.spring2 = [[] for _ in range(len(balls))]

    #dataframe of springs
    #ID, X1, Y1, Z1, X2, Y2, Z2, k, Natural length, Extension in length, Ball1, Ball2
    springs_colnames=['ID', 'x1', 'y1', 'z1', 'x2', 'y2', 'z2', 
                      'k', 'l0', 'l1','dl', 
                      'ball1', 'ball2','type','viscoelastic_coeff']
    springs=pd.DataFrame(0, index=range(0), columns=range(len(springs_colnames)))
    springs.columns=springs_colnames
    
    
    #return(balls, springs)


    for point in balls.ID:

        #get an array of all neighbours
        #find all balls for which the difference in row is less than 2 AND difference in column is less than 2 AND difference in stack is less than 2
        #[row,col,stack] = [balls.loc[point,'row'], balls.loc[point,'column'], balls.loc[point,'stack']]
        ID = balls.loc[point,'ID']

        neighbours = [ID+1, ID-1]

        if ID == 0:
            if close_loop:
                neighbours = [len(balls),1]
            else:
                neighbours = [1]

        if ID == len(balls)-1:
            if close_loop:
                neighbours = [0,ID-1]
            else:
                neighbours = [ID-1]

        #print(point)
        #print(neighbours)

        #when stack starts working then try this:
        #neighbours = balls.loc[((np.abs(balls.row - row) < 2) & (np.abs(balls.column - col) < 2)) & (np.abs(balls.stack - stack) < 2),'ID'].values
        for neighbour in neighbours:
            #check whether the connection has already been saved, 
            #if not, add it as a spring
            if neighbour == point:
                continue
            if( balls.loc[point, 'neighbours'].count(neighbour) == 0):
                balls.loc[point, 'neighbours'].append(neighbour)
                balls.loc[neighbour, 'neighbours'].append(point)


                spring_id = springs.shape[0]

                balls.loc[point, 'spring1'].append(spring_id)
                balls.loc[neighbour, 'spring2'].append(spring_id)


                length = np.sqrt(sum(np.square(balls.loc[point, ['x', 'y', 'z']].values - 
                                               balls.loc[neighbour, ['x', 'y', 'z']].values)))

                row=pd.DataFrame([[spring_id] + 
                                  list(balls.loc[point, ['x', 'y', 'z']].values) + 
                                  list(balls.loc[neighbour, ['x', 'y', 'z']].values) + 
                                 [k, length, length, 0, point, neighbour, 0, ns]])

                row.columns=springs.columns
                springs=pd.concat([springs,row])
                springs.index = range(springs.shape[0])
                
    if mitotic_domains_bool:
        if md_bool is None:
            #balls['mitotic_domain'] = np.where(balls.ID.isin(mitotic_domains), True, False)
            #springs['mitotic_domain'] = np.where((springs.ball1.isin(mitotic_domains)) & (springs.ball2.isin(mitotic_domains)), True, False)

            #make a vector of uniform values - 
            #allows to find the coordinates of 
            #first and last particles in mitotic domains
            if R_vert_frac == 0:
                comp = np.linspace(0,R,Nballs)
                ref = np.linspace(0,R,Nballs)

            else: 
                if R_vert_frac == 1: 
                    phi_unif = np.arange(0, np.pi, np.pi/Nballs)

                    phi_unif[::-1].sort()

                    indices = np.arange(0, len(phi_unif))
                    comp = indices*(1/(len(phi_unif)-1))

                    ref = phi_unif


                else: 
                    phi_unif, ellipse_len = angles_in_ellipse(2*Nballs, R_vert_frac, 1)
                    phi_unif -= np.pi/2
                    phi_unif = phi_unif[(phi_unif <= np.pi) & (phi_unif >= 0)]

                    phi_unif[::-1].sort()

                    indices = np.arange(0, len(phi_unif))
                    comp = indices*(1/(len(phi_unif)-1))
                    ref = phi_unif

            md_border_coords = []
            md_bool = np.zeros(len(xs))

            for md in mitotic_domains:
                md_coords = ref[(comp > md[0] - 0.5*1/Nballs ) & (comp < md[1] + 0.5*1/Nballs )]

                md_bool = np.where((phi >= min(md_coords)) & (phi <= max(md_coords)), 1, md_bool)

        balls['mitotic_domain'] = md_bool
        springs['mitotic_domain'] = np.where((springs.ball1.isin(mitotic_domains)) & (springs.ball2.isin(mitotic_domains)), True, False)
    
    else: 
        balls['mitotic_domain'] = False
        springs['mitotic_domain'] = False

    balls['mitotic_domain'] = balls['mitotic_domain'].astype(float)

    #return(balls, springs)
       
    #set stretch and bending constants appropriately
    #L lives on springs, K lives on particles
    
    balls['L+1'] = np.where((balls.mitotic_domain).astype(bool) & (balls.iloc[(balls.ID + len(balls) + 1)%int(len(balls))].mitotic_domain.values).astype(bool), L_md, L)
    balls['L-1'] = np.where((balls.mitotic_domain).astype(bool) & (balls.iloc[(balls.ID + len(balls) - 1)%int(len(balls))].mitotic_domain.values).astype(bool), L_md, L)
    balls['K'] = np.where(balls.mitotic_domain, K_md, K)
    balls['K+1'] = np.where(balls.iloc[(balls.ID + len(balls) + 1)%int(len(balls))].mitotic_domain.values, K_md, K)
    balls['K-1'] = np.where(balls.iloc[(balls.ID + len(balls) - 1)%int(len(balls))].mitotic_domain.values, K_md, K)
    
    if intrinsic_curvature:
        balls['ko'] = 1/R
        balls["ko_i-1"] = np.array(balls["ko"][(np.array(balls.ID)+len(balls)-1)%int(len(balls))])
        balls["ko_i+1"] = np.array(balls["ko"][(np.array(balls.ID)+1)%int(len(balls))])
        
    else:
        balls['ko'] = 0
        balls["ko_i-1"] = 0
        balls["ko_i+1"] = 0
        
        
    balls['dr0'] = 1/(Nballs - 1)
    balls['dr0_i-1'] = 1/(Nballs - 1)
    
    
    #if vitelline_shape == 'circle':
    #    balls['dr0'] = 2*np.pi*R/Nballs #assuming the same natural length of each spring
    #    balls['dr0_i-1'] = 2*np.pi*R/Nballs
    #elif vitelline_shape == 'line':
    #    balls['dr0'] = R/Nballs
    #    balls['dr0_i-1'] = R/Nballs
    
    
    #update spring types as either face (default) or inplane
    #if z1 == z2
    springs['type'] = 'face'
    springs.loc[(springs.z1 == springs.z2), 'type'] = 'inplane'
    
    return([balls, springs])
"""

def get_tissue(Nballs = 100,
               K = 1,
               L = 1,
               vitelline_shape = 'line',
               R = 1,
               R_vert_frac = 0.4,
               vit_length = 1,
               germ_band_push = 0.2,
               mitotic_domains_bool = False,
               float_precision = 20,
               K_md = None, L_md = None,
               intrinsic_curvature = False,
               int_curvature_MD = None,
               int_curvature_CF = None,
               number_of_mitotic_domains = 3
              ):
    
    #from sympy import Ellipse, Point, Rational, Segment, Line, Circle, N
    #import scipy as sp
    
    # Example for using:
    
    
    # -------
    
    #get the vitelline shape

    #compute the initial and preferred lengths of the springs

    #start from the anterior pole of line and populate points 
    
    #stop when total desired initial length of tissue is achieved
    
    #if mitotic domains is true: then make function to add a neighbour between two points

    # --------
    
    #get the vitelline shape
    #vit  = Ellipse(Point(0,0), R, R_vert_frac)
    #try:
    #    vit_length = vit.circumference/2
    #except:
    #    vit_length = vit.length
        
    if R_vert_frac == 0:
        vit_length = 2*R
    else:
        m = 1 - (R/R_vert_frac)**2
        vit_length = R_vert_frac*sp.special.ellipeinc(np.pi, m)

    #compute the initial and preferred lengths of the springs
    #We will compute the preferred length as vit_length/Nballs
    #Then we will try to place the points on the ellipse
    #we will only add as many particles as will span the vitelline membrane
    #the number of particles may end up not being exactly Nballs
    
    #in unsqueezed form
    preferred_length = vit_length/Nballs
    #to squeeze the particles
    initial_length = preferred_length*(1 - germ_band_push)
    
    flag = 1
    
    #initialize a dataframe to add particles
    balls_colnames=['x', 'y', 'arclength'] #we can add neighbours later
    balls=pd.DataFrame(0, index=range(0), columns=range(len(balls_colnames)))
    balls.columns=balls_colnames
    
    #springs = None
    
    #start from the anterior pole of vit_membrane and populate points 
    
    #[x1, y1] = [x0, y0]
    
    Point0 = None
    Point1 = [-R, 0]
    
    while flag:
        
        #print('step', len(balls))
        
        
        if len(balls) == 0:
            
            #first particle has not been added yet
            #row=pd.DataFrame([list(Point1.coordinates)])
            #row = pd.DataFrame([[float(N(Point1[0] , float_precision)), float(N(Point1[1] , float_precision)), 0]])
            row = pd.DataFrame([[Point1[0], Point1[1], 0]])
            row.columns=balls.columns
            balls=pd.concat([balls,row])
            balls.index = range(balls.shape[0])
            
        #make a circle centered at Point1
        
        #find the points where this circle intersects with the vit object
        
        #select the point which is not already accounted for, this point is already saved as Point0
        
        #save this new point as Point2
        
        #next, save Point1 as Point0 and Point2 as Point1
        #print('checkpoint1')
        #make a circle centered at x1, y1
        
        if R_vert_frac == 0:
            #Point2 = Point(float(N(Point1[0], float_precision)) + preferred_length, 0)
            Point2 = [Point1[0] + initial_length, 0]
        else:
            #a = float(N(Point1[0], float_precision))
            a = Point1[0]
            #b = float(N(Point1[1], float_precision))
            b = Point1[1]
            c = initial_length
            
            #phi_guess = phi of current point + minor step in the anticlockwise direction
            phi_guess = np.arccos(-a/R) + np.pi/10
            
            [x_guess, y_guess] = [-R*np.cos(phi_guess), R*np.sin(phi_guess)]
            
            def func(x):
                #return [ellipse equation, circle equation]
                return( (x[0]/R)**2 + (x[1]/R_vert_frac)**2 - 1, (x[0] - a)**2 + (x[1] - b)**2 - c**2)
            
            root = fsolve(func, [x_guess, y_guess])
            #Point2 = Point(root[0], root[1])
            Point2 = root
            
        #check if Point2 is beyond desired tissue length
        #but to save time, do this check only when we are close to the number of desired particles
        if len(balls)/Nballs > 0.8:
            #get arclength of vit membrane from starting till the Point2
            #if the arclength > (1 - germ_band_push)*vit_length, break the loop
            if R_vert_frac == 0:
                #arclength = Segment(Point(-R, 0), Point2).length
                #arclength = float(N(Point2[0], float_precision)) + R
                arclength = Point2[0] + R
            elif R_vert_frac == 1:
                #circle = (-Rcos(t), Rsin(t))
                # => t = arccos(-x/R)
                #segment_angle = np.arccos(N(-Point2[0], float_precision)/R)
                #segment_angle = np.arccos(float(N(-Point2[0],float_precision))/R)
                segment_angle = np.arccos(-Point2[0]/R)
                arclength = segment_angle*R
                
                #print('segment_angle', int(100*segment_angle)/100)
                #print('arclength in terms of pi', int(100*segment_angle)/(np.pi*100))
            else:
                #ellipse = (R.cos(p), R_ver_frac.sin(p)) => Here x and y are swapped on purpose
                #arclength = incomplete elliptic integral of second kind going from p = 0 till the p

                #Point_reflected = [-Point2[0], Point2[1]] #swap x to -x
                #phi = np.arccos(Point_reflected[0]/R)
                #m = 1 - (R/R_vert_frac)**2
                #arclength = R_vert_frac*sp.special.ellipeinc(phi, m)
                
                #Point_reflected = [-float(N(Point2[0],float_precision)), float(N(Point2[1],float_precision))]
                Point_reflected = [-Point2[0], Point2[1]]
                phi = np.arccos(Point_reflected[0]/R)
                m = 1 - (R/R_vert_frac)**2
                arclength = R_vert_frac*sp.special.ellipeinc(phi, m)
                
                #print('arclength', arclength)
                #print('phi', phi)

            #stop when total desired initial length of tissue is achieved
            if arclength > (1 - germ_band_push)*vit_length:
                print("length achieved")
                break
                
            if Point2[1] < 0:
                print("reached x axis")
                break
                
            if len(balls)>1.1*Nballs:
                print('too many balls? Something does not seem right!')
                break
            
        #if len(balls) == Nballs - 1:
        #    flag = 0
        #    print("Number of balls achieved")
        
        
        #add new particle
        
        #[x_new, y_new] =  [float(N(Point2[0], float_precision)), float(N(Point2[1], float_precision))]
        [x_new, y_new] =  [Point2[0], Point2[1]]
        arclength_new = balls.loc[len(balls)-1, 'arclength'] + np.sqrt( (balls.loc[len(balls)-1,'x'] - x_new)**2 + (balls.loc[len(balls)-1,'y'] - y_new)**2 )
        row=pd.DataFrame([[x_new, y_new, arclength_new]]) #N(value, precision) converts the output of sympy to float
        row.columns=balls.columns
        balls=pd.concat([balls,row])
        balls.index = range(balls.shape[0])
        
        #next, save Point1 as Point0 and Point2 as Point1
        #shift point0 and point1
        #Point0 = Point(float(N(Point1[0], float_precision)), float(N(Point1[1], float_precision)))
        Point0 = [Point1[0], Point1[1]]
        #Point1 = Point(float(N(Point2[0], float_precision)), float(N(Point2[1], float_precision)))
        Point1 = [Point2[0],Point2[1]]
    
    balls['arclength_frac'] = balls['arclength']/balls.loc[len(balls)-1, 'arclength']
    balls['mitotic_domain_id'] = -1
    balls['mitotic_domain'] = False
    
    #if mitotic domains is true: then make function to add a neighbour between two points
    if number_of_mitotic_domains == 0:
        mitotic_domains_bool = False
        
    if mitotic_domains_bool:
        
        if number_of_mitotic_domains == 3:
            md_locations = [(0.1, 0.15), (0.2,0.25), (0.3,0.35)]
        elif number_of_mitotic_domains == 4:
            md_locations = [(0.1, 0.15), (0.2,0.25), (0.3,0.35), (0.4, 0.45)]
        
        #for md_location in md_locations:
            
            # find the first row which is bigger than start location
            # find the first row which is bigger than end location
            # for each row between these locations, add a row with the mid point
            
        #we need the location in terms of arclength where the mitotic domains will be
        
        row_pointer = -1
        in_md = False
        md_pointer = 0
        md_location = md_locations[md_pointer]
        
        while md_pointer < len(md_locations):
            
            row_pointer += 1
            if row_pointer >= len(balls):
                break
            
            arclength_frac = balls.loc[row_pointer, 'arclength_frac']
            
            md_location = md_locations[md_pointer]
            if arclength_frac >= md_location[0]:
                if arclength_frac <= md_location[1]:
                    in_md = True
                else:
                    #we have come outside the mitotic domain
                    in_md = False
                    md_pointer += 1
                    
            if in_md:
                #add a row
                balls.loc[row_pointer-1, 'mitotic_domain_id'] = md_pointer
                balls.loc[row_pointer-1, 'mitotic_domain'] = True
                
                [x_new, y_new] =  [(balls.loc[row_pointer-1,'x'] + balls.loc[row_pointer,'x'])/2, (balls.loc[row_pointer-1,'y'] + balls.loc[row_pointer,'y'])/2]
                arclength_new = balls.loc[row_pointer-1, 'arclength'] + np.sqrt( (balls.loc[row_pointer-1,'x'] - x_new)**2 + (balls.loc[row_pointer-1,'y'] - y_new)**2 )
                arclength_frac_new = arclength_new/balls.loc[len(balls)-1, 'arclength']
                mitotic_domain_id_new = md_pointer
                #here we add a new row at a fractional index which is converted to integer later
                balls.loc[row_pointer - 0.5] = x_new, y_new, arclength_new, arclength_frac_new, mitotic_domain_id_new, True
                balls = balls.sort_index().reset_index(drop=True)
                #increase row pointer
                row_pointer += 1
        
    balls['ID'] = balls.index
    # assigning preferred lengths
    balls['dr0'] = preferred_length
    balls['dr0_i-1'] = preferred_length
    balls['dr0_i+1'] = preferred_length
    
    # assigning elastic moduli
    
    balls['mitotic_domain'] = balls['mitotic_domain'].astype(float)
       
    #set stretch and bending constants appropriately
    #L lives on springs, K lives on particles
    
    if K_md is None:
        K_md = K
    if L_md is None:
        L_md = L
    balls['L+1'] = np.where((balls.mitotic_domain).astype(bool) & (balls.iloc[(balls.ID + len(balls) + 1)%int(len(balls))].mitotic_domain.values).astype(bool), L_md, L)
    balls['L-1'] = np.where((balls.mitotic_domain).astype(bool) & (balls.iloc[(balls.ID + len(balls) - 1)%int(len(balls))].mitotic_domain.values).astype(bool), L_md, L)
    balls['K'] = np.where(balls.mitotic_domain, K_md, K)
    balls['K+1'] = np.where(balls.iloc[(balls.ID + len(balls) + 1)%int(len(balls))].mitotic_domain.values, K_md, K)
    balls['K-1'] = np.where(balls.iloc[(balls.ID + len(balls) - 1)%int(len(balls))].mitotic_domain.values, K_md, K)
    
    if intrinsic_curvature:
        balls['ko'] = 1/R
        balls["ko_i-1"] = np.array(balls["ko"][(np.array(balls.ID)+len(balls)-1)%int(len(balls))])
        balls["ko_i+1"] = np.array(balls["ko"][(np.array(balls.ID)+1)%int(len(balls))])
        
    else:
        balls['ko'] = 0
        balls["ko_i-1"] = 0
        balls["ko_i+1"] = 0

    if int_curvature_MD is not None:

        balls.loc[balls['mitotic_domain'] == 1.0, 'ko'] = int_curvature_MD
        balls["ko_i-1"] = np.array(balls["ko"][(np.array(balls.ID)+len(balls)-1)%int(len(balls))])
        balls["ko_i+1"] = np.array(balls["ko"][(np.array(balls.ID)+1)%int(len(balls))])

    if int_curvature_CF is not None:
        balls.loc[ (balls['arclength_frac']>0.36) & (balls['arclength_frac']<0.38) , 'ko'] = int_curvature_CF
        balls["ko_i-1"] = np.array(balls["ko"][(np.array(balls.ID)+len(balls)-1)%int(len(balls))])
        balls["ko_i+1"] = np.array(balls["ko"][(np.array(balls.ID)+1)%int(len(balls))])

        
    springs = pd.DataFrame()

    for index, row in balls.iterrows():

        if index == len(balls)-1:
            break

        (ball1, ball2) = (index, index+1) #ball ID is same as index
        (x1, y1) = (row['x'], row['y'])
        (x2, y2) = (balls.loc[index+1, 'x'], balls.loc[index+1, 'y'])
        l1 = np.sqrt((x1-x2)**2 + (y1-y2)**2)

        if row['mitotic_domain'] == 1.0 and balls.loc[index+1, 'mitotic_domain'] == 1.0:
            md_status = True
        else:
            md_status = False

        springs_row = pd.DataFrame({'ball1':ball1,
                                    'ball2':ball2,
                                    'x1':[x1],
                                    'y1':[y1],
                                    'z1':[0],
                                    'x2':[x2],
                                    'y2':[y2],
                                    'z2':[0],
                                    'mitotic_domain':[md_status]
                                   })
        springs = pd.concat([springs, springs_row], ignore_index=True)
    
    springs['l0'] = preferred_length
    balls['z'] = 0.0
    balls['ID'] = range(len(balls))
    
    return([balls, springs])

def angles_in_ellipse(num, a, b):
    assert(num > 0)
    assert(a < b)
    angles = 2 * np.pi * np.arange(num) / num
    if a != b:
        e = (1.0 - a ** 2.0 / b ** 2.0) ** 0.5
        tot_size = sp.special.ellipeinc(2.0 * np.pi, e)
        arc_size = tot_size / num
        arcs = np.arange(num) * arc_size
        res = sp.optimize.root(
            lambda x: (sp.special.ellipeinc(x, e) - arcs), angles)
        angles = res.x 
    return angles, tot_size

#function to update position of particles in the springs database
def update_springs(springs,ball_positions,compute_lo=False):
    springs_ball1s=ball_positions.loc[pd.Series.tolist(springs.ball1)]
    springs_ball1s.columns=['x1','y1','z1']
    springs_ball1s.reset_index(drop=True, inplace=True)
    springs.loc[:,['x1','y1','z1']]=springs_ball1s

    springs_ball2s=ball_positions.loc[pd.Series.tolist(springs.ball2)]
    springs_ball2s.columns=['x2','y2','z2']
    springs_ball2s.reset_index(drop=True, inplace=True)
    springs.loc[:,['x2','y2','z2']]=springs_ball2s

    #change the l1 and dls for the springs
    springs_ball1s.columns=['x','y','z']
    springs_ball2s.columns=['x','y','z']
    disp=[springs_ball2s.x-springs_ball1s.x,springs_ball2s.y-springs_ball1s.y,springs_ball2s.z-springs_ball1s.z]
    length=disp[0]**2+disp[1]**2+disp[2]**2
    length=length.apply(lambda row: math.sqrt(row))
    springs.l1=length
    
    if compute_lo:
        springs.l0=springs.l1
    springs.dl=springs.l1-springs.l0
    
    
    return(springs)



def plot_curve(balls,springs, scatter = False, color_print = False, plot_vitelline = True, vitelline_shape = 'circle', R1 = 1, R2 = 2, R_vit = 1,
              plot_mitotic_sections = True, mitotic_domains_bool = False, save_fig = False, filename = 'plot.png'):

    import matplotlib.pyplot as plt
    #fig = plt.figure(figsize = (8,8))
    #ax = fig.gca()
    #ax.set_aspect('equal')
    #ax.scatter(balls.x, balls.y)
    
    
    springs = update_springs(springs, balls.loc[:,['x','y','z']])
    
    ball1_coords = np.column_stack((springs.x1.values,springs.y1.values)) #these are the coordinates of the initial condition actually
    ball2_coords = np.column_stack((springs.x2.values,springs.y2.values)) #these are the coordinates of the initial condition actually
    spring_coords = np.stack((ball1_coords, ball2_coords), axis = 1)


    #springs['delta_l'] = (springs.l0 - springs.l1)/springs.l1
    springs['col'] = np.where(springs.mitotic_domain, 1, -1)
    #cols = springs.delta_l
    cols = springs.col

    #make sure colorbar centered at 0
    vmin=min(cols)
    vmax=max(cols)
    thr = max(abs(vmin), abs(vmax))

    norm = mcolors.Normalize(vmin=-thr, vmax=thr)
    c = cm.ScalarMappable(norm=norm, 
                          #cmap = cm.viridis)
                          #cmap=cm.PiYG,
                          cmap = cm.bwr)

    colmap = [c.to_rgba(col) for col in cols]

    #if color_print:
    #    lc = mc.LineCollection(spring_coords, color = colmap, linewidths=2)
    #else:
    #    lc = mc.LineCollection(spring_coords, linewidths=2, color = 'grey', alpha = 0.5)
    fig = plt.figure()


    if vitelline_shape == 'circle':
        fig = plt.figure(figsize = (8,8))
    elif vitelline_shape == 'line':
        #print('line')
        fig = plt.figure(figsize = (20,2))


     
        
    ax = fig.gca()
    #ax.add_collection(lc)
    #ax.plot(balls.x, balls.y, linewidth=2, color = 'blue')
    
    if plot_vitelline:
        
        if vitelline_shape == 'circle':
            #cc = plt.Circle(( 0 , 0 ), R_vit )
            #ax.add_artist( cc )
            thetas = np.arange(0, 2*np.pi, 2*np.pi/100)
            xs = [R_vit*np.cos(theta) for theta in thetas]
            ys = [R_vit*np.sin(theta) for theta in thetas]
            xs.append(xs[0])#to close the loop for plotting
            ys.append(ys[0])#to close the loop for plotting

            ax.plot(xs, ys, color = 'grey', linewidth = 2, alpha = 0.5)
            
        elif vitelline_shape == 'line':

            if R2 == 0: 
                xs = np.linspace(-R1, R1, 50)
                ys = np.zeros(50)

                ax.plot(xs, ys, color = 'grey', linewidth = 2, alpha = 0.8)

            else: 
                thetas = np.arange(0, np.pi, np.pi/100)
                xs = [R1*np.cos(theta) for theta in thetas]
                ys = [R2*np.sin(theta) for theta in thetas]

                ax.plot(xs, ys, color = 'grey', linewidth = 2, alpha = 0.8)


    

    if scatter:
        ax.scatter(balls.x, balls.y, s = 0.5, color =  'black')
        
    ax.set_aspect('equal')
    ax.autoscale()
    ax.margins(0.1)

    if vitelline_shape == 'line':
        plt.ylim((-0.1, 1.1))

    if color_print:
        fig.colorbar(c)
        
    if mitotic_domains_bool:
        plt.scatter(balls[balls.mitotic_domain].x, balls[balls.mitotic_domain].y, c = 'r')
        plt.scatter(balls[~balls.mitotic_domain].x, balls[~balls.mitotic_domain].y, c = 'black')
        
    if plot_mitotic_sections:
        lc = mc.LineCollection(spring_coords, color = colmap, linewidths=4)
        ax.add_collection(lc)
        
    #plt.frame(False)
    ax.axis('off')
        
    if save_fig:
        plt.savefig(filename, bbox_inches = 'tight')
        

def return_closest_ellipse_points(balls_outside, coords_vit):

    balls_outside_x_sign = np.where(balls_outside.x >= 0, 1, -1)

    balls_outside['x_abs'] = np.abs(balls_outside.x)

    coords_query = np.stack((balls_outside.x_abs, balls_outside.y), axis = 1)

    closest_points = coords_vit[np.argmin(distance.cdist(coords_query,coords_vit), axis = 1)]
    
    dists = np.min(distance.cdist(coords_query,coords_vit), axis = 1)
    
    closest_points[:,0] *= balls_outside_x_sign
    
    return(closest_points, dists)


def energy_minimize(balls, springs = None, t_final = 100, dt = 0.01, dynamic_dt = False,
                    vitelline_shape = 'line', 
                    R = 1, R_vert_frac = 0, coords_vit = None,
                    K = 1, L = 0.1, k_vit = 1, tau = None, tol = 1e-12, dynamic_t_final = True,
                    debug = False, output_intervals = 10, 
                    R_vit = 1, vitelline = False, elastic_vitelline = False,
                    bend_energy = True, stretch_energy = True, 
                    noise = False, noise_scale = 1, anneal_checkpoint = None, no_noise_mitotic_domain = False,
                    growth = False, desired_growth = 1.1, growth_intervals = 0.5, growth_rate = 1.005, 
                    mitotic_domains_growth = False, mitotic_domains_bool = False, int_curve_checkpoint = None,
                    output_dir = '', 
                    dt_after_noise = None, dt_upp_thr = 10**(-2), dt_low_thr = 10**(-4), variable_dt_start = 0.001):
                    
    #this function just includes 

    first_anneal_bool = True

    if tau is None:
        tau = 1/L
    
    R1 = R
    R2 = R_vert_frac
    R_vit = R #this is abitrary - this is a constant value in the dynamics eqation and should not matter
    
    balls =balls.copy(deep=True)
    #springs =springs.copy(deep=True)
    balls_timeseries = pd.DataFrame()
    
    t = 0
    times = np.array([])
    N_iter = -1
    
    t_anneal_checkpoint = anneal_checkpoint*t_final
    #for annealing, we will stop the noise and increase dt little by little
    
    cols = []
    balls[cols] = 0
    
    while(t<t_final):

        N_iter = N_iter + 1
        energy = 0
        
        balls["t"] = t
        balls["N_iter"] = N_iter
        
        #balls[["x_i-2", "y_i-2", "x_i-1", "y_i-1", "x_i+1", "y_i+1", "x_i+2", "y_i+2"]] = 0
        col_to_set_zero = ["x_i-2", "y_i-2", "x_i-1", "y_i-1", "x_i+1", "y_i+1", "x_i+2", "y_i+2", "curvature", "curvature_i-1", "curvature_i+1", "d_r_mod", "d_r_mod_i-1","d_r_mod_i+1","d_r_mod_alternate","d_r_mod_alternate_i-1","d_r_mod_alternate_i+1", "d_curvature_x", "d_curvature_y", "d_curvature_i-1_x","d_curvature_i+1_x","d_curvature_i-1_y","d_curvature_i+1_y","d_energy_x", "d_energy_y", "noise_x", "noise_y"]
        for col in col_to_set_zero:
            balls[col] = 0
        
        balls["x_i-2"][2:] += np.array(balls["x"][(np.array(balls.ID[2:])-2)])
        balls["y_i-2"][2:] += np.array(balls["y"][(np.array(balls.ID[2:])-2)])
        balls["x_i-1"][1:] += np.array(balls["x"][(np.array(balls.ID[1:])-1)])
        balls["y_i-1"][1:] += np.array(balls["y"][(np.array(balls.ID[1:])-1)])
        balls["x_i+1"][:-1] += np.array(balls["x"][(np.array(balls.ID)[:-1]+1)])
        balls["y_i+1"][:-1] += np.array(balls["y"][(np.array(balls.ID)[:-1]+1)])
        balls["x_i+2"][:-2] += np.array(balls["x"][(np.array(balls.ID)[:-2]+2)])
        balls["y_i+2"][:-2] += np.array(balls["y"][(np.array(balls.ID)[:-2]+2)])
        
        balls["r"] = np.array(np.sqrt((balls["x"])**2+(balls["y"])**2))
        
        # Measure curvature
        
        #balls[["curvature", "curvature_i-1", "curvature_i+1", "d_r_mod"]] = 0
        
        balls["curvature"][1:-1] += np.array((4* ((balls["x_i+1"] - balls["x_i-1"])*(balls["y_i+1"] - 2*balls["y"] + balls["y_i-1"]) - (balls["y_i+1"] - balls["y_i-1"])*(balls["x_i+1"] - 2*balls["x"] + balls["x_i-1"]))/(((balls["x_i+1"] - balls["x_i-1"])**2 + (balls["y_i+1"] - balls["y_i-1"])**2)**(3/2))))[1:-1]
        balls["curvature_i-1"][1:] += np.array(balls["curvature"][(np.array(balls.ID[1:])-1)])
        balls["curvature_i+1"][:-1] += np.array(balls["curvature"][(np.array(balls.ID[:-1])+1)])
        balls["d_r_mod"][:-1] += np.array(((balls["x_i+1"]-balls["x"])**2 + (balls["y_i+1"]-balls["y"])**2)**0.5)[:-1]
        
        # non boundary nodes
        #balls[["d_r_mod_i-1","d_r_mod_i+1","d_r_mod_alternate","d_r_mod_alternate_i-1","d_r_mod_alternate_i+1"]] = 0
        idx = ~np.isin(np.arange(len(balls)), [0,1,len(balls)-2, len(balls)-1])
        
        balls["d_r_mod_i-1"][idx] += np.array(((balls["x_i-1"]-balls["x"])**2 + (balls["y_i-1"]-balls["y"])**2)**0.5)[idx]
        balls["d_r_mod_i+1"][idx] += np.array(((balls["x_i+2"]-balls["x_i+1"])**2 + (balls["y_i+2"]-balls["y_i+1"])**2)**0.5)[idx]
        balls["d_r_mod_alternate"][idx] += np.array((((balls["x_i+1"]-balls["x_i-1"])**2 + (balls["y_i+1"]-balls["y_i-1"])**2)**0.5)/2)[idx]
        balls["d_r_mod_alternate_i-1"][idx] += np.array((((balls["x"]-balls["x_i-2"])**2 + (balls["y"]-balls["y_i-2"])**2)**0.5)/2)[idx]
        balls["d_r_mod_alternate_i+1"][idx] += np.array((((balls["x_i+2"]-balls["x"])**2 + (balls["y_i+2"]-balls["y"])**2)**0.5)/2)[idx]
        
        #compute curvature gradient
        #balls[["d_curvature_x", "d_curvature_y"]] = 0
        balls["d_curvature_x"][idx] += np.array((balls["y_i+1"] - balls["y_i-1"])/(balls["d_r_mod_alternate"]**3))[idx]
        balls["d_curvature_y"][idx] += np.array(-(balls["x_i+1"] - balls["x_i-1"])/(balls["d_r_mod_alternate"]**3))[idx]
        
        # adding curvature gradient of neighbours
        #balls[["d_curvature_i-1_x","d_curvature_i+1_x","d_curvature_i-1_y","d_curvature_i+1_y"]] = 0
        balls["d_curvature_i-1_x"][idx] += np.array(-(balls['y_i-1'] - balls['y_i-2'])/((balls["d_r_mod_alternate_i-1"])**3) - 3*(balls['x'] - balls['x_i-2'])*balls["curvature_i-1"]/((2*balls["d_r_mod_alternate_i-1"])**2))[idx]
        balls["d_curvature_i+1_x"][idx] += np.array(-(balls['y_i+2'] - balls['y_i+1'])/((balls["d_r_mod_alternate_i+1"])**3) + 3*(balls['x_i+2'] - balls['x'])*balls["curvature_i+1"]/((2*balls["d_r_mod_alternate_i+1"])**2))[idx]
        balls["d_curvature_i-1_y"][idx] += np.array((balls['x_i-1'] - balls['x_i-2'])/((balls["d_r_mod_alternate_i-1"])**3) - 3*(balls['y'] - balls['y_i-2'])*balls["curvature_i-1"]/((2*balls["d_r_mod_alternate_i-1"])**2))[idx]
        balls["d_curvature_i+1_y"][idx] += np.array((balls['x_i+2'] - balls['x_i+1'])/((balls["d_r_mod_alternate_i+1"])**3) + 3*(balls['y_i+2'] - balls['y'])*balls["curvature_i+1"]/((2*balls["d_r_mod_alternate_i+1"])**2))[idx]
        
        #compute energy gradient 
        #balls[["d_energy_x", "d_energy_y"]] = 0
        
        if stretch_energy:
          
            balls["d_energy_x"][idx] += np.array(balls['L+1']*(balls["d_r_mod"]/balls['dr0'] - 1)*(balls["x"] - balls["x_i+1"])/(balls["d_r_mod"]) + 
                                             balls['L-1']*(balls["d_r_mod_i-1"]/balls['dr0_i-1'] - 1)*(balls["x"] - balls["x_i-1"])/(balls["d_r_mod_i-1"]))[idx]
        
            balls["d_energy_y"][idx] += np.array(balls['L+1']*(balls["d_r_mod"]/balls['dr0'] - 1)*(balls["y"] - balls["y_i+1"])/(balls["d_r_mod"]) + 
                                             balls['L-1']*(balls["d_r_mod_i-1"]/balls['dr0_i-1'] - 1)*(balls["y"] - balls["y_i-1"])/(balls["d_r_mod_i-1"]))[idx]
    
        #adding force due to bending energy
        if bend_energy:
            
            balls["d_energy_x"][idx] += np.array(2*balls['K-1']*(balls["curvature_i-1"]-balls["ko_i-1"])*balls["d_curvature_i-1_x"]*balls["dr0_i-1"] + 
                                            2*balls['K']*(balls["curvature"] - balls["ko"])*balls["d_curvature_x"]*balls["dr0"] + 
                                            2*balls['K+1']*(balls["curvature_i+1"] - balls["ko_i+1"])*balls["d_curvature_i+1_x"]*balls["dr0_i+1"]
                                            )[idx]

            balls["d_energy_y"][idx] += np.array(2*balls['K-1']*(balls["curvature_i-1"]-balls["ko_i-1"])*balls["d_curvature_i-1_y"]*balls["dr0_i-1"] + 
                                            2*balls['K']*(balls["curvature"] - balls["ko"])*balls["d_curvature_y"]*balls["dr0"] + 
                                            2*balls['K+1']*(balls["curvature_i+1"] - balls["ko_i+1"])*balls["d_curvature_i+1_y"]*balls["dr0_i+1"]
                                            )[idx]
        
        #compute displacement
        balls["dx_t"] = -(balls["d_energy_x"])*((R_vit)/(tau * L))
        balls["dy_t"] = -(balls["d_energy_y"])*((R_vit)/(tau * L)) #tau is set so that tau*L is 1
        
        if (dt_after_noise is not None) and (t > t_anneal_checkpoint):

            dt_step = 0.0003 #typically we keep dt = 0.0001 and increase by 0.0003 until we reach 0.001
            if first_anneal_bool:
                #Just after we switch off noise, we don't want to increase dt
                #after t_anneal_step we increase dt (basically next time)
                print('switching off noise but not increasing dt')
                dt_step = 0
                first_anneal_bool = False
                noise = False #now we switch off noise

            t_anneal_step = 1000

            print('changing dt')
            dt = dt + dt_step
            print('new dt :' + str(dt))
            t_anneal_checkpoint = t_anneal_checkpoint + t_anneal_step*dt 
            print('new anneal checkpoint', t_anneal_checkpoint)
            #this will ensure that this if block is accessed only when t_anneal_step number of iterations have passed

            if dt >= dt_after_noise:
                print('target dt reached so no more increasing dt')
                dt_after_noise = None #this will stop the process of increasing dt

        balls["x_t+1"] = balls["x"] + dt*balls["dx_t"] 
        balls["y_t+1"] = balls["y"] + dt*balls["dy_t"]

        
        #add noise
        if noise and t< t_anneal_checkpoint:
            print('noise!')
            #noise_arr = np.abs(np.random.normal(loc = 0, scale = noise_scale, size = len(balls)))
            #balls['x'] = balls['x'] - np.sign(balls['x']) * noise_arr
            #balls['y'] = balls['y'] - np.sign(balls['y']) * noise_arr
            balls['noise_x'][idx] += np.random.normal(loc = 0, scale = noise_scale, size = len(balls))[idx]
            balls['noise_y'][idx] += np.random.normal(loc = 0, scale = noise_scale, size = len(balls))[idx]
            #make the noise for mitotic domains zero
            if no_noise_mitotic_domain:
                #find the indices of mitotic domains and make the noise zero
                balls.loc[balls['mitotic_domain'] == 1.0,'noise_x'] = 0
                balls.loc[balls['mitotic_domain'] == 1.0,'noise_y'] = 0

            balls["x_t+1"] = balls["x_t+1"] + balls['noise_x']
            balls["y_t+1"] = balls["y_t+1"] + balls['noise_y']
            
            
        #adding vitelline repulsion
        
        if R_vert_frac == 0:
            # find points that are outside vitelline
            balls_outside = balls[balls['y_t+1']>0]
            
            # intersection of two straight lines
            for index, row in balls_outside.iterrows():
                
                P1 = [row['x'], row['y']]
                P2 = [row['x_t+1'], row['y_t+1']]
                
                slope = (P2[1] - P1[1])/(P2[0] - P1[0]) #slope
                intercept = P1[1] - slope*P1[0]

                #[x_guess, y_guess] = P1
                
                def func(x):
                    #return [line equation, line equation]
                    return( (x[1] - slope*x[0] - intercept , x[1]) )
                
                root = fsolve(func, P1)
                
                balls.loc[index, 'x_t+1'] = root[0]
                balls.loc[index, 'y_t+1'] = root[1]
            
            
        else:
            # find points that are outside vitelline
            balls_outside = balls[balls['x_t+1']**2/R1**2 + balls['y_t+1']**2/R2**2 > 1]
            
            # intersection of ellipse with line
            for index, row in balls_outside.iterrows():
                
                P1 = [row['x'], row['y']]
                P2 = [row['x_t+1'], row['y_t+1']]
                
                slope = (P2[1] - P1[1])/(P2[0] - P1[0]) #slope
                intercept = P1[1] - slope*P1[0]

                #[x_guess, y_guess] = P1
                
                def func(x):
                    #return [line equation, ellipse equation]
                    return( (x[1] - slope*x[0] - intercept , (x[0]/R1)**2 + (x[1]/R2)**2 - 1) )
                
                root = fsolve(func, P1)
                
                balls.loc[index, 'x_t+1'] = root[0]
                balls.loc[index, 'y_t+1'] = root[1]
                
                
        #storing final position of each particle
        balls['x'] = balls['x_t+1']
        balls['y'] = balls['y_t+1']
        

        if N_iter%output_intervals == 0: #change output intervals to output more or less number of timepoints in the dataframe
            print("timestep : "+str(t))
            
            balls_timeseries = pd.concat([balls_timeseries, balls], ignore_index=True)
            #save
            filehandler = open(output_dir+'database.pickle', 'wb') 
            pickle.dump(balls_timeseries, filehandler)
            #plot
            #plot_curve(balls, springs, scatter = True, vitelline_shape = vitelline_shape,
            #            R1 = R1, R2 = R2, save_fig=True, mitotic_domains_bool=False, 
            #            filename=output_dir+'niter_'+str(N_iter)+'.png',)
            

        t = t + dt
        
    print("exiting with timestep : "+str(t))

    balls_timeseries = pd.concat([balls_timeseries, balls], ignore_index=True)
    plot_curve(balls,springs,scatter = False,save_fig=True, mitotic_domains_bool=False, filename=output_dir+'niter_'+str(N_iter)+'.png')
    filehandler = open(output_dir+'database.pickle', 'wb') 
    pickle.dump(balls_timeseries, filehandler)
    
    return(balls_timeseries, np.sum(balls["d_r_mod_i-1"]))
             
def get_subplot(fig, ax, x, y, 
                xlabel_font = 20, xlabel = "N_iter", xlabelpad = 10,
                ylabel_font = 20, ylabel = "Prop", ylabelpad = 5,
                xlim = None, ylim = None,
                lw = 2, color = "blue",
               ):

    ax.plot(x,y,lw = lw, color = color)
    ax.set_xlabel(xlabel, fontsize = xlabel_font)
    ax.set_ylabel(ylabel, fontsize = ylabel_font, rotation = 'vertical', labelpad = ylabelpad)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    #ax.tick_params(axis='both', which='major', labelsize=15)
    #ax.set_title(title, fontsize = 20, pad = 25)
    #ax.grid()
    return(fig,ax)


def load_timepoints(glob_path, last_timepoint_only = False, npart = 120):

    for i, path in enumerate(glob.glob(glob_path)):

        try: 
            #df_temp = pickle.load(open(glob.glob(path + 'database.pickle')[0], "rb"))
            df_temp = pd.read_csv(path + '/final.csv')
            initial_df = pd.read_csv(path + '/init_balls.csv')
        except: 
            print('files corrupted or empty path: ' + path)
            #print(path)
            continue
            
        initial_df_common = initial_df[["ID", "mitotic_domain", "mitotic_domain_id","dr0", "K", "L+1", "L-1","ko","arclength", "arclength_frac"]]
        df_temp = pd.concat([initial_df_common, df_temp], axis = 1)
        df_temp["init_x"] = initial_df["x"]
        df_temp["init_y"] = initial_df["y"]

        #split path name (path name will be split by numerical values)
        #to get all parameter names
        param_names = re.split(r'_[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?_*', path.split('/')[-2])
        #usually get an empty string as a result of splitting
        if '' in param_names:
            param_names.remove('')
        #for each param name, grab its value
        for param in param_names:   
            param_val = re.findall(param + '_[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?_*', path)
            val = float(re.findall('[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', param_val[0])[-1])
            #for each param, make a column in the dataframe and save the parameter value there
            df_temp[param] = val

        df_temp["path"] = path
        
        if i == 0:
            df = df_temp
        else: 
            df = pd.concat(
                            [df, df_temp],
                            axis=0,
                            join="outer",
                            ignore_index=True,
                            keys=None,
                            levels=None,
                            names=None,
                            verify_integrity=False,
                            copy=True,
                        )
            
    return df


def get_closest_ellipse_pt(semi_major, semi_minor, p): 
    #from https://stackoverflow.com/questions/22959698/distance-from-given-point-to-given-ellipse
    
    px = abs(p[0])
    py = abs(p[1])

    tx = 0.707
    ty = 0.707

    a = semi_major
    b = semi_minor

    for x in range(0, 3):
        x = a * tx
        y = b * ty

        ex = (a*a - b*b) * tx**3 / a
        ey = (b*b - a*a) * ty**3 / b

        rx = x - ex
        ry = y - ey

        qx = px - ex
        qy = py - ey

        r = math.hypot(ry, rx)
        q = math.hypot(qy, qx)

        tx = min(1, max(0, (qx * r / q + ex) / a))
        ty = min(1, max(0, (qy * r / q + ey) / b))
        t = math.hypot(ty, tx)
        tx /= t 
        ty /= t 

    return (math.copysign(a * tx, p[0]), math.copysign(b * ty, p[1]))
    

















