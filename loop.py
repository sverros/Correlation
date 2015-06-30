import numexpr as ne
import numpy as np
from scipy import linalg
import math
import matplotlib.pyplot as plt
import datetime
from neicio.readstation import readStation
from neicio.shake import ShakeGrid
from openquake.hazardlib.correlation import JB2009CorrelationModel
from openquake.hazardlib.correlation import BaseCorrelationModel
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo import Point
from openquake.hazardlib.geo.geodetic import geodetic_distance
from openquake.hazardlib.imt import from_string
import time
from matplotlib import cm
from neicio.gmt import GMTGrid
import time
import sys


def main(var, r, voi, rand, vscorr, intensity_factor):
    np.set_printoptions(linewidth = 200)
    #####
    # Main program for computing spatial correlation
    # IN: var- variables dictionary from initialize function. Contains M,N,K,site_collection_SM, site_collection_station
    #          uncertaintydata, data, location_lon/lat_g
    #     r - radius
    #     voi- variable of interest, i.e. PGA
    #     rand- normally distributed random variable of size var['M']*var['N']
    #     vscorr- boolean, whether vs30 is correlated or uncorrelated
    #OUT: cor- grid of spatially correlated epsilon
    #     data- grid of ShakeMap data
    #     data_new- data with added spatial correlation
    #     grid_arr- array for storing grid indices for multiple realizations
    #     mu_arr- array for storing Sig21.T*Sig11inv for multiple realizations
    #     sigma_arr- array for storing sigma for multiple realizations
    #####

    start = time.time()    
    OL_time = 0
    IL_time = 0
    full_dist_time = 0
    out_t = 0
    red_dist = 0
    base_time = 0
    base_t = 0
    other_time = 0
    other_t = 0
    corr_time = 0
    M = var['M']
    N = var['N']
    K = var['K']


    JB_cor_model = JB2009CorrelationModel(vscorr)

    # Initialize vector where data_new will be stored
    X = np.zeros([M*N,1])

    # Initialize vectors for storing data for multiple realizations
    grid_arr = [None] * (M*N)
    mu_arr = [None] * (M*N)
    sigma_arr = np.zeros([M*N,1])
    rand_arr = np.zeros([M*N,1])
    
    # Get spcing of horozontal and vertical points
    ld  = set_up_grid_dist(M,N,var['site_collection_SM'])
    
    pre_loop_time = time.time() - start

    for i in range(0,M):
        OL_start = time.time()
        
        # Find the number of points in radius horozontally and vertically for each row
        vhva = calc_vert_hor(i, r, ld['l'], ld['d'])
        full_dist_t = time.time()

        # Calculate the full distance matrix for each row
 
        dist = calc_full_dist(vhva['vert'], vhva['hor'], N, var['site_collection_SM'])
        full_dist_time += time.time() - full_dist_t
        first_time_per_row = 1
        OL_time += time.time() - OL_start
        for j in range(0,N):
            IL_start = time.time()
            num = i*N+j
        
            t = time.time()
            # Find the reduced distance matrix 
            dist_calc = reduce_distance(j, vhva['vert'], vhva['hor'], N, dist['distance_matrix'], dist['grid_indices'], num)
            
            red_dist += time.time() - t
            
            t = time.time()
            # Include stations in distance matrix and find the indices of the points within the radius
            out = inc_stations(j, i, N, K, r, var['site_collection_SM'], var['site_collection_station'], 
                               dist_calc['dist_mat'], X, dist_calc['inc_indices'])
            
            out_t += time.time() - t
            
            if np.size(dist_calc['inc_indices']) == 1:
            
                # no conditioning points, correlation value is random                                                                                          
                X[num] = rand[num]
            
                # Store for multiple realizations                                                                                                                 
                grid_arr [num] = np.zeros(0)
                mu_arr   [num] = np.zeros(0)
                sigma_arr[num] = 1
                rand_arr [num] = X[num]
            else:
                # Check if reduced distance matrix is full distance matrix
#                if ((vhva['vert'] == 1 and dist_calc['num_indices'] == vhva['hor']+1)or(vhva['vert'] != 1 and \
#                       dist_calc['num_indices'] == 2*vhva['hor']+1)) and (np.size(out['inc_sta_indices']) == 0):
#                    # If this is the first full distance matrix per row, calculate base case
#                    if first_time_per_row == 1:
#                        t = time.time()
#                        base = calculate_corr(out['dist_mat'], voi, JB_cor_model, var, out['inc_sta_indices'], intensity_factor)
#                        first_time_per_row = 0
#                        base_time += time.time() - t
#                    t = time.time()
#                    mu  = base['Sig12'].T*base['Sig11inv']*(out['x'])
#                    
#                    rand_num = rand[num]
#                    X[num] = mu+rand_num*base['R']
#            
#                    # Store for multiple realizations                                                                                                       
#                    grid_arr [num] = dist_calc['inc_ind'][0:-1]
#                    mu_arr   [num] = base['Sig12'].T*base['Sig11inv']
#                    sigma_arr[num] = base['R']
#                    rand_arr [num] = rand_num
#                    base_t += time.time() - t
#                else:
                t = time.time()
#                other = calculate_corr(out['dist_mat'], voi, JB_cor_model, var, out['inc_sta_indices'], intensity_factor)
                other = calculate_corr(dist_calc['dist_mat'], voi, JB_cor_model, var, [], intensity_factor)
                other_time += time.time() - t
                    #mu = other['Sig12'].T*other['Sig11inv']*(out['x']- np.mean(X[0:i*N+j]))
                t = time.time()
#                mu = other['Sig12'].T*other['Sig11inv']*(out['x'])
                x = X[dist_calc['inc_indices'][0:-1]]
                mu = other['Sig12'].T*other['Sig11inv']*x

                rand_num = rand[num]
                X[num] = mu+rand_num*other['R']
                
                    # Store for multiple realizations                                                                             
                grid_arr [num] = dist_calc['inc_indices'][0:-1]
                mu_arr   [num] = other['Sig12'].T*other['Sig11inv']
                sigma_arr[num] = other['R']
                rand_arr [num] = rand_num
                other_t += time.time()-t
            IL_time += time.time() - IL_start
            if np.mod(i*N+j,5000) == 0:
                print 'Finishing step:', i*N+j
                sys.stdout.flush()

    DATA = var['data']
    COR = np.reshape(X, [M,N]) #units epsilon
    #Multiply by uncertainty                                                                                                   
    X = np.multiply(COR, var['uncertaintydata']) # ln(pctg)
                
    DATA_NEW = np.log(DATA) + X
    DATA_NEW = np.exp(DATA_NEW)

    end = time.time() - start
    print 'Total Time', end
    sys.stdout.flush()
    print 'Pre loop Time', pre_loop_time
    print 'Inner loop time', IL_time
    print 'Outer loop time', OL_time
    print 'Reduced Dist', red_dist
    print 'out fn', out_t
    print 'base_time', base_time, base_t
    print 'other time', other_time, other_t
    return {'cor':COR, 'data':DATA, 'data_new':DATA_NEW, 'grid_arr':grid_arr, 'mu_arr':mu_arr, 'sigma_arr':sigma_arr, 'rand_arr':rand_arr}



def calculate_corr(dist_mat, voi, JB_cor_model, var, inc_sta_indices, intensity_factor = 0.9):
    #####
    # Calculates correlation model for distance matrix and voi
    # IN: dist_mat- reduced distance matrix
    #     voi- variable of interest
    #     JB_cor_model- correlation model from correlation in oq-hazardlib
    #OUT: Sig12, Sig11inv- partitions of correlation matrix
    #     R - Sqrt of sigma
    #####
    correlation_model = JB_cor_model._get_correlation_model(dist_mat, from_string(voi))
    
    intensity = var['intensity'][inc_sta_indices]
    if np.size(intensity) != 0:
        for i in range(0,np.size(intensity)):
            if intensity[i] == 1:
                correlation_model[i,i+1:] = correlation_model[i,i+1:].copy()*intensity_factor
                correlation_model[i+1:,i] = correlation_model[i+1:,i].copy()*intensity_factor

    Sig11 = np.mat(correlation_model[0:-1, 0:-1])
    Sig12 = np.mat(correlation_model[0:-1, -1]).T
    Sig22 = np.mat(correlation_model[-1,-1])

    Sig11inv = np.mat(np.linalg.pinv(Sig11))
    sigma = Sig22 - (Sig12.T*Sig11inv*Sig12)
    R = np.sqrt(sigma)
    
    return {'Sig12':Sig12, 'Sig11inv':Sig11inv, 'R':R}

def inc_stations(j,i,N,K,r,site_collection_SM, site_collection_station, dist_mat, X, inc_indices):
    #####
    # If there are stations included within the radius for a point, this function will add those stations to the 
    # distance matrix and determine the array of points included in the radius, x
    # IN: i,j- current points row and column 
    #     N,K - number of points in row and total number of stations
    #     r- radius 
    #     site_collection_SM/station- site collections for ShakeMap and station data
    #     dist_mat- reduced distance matrix
    #     X- array of previously calculated correlation values
    #     inc_ind- indices of included points
    #     inc_indices- total number of points in the top most row of distance matrix
    #OUT: dist_mat- reduced distance matrix, modified to include stations
    #     x- array of points in X included in radius and stations included 
    #     inc_sta_indices- indices of stations included in the radius
    #####

    num = i*N+j
    
    # Compute the distances for all stations to the grid point we're looking at                                                                                                                
    dist_sta_sit = np.array(geodetic_distance(site_collection_SM.lons[j+i*N], site_collection_SM.lats[j+i*N],
                                              site_collection_station.lons[0:K], site_collection_station.lats[0:K]))
        
    # Find which of those stations are in the radius we are considering
    inc_sta_indices = np.where(dist_sta_sit < r)
    if np.size(inc_sta_indices) != 0:
            
        station_distance_matrix = np.zeros([np.size(inc_sta_indices), np.size(inc_sta_indices)+np.size(inc_indices)])
        # Calculate distance between each included station and all included grid points, then calculate the distance
        # from each included station to every other included station
        for eta in range(0, np.size(inc_sta_indices)):
            for beta in range(0,np.size(inc_indices)):
                station_distance_matrix[eta,np.size(inc_sta_indices) + beta] = geodetic_distance(
                    site_collection_station.lons[inc_sta_indices[0][eta]], site_collection_station.lats[inc_sta_indices[0][eta]], 
                    site_collection_SM.lons[inc_indices[beta]], site_collection_SM.lats[inc_indices[beta]])
            for beta in range(0, np.size(inc_sta_indices)):
                station_distance_matrix[eta, beta] = geodetic_distance(
                    site_collection_station.lons[inc_sta_indices[0][eta]], site_collection_station.lats[inc_sta_indices[0][eta]],
                    site_collection_station.lons[inc_sta_indices[0][beta]], site_collection_station.lats[inc_sta_indices[0][beta]])
            
        # Concatenate the station distance matrix with the modified distance matrix, dist_mat
        dist_mat = np.concatenate((station_distance_matrix[:, np.size(inc_sta_indices):], dist_mat), axis=0)
        dist_mat = np.concatenate((station_distance_matrix.T, dist_mat), axis=1)

        # x: vector of previously calculated covariance values
        x = np.concatenate((np.zeros([np.size(inc_sta_indices),1]),X[inc_indices]), axis = 0)
        x = np.mat(x[0:-1])
            
    else:
        # x: vector of previously calculated covariance values
        x = X[inc_indices]
        x = np.mat(x[0:-1])
    
    return {'dist_mat':dist_mat, 'x':x, 'inc_sta_indices':inc_sta_indices}

def reduce_distance(j, vert, hor, N, distance_matrix, grid_indices, num):
    # Find which columns/rows in the distance matrix to keep
    # IN: j- points column
    #     vert- number of rows included in the radius
    #     hor- number of columns included in radius
    #     added_vert- number of rows in between first row and first included row 
    #     N- number of points in row
    #     distance_matrix- full distance matrix
    #     grid_indices- indices included in full distance matrix
    #OUT: dist_mat- reduced distance matrix
    #     inc_indices- number of points in top most row of dist_mat
    #     inc_ind - indices of points included in dist_mat
    #     num_indices- number of points in top most row of distance_matrix

    n_grid_indices = 0

    grid_indices = np.concatenate([grid_indices, np.zeros([hor,1])])
    grid_indices = np.reshape(grid_indices, [vert+1, 2*hor+1])

    keep_ind = np.reshape(range(0,np.size(grid_indices)), [vert+1, 2*hor+1])

    inc_col = range(0, 2*hor+1)
    zer = hor

    if j < hor:
        for k in range(0, hor-j):
            inc_col.pop(0)

    if j > N-hor-1:
        for k in range(N, j+hor+1):
            inc_col.pop(-1)
            zer -= 1
    diff = num - (vert)*N-hor

    inc_grid_indices = grid_indices[:,inc_col].flatten()
    keep_ind = keep_ind[:,inc_col].flatten()
    
    while zer != 0:
        inc_grid_indices = inc_grid_indices[0:-1]
        keep_ind = keep_ind[0:-1]
        zer -= 1

    inc_grid_indices = [int(ind) + diff for ind in inc_grid_indices]
    keep_ind = [int(ind) for ind in keep_ind]

    dist_mat = distance_matrix[:,keep_ind]
    dist_mat = dist_mat[keep_ind,:]

    return {'dist_mat':dist_mat, 'inc_indices':inc_grid_indices}

def calc_full_dist(vert, hor, N, site_collection_SM):
    #####
    # Calculates full distance matrix. Called once per row.
    # IN: vert- number of included rows
    #     hor- number of columns within radius 
    #     N- number of points in row
    #     site_collection_SM- site collection for ShakeMap data
    #OUT: grid_indices- indices of points included in distance matrix
    #     distance_matrix- full distance matrix
    #####

    # gathers indices for full distance matrix for each row
    grid_indices = [None]*((vert+1)*(2*hor+1))
    n_grid_indices = 0
    for k in range(0, vert+1):
        if k == vert:
            for j in range(0,hor+1):
                grid_indices[n_grid_indices] = j + N*k
                n_grid_indices += 1
        else:
            for j in range(0,2*hor+1):
                grid_indices[n_grid_indices] = j + N*k
                n_grid_indices += 1 
    del grid_indices[n_grid_indices:]

    distance_matrix = np.zeros([np.size(grid_indices), np.size(grid_indices)])
    grid_indices = np.vstack(grid_indices)
    
    # Create full distance matrix for row
    for k in range(0, np.size(grid_indices)):
        distance_matrix[k, k:] = geodetic_distance(
                   site_collection_SM.lons[grid_indices[k ]], site_collection_SM.lats[grid_indices[k]],
                   site_collection_SM.lons[grid_indices[k:]], site_collection_SM.lats[grid_indices[k:]]).flatten()
    
    distance_matrix = distance_matrix + distance_matrix.T

    return {'grid_indices':grid_indices, 'distance_matrix':distance_matrix}

def calc_vert_hor(i, r, l, d):
    ######
    # Calculates the number of vertical of horozontal points in the full distance matrix
    # IN: i- current row
    #     r- radius
    #     l,d- vectors of distances between points vertically and horozontally
    #OUT: vert- number of rows in radius above current point
    #     hor- number of columns in radius to the left of current point
    #     added_vert- number of rows in between first row and the first row included in radius
    ######
    hor = int(np.floor(r/l[1]))
    
    # number of data points in vertical direction
    k = i
    vert_r = 0
    vert = 0
    while vert_r+d[k] <= r:
        vert_r += d[k]
        k -= 1
        vert += 1
    
    # adjusts for the unincluded rows
    if i < vert:
        vert = i

    return {'vert':vert, 'hor':hor}

def set_up_grid_dist(M,N, site_collection_SM):
    ######
    # Calculates the vertical and horozontal spacing between points for each row
    # IN: M,N- number of points in grid vertically and horozontally
    #     site_collection_SM- site collection for ShakeMap data
    # OUT: l,d- vectors of distances between points vertically and horozontally
    ######

    l = np.zeros([M-1])
    d = np.zeros([M])
    
    # Calculate vertical and horozonal spacing between points for each row
    l[:] = geodetic_distance(site_collection_SM.lons[range(N,N*M,N)],   site_collection_SM.lats[range(N,N*M,N)],
                             site_collection_SM.lons[range(0,N*M-N,N)], site_collection_SM.lats[range(0,N*M-N,N)])
    d[:] = geodetic_distance(site_collection_SM.lons[range(0, M*N, N)], site_collection_SM.lats[range(0, M*N, N)],
                             site_collection_SM.lons[range(1, M*N, N)], site_collection_SM.lats[range(1, M*N, N)])
    return {'l':l, 'd':d}
