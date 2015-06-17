import numpy as np
from scipy import linalg
import math
import sys
import random
from neicio.readstation import readStation
from neicio.shake import ShakeGrid
from mpl_toolkits.mplot3d import Axes3D
from openquake.hazardlib.correlation import JB2009CorrelationModel
from openquake.hazardlib.correlation import BaseCorrelationModel
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo import Point
from openquake.hazardlib.geo.geodetic import geodetic_distance
from openquake.hazardlib.imt import from_string
import time
from matplotlib import cm
from neicio.gmt import GMTGrid


def main(var, r, voi, rand, vscorr, intensity_factor):

    np.set_printoptions(linewidth = 200)
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
    
    # Get spcing of horozontal and vertical points
    ld  = set_up_grid_dist(M,N,var['site_collection_SM'])

    full_dist_arr = [None]*M
    full_grid_indices = [None]*M
    vhva = [None]*M
        
    for i in range(0,M):
        # Find the number of points in radius horozontally and vertically for each row
        vhva[i] = calc_vert_hor(i, r, ld['l'], ld['d'], M)
        
       # Calculate the full distance matrix for each row
        dist = calc_full_dist(vhva[i]['vert_u'], vhva[i]['vert_d'], vhva[i]['hor'], N, var['site_collection_SM'])

        full_dist_arr[i] = dist['distance_matrix']
        full_grid_indices[i] = dist['grid_indices']

    order = range(0,M*N)
    random.shuffle(order)

    for i in range(0,M*N):
        num = order[i]
        row = int(np.floor(num/N))
        col = int(np.mod(num,N))

        # Find the reduced distance matrix 
        dist_calc = reduce_distance(col, vhva[row]['vert_u'], vhva[row]['vert_d'], vhva[row]['hor'],
                                    N, full_dist_arr[row], full_grid_indices[row], num, order, i, X)

        # Include stations in distance matrix and find the indices of the points within the radius
        out = inc_stations(col, row, N, K, r, var['site_collection_SM'], var['site_collection_station'], 
                           dist_calc['dist_mat'], X, dist_calc['inc_indices'])

        if np.size(dist_calc['inc_indices']) == 1:
        
            # no conditioning points, correlation value is random
            X[num] = rand[num]

            # Store for multiple realizations
            grid_arr [num] = np.zeros(0)
            mu_arr   [num] = np.zeros(0)
            sigma_arr[num] = 1

        else:

            other = calculate_corr(out['dist_mat'], voi, JB_cor_model, var, out['inc_sta_indices'], intensity_factor)
            mu = other['Sig12'].T*other['Sig11inv']*out['x']
            rand_num = rand[num]
            X[num] = mu+rand_num*other['R']
            # Store for multiple realizations                                                                             
            grid_arr [num] = dist_calc['inc_indices'][0:-1]
            mu_arr   [num] = other['Sig12'].T*other['Sig11inv']
            sigma_arr[num] = other['R']

        if np.mod(i, 5000) == 0:
            print 'Finishing step:', i
            sys.stdout.flush()

    DATA = var['data']
    COR = np.reshape(X, [M,N])

    #Multiply by uncertainty                                                                                                   
    X = np.multiply(COR, var['uncertaintydata'])
                
    DATA_NEW = np.log(DATA) + X
    DATA_NEW = np.exp(DATA_NEW)

    return {'cor':COR, 'data':DATA, 'data_new':DATA_NEW, 'grid_arr':grid_arr, 'mu_arr':mu_arr, 'sigma_arr':sigma_arr}


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

def reduce_distance(j, vert_u, vert_d, hor, N, distance_matrix, grid_indices, num, order, i, X):
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
    
    grid_indices = np.reshape(grid_indices, [vert_u + vert_d+1, 2*hor+1])
    keep_ind = np.reshape(range(0,np.size(grid_indices)), [vert_u + vert_d+1, 2*hor+1])

    inc_col = range(0,2*hor+1)

    if j < hor:
        for k in range(0,hor-j):
            inc_col.pop(0)

    if j > N-hor-1:
        for k in range(N,j+hor+1):
            inc_col.pop(-1)

    diff = num - (vert_u)*N-hor
    
    inc_grid_indices = grid_indices[:,inc_col].flatten()
    inc_grid_indices = [ind + diff for ind in inc_grid_indices]

    cl = inc_grid_indices.index(num)

    inc_grid_indices.pop(cl)
    inc_grid_indices.append(num)

    keep_ind = keep_ind[:,inc_col].flatten()
    keep_ind = [ind for ind in keep_ind]

    a = keep_ind.pop(cl)
    keep_ind.append(a)
    
    inc_indices = inc_grid_indices[:]
    
    for i in range(np.size(inc_indices)-2, -1,-1):
        index = X[inc_indices[i]]
        if index == 0:
            inc_indices.pop(i)
            keep_ind.pop(i)

    dist_mat = distance_matrix[:,keep_ind]
    dist_mat = dist_mat[keep_ind,:]

    return {'dist_mat':dist_mat, 'inc_indices':inc_indices}

 
def calc_full_dist(vert_u, vert_d, hor, N, site_collection_SM):
    #####
    # Calculates full distance matrix. Called once per row.
    # IN: vert_u- number of included rows above point
    #     vert_d- number of included rows below point
    #     hor- number of columns within radius 
    #     N- number of points in row
    #     site_collection_SM- site collection for ShakeMap data
    #OUT: grid_indices- indices of points included in distance matrix
    #     distance_matrix- full distance matrix
    #####

    # gathers indices for full distance matrix for each row
    grid_indices = [None]*((vert_u+vert_d+1)*(2*hor+1))
    n_grid_indices = 0
    
    for k in range(0, vert_u+vert_d+1):
        for j in range(0,2*hor+1):
            grid_indices[n_grid_indices] = j + N*k
            n_grid_indices +=1
                
    del grid_indices[n_grid_indices:]

    mod_grid_indices = grid_indices[:]
    mod_grid_indices.remove(N*vert_u + hor)
    mod_grid_indices.append(N*vert_u + hor)
    
    distance_matrix = np.zeros([np.size(mod_grid_indices), np.size(mod_grid_indices)])
    grid_indices = np.vstack(grid_indices)
    
    # Create full distance matrix for row
    for k in range(0, np.size(grid_indices)):
        distance_matrix[k, k:] = geodetic_distance(
                   site_collection_SM.lons[grid_indices[k ]], site_collection_SM.lats[grid_indices[k]],
                   site_collection_SM.lons[grid_indices[k:]], site_collection_SM.lats[grid_indices[k:]]).flatten()
    
    distance_matrix = distance_matrix + distance_matrix.T
    
    return {'grid_indices':grid_indices, 'distance_matrix':distance_matrix}

def calc_vert_hor(i, r, l, d, M):
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

    k = i
    vert = 0
    vert_r = 0
    while vert_r + d[k] <= r:
        vert_r += d[k]
        k -= 1
        vert += 1

    if vert > i:
        vert_u = i
    else:
        vert_u = vert

    if vert > M-i-1:
        vert_d = M-i-1
    else:
        vert_d = vert

    return {'vert_u':vert_u, 'vert_d': vert_d, 'hor':hor}

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
