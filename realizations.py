import numpy as np

def realizations(num_realizations, radius, N,M, grid_arr, mu_arr, sigma_arr, uncertaintydata, DATA):

    if num_realizations == 0:
        return
    else:
        h = open('random_arrays_LP_100', 'r')
        
        g = open(('workfile_LP_radius_R%i_100' % radius), 'w')

        s = str(num_realizations)
        g.write(s+'\n')

        rand_arr = np.zeros(M*N)

        for j in range(0, num_realizations):

            for i in range(0,M*N):
                rand_arr[i] = h.readline()

            X = np.zeros([M*N,1])

            for i in range(0,M*N):
                nzeros = np.size(mu_arr[i]) - np.size(grid_arr[i])
                x = np.append(np.zeros(nzeros), X[np.array(grid_arr[i], dtype = 'i').squeeze()])
                mu = np.dot(mu_arr[i], x)
                X[i] = mu + rand_arr[i] * sigma_arr[i]
                s = str(X[i])
                g.write(s+'\n')


#        COR = np.reshape(X, [M,N])
#        X = np.multiply(COR, uncertaintydata)
#        DATA_NEW = DATA * np.exp(X)

#        if j == 0:
#            ACCUM_ARRAY = DATA_NEW.copy()
#        else:
#            ACCUM_ARRAY += DATA_NEW
        
            if np.mod(j+1, 25) == 0:
                print "Done with", j+1, "of", num_realizations, "iterations."

        g.close()
        h.close()

        return
