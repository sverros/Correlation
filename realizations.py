import numpy as np

def realizations(num_realizations, N,M, grid_arr, mu_arr, sigma_arr, uncertaintydata, DATA):

    if num_realizations == 0:
        return
    else:
        g = open('workfileNR0_realizations_1000', 'w')

        s = str(num_realizations)
        g.write(s+'\n')

        for j in range(0, num_realizations):

            X = np.zeros([M*N,1])
            rand_arr = np.random.randn(M*N)
            X[0] = rand_arr[0]
            for i in range(1,M*N):
                nzeros = np.size(mu_arr[i]) - np.size(grid_arr[i])
                x = np.append(np.zeros(nzeros), X[np.array(grid_arr[i]).squeeze()])
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

        return
