import numpy as np

def realizations(num_realizations, N,M, grid_arr, mu_arr, sigma_arr, uncertaintydata, DATA):
    for j in range(0, num_realizations):

        X = np.zeros([M*N,1])
        rand_arr = np.random.randn(M*N)
        X[0] = rand_arr[0]
        for i in range(1,M*N):
            nzeros = np.size(mu_arr[i]) - np.size(grid_arr[i])
            x = np.append(np.zeros(nzeros), X[np.array(grid_arr[i]).squeeze()])
            mu = np.dot(mu_arr[i], x - np.mean(X[0:i]))
            X[i] = mu + rand_arr[i] * sigma_arr[i]

        COR = np.reshape(X, [M,N])
        X = np.multiply(COR, uncertaintydata)
        DATA_NEW = DATA * np.exp(X)
        #    pow2 = np.sum(DATA_NEWnew.flatten())
        #    DATA_NEWnew = DATA_NEWnew * pow1 / pow2
        if j == 0:
            ACCUM_ARRAY = DATA_NEW.copy()
        else:
            ACCUM_ARRAY += DATA_NEW
        if num_realizations < 10:
            print "Done with", j+1, "of", num_realizations, "iterations."
        else:
            if np.mod(j+1, 25) == 0:
                print "Done with", j+1, "of", num_realizations, "iterations."

    ACCUM_ARRAY = ACCUM_ARRAY / num_realizations

    return ACCUM_ARRAY
