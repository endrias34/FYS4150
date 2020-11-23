

import numpy as np 
import matplotlib.pyplot as plt




L = [20, 30, 40, 50, 60, 70, 80, 90, 100]

parallel_if_averaging = [7.094106, 17.507081, 34.074822, 52.338674, 76.492429, 103.331404, 135.892967, 176.009939, 220.654341]
parallel_averaging = [6.865524, 15.818339, 27.414457, 45.793057, 70.254156, 99.516360, 139.731916, 177.306885, 225.547750]
parallel_locks = [13.121616, 29.486451, 54.196086, 83.236027, 120.018317, 174.713468, 239.429709, 306.842964, 391.434286]
parallel_non_lock = [6.081202, 14.576147, 25.981612, 42.052074, 61.284745, 82.962024, 115.292968, 153.536055, 199.991103]
one_thread = [19.835734, 45.516689, 81.656112, 125.133267, 182.664597, 244.037707, 326.672486, 414.778282, 531.840172]




from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
poly_model = make_pipeline(PolynomialFeatures(2),
                           LinearRegression())

xfit = np.linspace(20, 100, 1000)

poly_model.fit(np.array(L)[:, np.newaxis], np.array(parallel_if_averaging))
yfit = poly_model.predict(xfit[:, np.newaxis])
plt.scatter(np.array(L), np.array(parallel_if_averaging))
plt.plot(xfit, yfit, label='Parallel with if else averaging');




poly_model.fit(np.array(L)[:, np.newaxis], np.array(parallel_averaging))
yfit = poly_model.predict(xfit[:, np.newaxis])
plt.scatter(np.array(L), np.array(parallel_averaging))
plt.plot(xfit, yfit, label='Parallel with averaging');







poly_model.fit(np.array(L)[:, np.newaxis], np.array(parallel_non_lock))
yfit = poly_model.predict(xfit[:, np.newaxis])
plt.scatter(np.array(L), np.array(parallel_non_lock))
plt.plot(xfit, yfit, label='Parallel with bias');
plt.legend(loc='upper left', borderaxespad=0.)
plt.legend(loc='upper left', borderaxespad=0., fontsize=20)
plt.ylabel("Seconds", fontsize=20)
plt.xlabel("Lattice size", fontsize=20)

ax = plt.gca()
ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
ax.tick_params(axis = 'both', which = 'minor', labelsize = 18)

plt.show()





poly_model.fit(np.array(L)[:, np.newaxis], np.array(parallel_locks))
yfit = poly_model.predict(xfit[:, np.newaxis])
plt.scatter(np.array(L), np.array(parallel_locks), c="purple")
plt.plot(xfit, yfit, c="purple", label='Parallel with locks');




poly_model.fit(np.array(L)[:, np.newaxis], np.array(one_thread))
yfit = poly_model.predict(xfit[:, np.newaxis])
plt.scatter(np.array(L), np.array(one_thread), c="pink")
plt.plot(xfit, yfit, c="pink", label='One thread');
plt.legend(loc='upper left', borderaxespad=0., fontsize=20)
plt.ylabel("Seconds", fontsize=20)
plt.xlabel("Lattice size", fontsize=20)

ax = plt.gca()
ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
ax.tick_params(axis = 'both', which = 'minor', labelsize = 18)
plt.show()






############################# DATA #############################

'''
Run of method 3 with lattice size 100x100 and 1000000 MC samples with temperature 1.000000 took 225.547750 seconds

Run of method 3 with lattice size 90x90 and 1000000 MC samples with temperature 1.000000 took 177.306885 seconds

Run of method 3 with lattice size 80x80 and 1000000 MC samples with temperature 1.000000 took 139.731916 seconds

Run of method 3 with lattice size 70x70 and 1000000 MC samples with temperature 1.000000 took 99.516360 seconds

Run of method 3 with lattice size 60x60 and 1000000 MC samples with temperature 1.000000 took 70.254156 seconds

Run of method 3 with lattice size 50x50 and 1000000 MC samples with temperature 1.000000 took 45.793057 seconds

Run of method 3 with lattice size 40x40 and 1000000 MC samples with temperature 1.000000 took 27.414457 seconds

Run of method 3 with lattice size 30x30 and 1000000 MC samples with temperature 1.000000 took 15.818339 seconds

Run of method 3 with lattice size 20x20 and 1000000 MC samples with temperature 1.000000 took 6.865524 seconds 


'''

'''
Run of method 0 with lattice size 100x100 and 1000000 MC samples with temperature 1.000000 took 531.840172 seconds

Run of method 0 with lattice size 90x90 and 1000000 MC samples with temperature 1.000000 took 414.778282 seconds

Run of method 0 with lattice size 80x80 and 1000000 MC samples with temperature 1.000000 took 326.672486 seconds

Run of method 0 with lattice size 70x70 and 1000000 MC samples with temperature 1.000000 took 244.037707 seconds

Run of method 0 with lattice size 60x60 and 1000000 MC samples with temperature 1.000000 took 182.664597 seconds 

Run of method 0 with lattice size 50x50 and 1000000 MC samples with temperature 1.000000 took 125.133267 seconds

Run of method 0 with lattice size 40x40 and 1000000 MC samples with temperature 1.000000 took 81.656112 seconds

Run of method 0 with lattice size 30x30 and 1000000 MC samples with temperature 1.000000 took 45.516689 seconds

Run of method 0 with lattice size 20x20 and 1000000 MC samples with temperature 1.000000 took 19.835734 seconds

'''

'''

Run of method 1 with lattice size 100x100 and 1000000 MC samples with temperature 1.000000 took 391.434286 seconds

Run of method 1 with lattice size 90x90 and 1000000 MC samples with temperature 1.000000 took 306.842964 seconds

Run of method 1 with lattice size 80x80 and 1000000 MC samples with temperature 1.000000 took 239.429709 seconds

Run of method 1 with lattice size 70x70 and 1000000 MC samples with temperature 1.000000 took 174.713468 seconds

Run of method 1 with lattice size 60x60 and 1000000 MC samples with temperature 1.000000 took 120.018317 seconds

Run of method 1 with lattice size 50x50 and 1000000 MC samples with temperature 1.000000 took 83.236027 seconds

Run of method 1 with lattice size 40x40 and 1000000 MC samples with temperature 1.000000 took 54.196086 seconds

Run of method 1 with lattice size 30x30 and 1000000 MC samples with temperature 1.000000 took 29.486451 seconds

Run of method 1 with lattice size 20x20 and 1000000 MC samples with temperature 1.000000 took 13.121616 seconds

'''


'''

Run of method 2 with lattice size 100x100 and 1000000 MC samples with temperature 1.000000 took 199.991103 seconds

Run of method 2 with lattice size 90x90 and 1000000 MC samples with temperature 1.000000 took 153.536055 seconds

Run of method 2 with lattice size 80x80 and 1000000 MC samples with temperature 1.000000 took 115.292968 seconds

Run of method 2 with lattice size 70x70 and 1000000 MC samples with temperature 1.000000 took 82.962024 seconds

Run of method 2 with lattice size 60x60 and 1000000 MC samples with temperature 1.000000 took 61.284745 seconds

Run of method 2 with lattice size 50x50 and 1000000 MC samples with temperature 1.000000 took 42.052074 seconds

Run of method 2 with lattice size 40x40 and 1000000 MC samples with temperature 1.000000 took 25.981612 seconds

Run of method 2 with lattice size 30x30 and 1000000 MC samples with temperature 1.000000 took 14.576147 seconds

Run of method 2 with lattice size 20x20 and 1000000 MC samples with temperature 1.000000 took 6.081202 seconds

'''

# IFFF EEEEELLLSE
'''
Run of method 3 with lattice size 100x100 and 1000000 MC samples with temperature 1.000000 took 220.654341 seconds

Run of method 3 with lattice size 90x90 and 1000000 MC samples with temperature 1.000000 took 176.009939 seconds

Run of method 3 with lattice size 80x80 and 1000000 MC samples with temperature 1.000000 took 135.892967 seconds

Run of method 3 with lattice size 70x70 and 1000000 MC samples with temperature 1.000000 took 103.331404 seconds

Run of method 3 with lattice size 60x60 and 1000000 MC samples with temperature 1.000000 took 76.492429 seconds

Run of method 3 with lattice size 50x50 and 1000000 MC samples with temperature 1.000000 took 52.338674 seconds

Run of method 3 with lattice size 40x40 and 1000000 MC samples with temperature 1.000000 took 34.074822 seconds

Run of method 3 with lattice size 30x30 and 1000000 MC samples with temperature 1.000000 took 17.507081 seconds

Run of method 3 with lattice size 20x20 and 1000000 MC samples with temperature 1.000000 took 7.094106 seconds

'''