import os
def mkwdir(where='./output_data'):
    counter = 0
    while os.path.isdir(where + '/wdir' + str(counter)): #improvement: use os.path.join()
        counter += 1
    path = where + '/wdir' + str(counter)
    os.makedirs(where + '/wdir' + str(counter))
    return path
    

def swim(wdir, iteration, action, nswim=100):
    os.system('./runs/bin/swim -dir ' + wdir + ' -it ' + str(iteration) + ' -a ' + str(action) + '-nsw' + str(nswim))

wdir = mkwdir()

for i in range(10000):
    swim(wdir, i, 6)