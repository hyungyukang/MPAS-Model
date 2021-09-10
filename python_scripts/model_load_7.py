import numpy as np
from numpy import savetxt
import os
import sys
import keras
from keras import models
from keras import layers
import matplotlib.pyplot as plt
from sklearn import preprocessing



ncpu = 32 

nx = 220
nz = 7
ny = 220
ncpu_train = 24
ncpu_test = 4
ncpu_val = 4

x_train = np.zeros([ncpu_train*nx,nz])
x_test  = np.zeros([ncpu_test*nx,nz])
x_val   = np.zeros([ncpu_val*nx,nz])

y_train = np.zeros([ncpu_train*nx,ny])
y_test  = np.zeros([ncpu_test*nx,ny])
y_val   = np.zeros([ncpu_val*nx,ny])


isum1 = 0
isum2 = 0
isum3 = 0

for n in range(ncpu):
    afile = "{0}{1:03d}{2}".format("data_new_1/mata_",n,".bin")
    bfile = "{0}{1:03d}{2}".format("data_new_1/matb_",n,".bin")

    mata = np.fromfile(afile)
    matb = np.fromfile(bfile)

    if ( n < ncpu_train ):
       for j in range(nx):
           x_train[isum1,:] = mata[j*nx:j*nx+nz]
           y_train[isum1,:] = matb[j*nx:j*nx+ny]
           isum1 = isum1 + 1
    elif ( n > ncpu_train-1 and n < ncpu_train+ncpu_test ):
       for j in range(nx):
           x_test[isum2,:] = mata[j*nx:j*nx+nz]
           y_test[isum2,:] = matb[j*nx:j*nx+ny]
           isum2 = isum2 + 1
    else:
       for j in range(nx):
           x_val[isum3,:] = mata[j*nx:j*nx+nz]
           y_val[isum3,:] = matb[j*nx:j*nx+ny]
           isum3 = isum3 + 1


#print(x_train[49,:])
#print(y_train[49,:])

# Y normalization ----------------------------------------------
y_train_min = np.min(y_train)
y_test_min = np.min(y_test)
y_val_min = np.min(y_val)
y_min = min(y_train_min,y_test_min,y_val_min)



y_train_max = np.max(y_train)
y_test_max = np.max(y_test)
y_val_max = np.max(y_val)
y_max = max(y_train_max,y_test_max,y_val_max) + (-y_min*2.0)

print(y_min,y_max)

y_train_max = np.max(abs(y_train))
y_test_max = np.max(abs(y_test))
y_val_max = np.max(abs(y_val))
y_max = max(y_train_max,y_test_max,y_val_max)



#y_train = y_train + abs(y_min*2.0)
#y_train = y_train / y_max
#y_test  = y_test  + abs(y_min*2.0)
#y_test  = y_test  / y_max
#y_val   = y_val   + abs(y_min*2.0)
#y_val   = y_val   / y_max


y_train = abs(y_train) / y_max
y_test  = abs(y_test)  / y_max
y_val   = abs(y_val)   / y_max

omit_value = abs(y_min*2.0) / y_max

for j in range(ncpu_train*nx):
    for i in range(ny):
        if ( y_train[j,i] == omit_value ):
           y_train[j,i] = 0.0 


# X normalization ----------------------------------------------
x_train = -x_train
x_test = -x_test
x_val = -x_val

x_train_min = np.min(x_train)
x_test_min = np.min(x_test)
x_val_min = np.min(x_val)
x_min = min(x_train_min,x_test_min,x_val_min)

x_train_max = np.max(x_train)
x_test_max = np.max(x_test)
x_val_max = np.max(x_val)
x_max = max(x_train_max,x_test_max,x_val_max) + (x_min*2.0)


x_train_max = np.max(abs(x_train))
x_test_max = np.max(abs(x_test))
x_val_max = np.max(abs(x_val))
x_max = max(x_train_max,x_test_max,x_val_max)


#x_train = x_train + abs(x_min*2.0)
#x_train = x_train / x_max
#x_test  = x_test  + abs(x_min*2.0)
#x_test  = x_test  / x_max
#x_val   = x_val   + abs(x_min*2.0)
#x_val   = x_val   / x_max


x_train = abs(x_train) / x_max
x_test  = abs(x_test)  / x_max
x_val   = abs(x_val)   / x_max


omit_value = abs(x_min*2.0) / x_max

for j in range(ncpu_train*nx):
    for i in range(nz):
        if ( x_train[j,i] == omit_value ):
           x_train[j,i] = 0.0 


#model = models.Sequential()
model = keras.models.load_model("my_model")


isum1 = 0
isum2 = 0
isum3 = 0


#predictions_train = model.predict(x_train)
#predictions_test = model.predict(x_test)
#predictions_val = model.predict(x_val)

for n in range(ncpu):
    afile = "{0}{1:03d}{2}".format("output_new_1/matp_",n,".csv")

    if ( n < ncpu_train ):
       predictions = model.predict(x_train[n*nx:(n+1)*nx])
       savetxt(afile,predictions,delimiter=' ')
    elif ( n > ncpu_train-1 and n < ncpu_train+ncpu_test ):
       predictions = model.predict(x_test[(n-ncpu_train)*nx:(n-ncpu_train+1)*nx])
       savetxt(afile,predictions,delimiter=' ')
    else:
       predictions = model.predict(x_val[(n-ncpu_train-ncpu_test)*nx:(n-ncpu_train-ncpu_test+1)*nx])
       savetxt(afile,predictions,delimiter=' ')

#predictions = model.predict(x_train)
#print(predictions[100,:]*y_max+y_min*2.0)
#xplot = predictions[1,:]*x_max+x_min*2.0

#c0 = model.layers[0].get_weights()
#w0 = c0[0]
#b0 = c0[1]
#c1 = model.layers[1].get_weights()
#w1 = c1[0]
#b1 = c1[1]
#c2 = model.layers[2].get_weights()
#w2 = c2[0]
#b2 = c2[1]
#c3 = model.layers[3].get_weights()
#w3 = c3[0]
#b3 = c3[1]
#c4 = model.layers[4].get_weights()
#w4 = c4[0]
#b4 = c4[1]
#
##for i in range(nx):
#y0 = np.matmul(y_val[0,:],w0)
#y0 = y0 + b0
#
##y1 = np.matmul(w1,y0)
#y1 = np.matmul(y0,w1)
#y1 = y1 + b1
##y2 = np.matmul(w2,y1)
#y2 = np.matmul(y1,w2)
#y2 = y2 + b2
##y3 = np.matmul(w3,y2)
#y3 = np.matmul(y2,w3)
#y3 = y3 + b3
##y4 = np.matmul(w4,y3)
#y4 = np.matmul(y3,w4)
#y4 = y4 + b4

#print(y4)

#print(w1[0].shape)
#print(w1[1].shape)

#aa = w1[0]

#print(aa[:,2])

#results = model.evaluate(y_val,x_val)
#print(results)

ii = 400

i = 0
js = 0
je = 3

predictions = model.predict(x_val)

xplot = predictions[i,js:je]
yplot = y_val[i,js:je]
plt.plot(yplot,'b',lw=2,label='true')
for i in range(ii):
    yplot = y_val[i,js:je]
    plt.plot(yplot,'b',lw=2)

i = 0
plt.plot(xplot,'r',label='predict')
for i in range(ii):
    xplot = predictions[i,js:je]
    plt.plot(xplot,'r')

#plt.legend()
#plt.show()

plt.legend()
#plt.show()
plt.savefig('test.png',dpi=150)

#print(xplot[0:10])
#print(yplot[0:10])
#print(xplot)
#print(yplot)
#print(y_max,y_min)
#print(y_max,y_min)

#plt.plot(epochs,loss,'bo',label='Traning acc')
#plt.plot(epochs,val_loss,'b',label='Validation acc')
#plt.xlabel('Epochs')
#plt.ylabel('Loss')
#plt.legend()
#plt.show()

#
#print(y_train[:,600])

#print(isum1)
#print(isum2)
#print(isum3)
#
#    for i in range(120):
#        print(mata[i*120])

#print(mata.shape)
#print(mata[0:122])
#print(matb[0:122])
