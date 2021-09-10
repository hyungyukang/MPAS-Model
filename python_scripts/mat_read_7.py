import numpy as np
import os
import sys
import keras
from keras import models
from keras import layers
import matplotlib.pyplot as plt
from sklearn import preprocessing
import pandas as pd


#def clip_relu (x):
#    return keras.activations.relu(x,threshold=1.e-8)


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

print(x_train[0,:])
print(x_val[0,:])
print(x_val[0,:])

print(x_train[49,:])
print(y_train[49,:])

# Y normalization ----------------------------------------------
y_train_min = np.min(y_train)
y_test_min = np.min(y_test)
y_val_min = np.min(y_val)
y_min = min(y_train_min,y_test_min,y_val_min)

y_train_max = np.max(y_train)
y_test_max = np.max(y_test)
y_val_max = np.max(y_val)
y_max = max(y_train_max,y_test_max,y_val_max) + (-y_min*2.0)


y_train_max = np.max(abs(y_train))
y_test_max = np.max(abs(y_test))
y_val_max = np.max(abs(y_val))
y_max = max(y_train_max,y_test_max,y_val_max)

print('yval',y_min,y_max)

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


#infp = float('inf')
#infm = float('-inf')
#x_train = np.where(x_train == infp,0.0,x_train)
#x_train = np.where(x_train == infm,0.0,x_train)
#x_test = np.where(x_test == infp,0.0,x_test)
#x_test = np.where(x_test == infm,0.0,x_test)
#x_val = np.where(x_val == infp,0.0,x_val)
#x_val = np.where(x_val == infm,0.0,x_val)
#
#x_train_min = np.min(x_train)
#mp = 11
#x_train = np.where(x_train > mp,0.0,x_train)
#x_test = np.where(x_test > mp,0.0,x_test)
#x_val = np.where(x_val > mp,0.0,x_val)


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




print('xval',x_min,x_max)


#x_train = x_train + abs(x_min*2.0)
#x_train = x_train / x_max
#x_test  = x_test  + abs(x_min*2.0)
#x_test  = x_test  / x_max
#x_val   = x_val   + abs(x_min*2.0)
#x_val   = x_val   / x_max


x_train = abs(x_train) / x_max
x_test  = abs(x_test)  / x_max
x_val   = abs(x_val)   / x_max




#print(x_train[0,:])
#print(x_min,x_max)

omit_value = abs(x_min*2.0) / x_max

for j in range(ncpu_train*nx):
    for i in range(nz):
        if ( x_train[j,i] == omit_value ):
           x_train[j,i] = 0.0 






#y_train_max = np.max(y_train)
#print(np.min(x_train))
#print(np.min(y_train))
#print(y_min)

#x_train = preprocessing.normalize(x_train)
#print(x_train[:,1])

#y_train = preprocessing.normalize(y_train)
#print(y_train[200,:])

model = models.Sequential()
model.add(layers.Dense(nz,activation='linear',input_shape=(nz,)))
model.add(layers.Dropout(0.05))
model.add(layers.Dense(ny,activation='relu'))
model.add(layers.Dropout(0.05))
model.add(layers.Dense(ny,activation='linear'))
model.add(layers.Dense(ny,activation='relu'))
model.add(layers.Dropout(0.05))
model.add(layers.Dense(ny,activation='relu'))
model.add(layers.Dropout(0.05))
model.add(layers.Dense(ny,activation='linear'))
model.add(layers.Dense(ny,activation='softplus'))

#model.add(layers.Dense(nx,activation='linear',input_shape=(ny,)))
#model.add(layers.Dense(nx,activation=clip_relu))
#model.add(layers.Dense(nx*2,activation='relu'))
#model.add(layers.Dense(nx*2,activation='relu'))
#model.add(layers.Dense(nx*2,activation='linear'))
#model.add(layers.Dense(nx,activation='relu'))
#model.add(layers.Dense(nx,activation='linear'))
#model.add(layers.Dense(nx,activation='linear'))


#model.add(layers.Dense(nx,activation='linear',input_shape=(ny,)))
#model.add(layers.Dense(nx,activation='relu',min_value=0))
#model.add(layers.Dense(nx,activation='linear',min_value=0))

#model.add(layers.Dense(nx,activation='softmax'))
#model.add(layers.Dense(7,activation='relu'))
#model.add(layers.Dense(10,activation='relu'))

model.compile(optimizer='rmsprop',loss='mean_squared_error',metrics=['accuracy'])
history=model.fit(x_train,y_train,epochs=1000,batch_size=64,validation_data=(x_test,y_test))

model.save("my_model")

loss=history.history['accuracy']
val_loss=history.history['val_accuracy']
epochs=range(1,len(loss)+1)

#w1 = model.layers[0].get_weights()
#print(w1[0].shape)
#print(w1[1].shape)
#results = model.evaluate(y_val,x_val)
#print(results)

predictions = model.predict(x_val)

print(predictions[0])

xplot = predictions[i,:]
yplot = y_val[i,:]
plt.plot(yplot,'b',lw=4,label='true')
for i in range(200):
    xplot = predictions[i,:]
    yplot = y_val[i,:]
    plt.plot(yplot,'b',lw=4)

xplot = predictions[i,:]
yplot = y_val[i,:]
plt.plot(xplot,'r',label='predict')
for i in range(200):
    xplot = predictions[i,:]
    yplot = y_val[i,:]
    plt.plot(xplot,'r')

plt.legend()
#plt.show()
plt.savefig('test.png')

#plt.plot(epochs,loss,'bo',label='Traning acc')
#plt.plot(epochs,val_loss,'b',label='Validation acc')
#plt.xlabel('Epochs')
#plt.ylabel('Loss')
#plt.legend()
#plt.show()

#plt.plot(yplot,'b',lw=4,label='true')
#plt.plot(xplot,'r',label='predict')
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
