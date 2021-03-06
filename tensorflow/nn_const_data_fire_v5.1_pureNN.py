import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'
import tensorflow.compat.v1 as tf
#from tensorflow.data import Dataset, Iterator
import numpy as np
#from numpy import genfromtxt
import sys
import argparse
import pandas as pd

tf.disable_v2_behavior()

#string regions_names[] = {"BONA (Boreal North America)",				// 1
#				          "TENA (Temperate North America)",				// 2
#				          "CEAM (Central America)",						// 3
#				          "NHSA (Northern Hemisphere South America)",	// 4
#				          "SHSA (Southern Hemisphere South America)",	// 5
#				          "EURO (Europe)",								// 6
#				          "MIDE (Middle East)",							// 7
#				          "NHAF (Northern Hemisphere Africa)",			// 8
#				          "SHAF (Southern Hemisphere Africa)",			// 9
#				          "BOAS (Boreal Asia)",							// 10
#				          "CEAS (Central Asia)",						// 11
#				          "SEAS (Southeast Asia)",						// 12
#				          "EQAS (Equatorial Asia)",						// 13
#				          "AUST (Australia and New Zealand)"};			// 14


#r1 = 6
#r2 = 7

parser = argparse.ArgumentParser()
parser.add_argument("-d", action="store", dest="output_dir", default="output")
parser.add_argument("-m", action="store", dest="model_dir", default="mod1")
result = parser.parse_args()

output_dir = result.output_dir  #'output_globe'
model_dir  = result.model_dir   #'output_globe'
print("output_dir = ", output_dir)
print("model_dir  = ", model_dir)

__learn_rate = 0.005
__batch_size = 5000
__n_steps = 2000

n_classes = 25

nn_h1 = 12

import sys
sys.path.insert(0, '../'+output_dir)
sys.path.insert(0, '../'+output_dir+'/'+model_dir)

from regions import *
from variables import *	# Import variable IDs as in training data csv

print("regions = ", r1, " ", r2)

ID_ft = list(range(ID_ftmap1, ID_ftmap11)) #+ [ID_ftmap12]	# this will exclude ftmap11 (cropland) and Barren (0)

X_ids = [ID_gpp, ID_gppl1, ID_ts, ID_cld, ID_vp,] + ID_ft

Y_id = ID_gfedclass
	
# functions to initialize weights and biases
def weight_variable(shape):
  initial = tf.truncated_normal(stddev=0.5, shape=shape)
  return tf.Variable(initial)

def bias_variable(shape):
  initial = tf.truncated_normal(stddev=0.5, shape=shape)
  return tf.Variable(initial)


#def create_dataset(filename, map_fun, batch_size, rep=1, buffer_size=0):
#  dat = tf.data.TextLineDataset(filename)
#  dat = dat.skip(1)
#  dat = dat.map(map_fun) 
#  dat = dat.repeat(rep)
#  if (buffer_size>0): 
#    dat = dat.shuffle(buffer_size)  # for each iteration, refills the buffer with new data and chooses random elements to put into batch <-- this behaviour is a bit counterintuitive
#  dat = dat.batch(batch_size)  
#  return dat  



# ~~~~~~ forward prop ~~~~~~~~~~
def denseNet(x, W1,b1,Wo,bo):
  y1 = tf.nn.elu(tf.matmul(x,W1) + b1)  # first layer neurons with sigmoid activation
#  y2 = tf.nn.sigmoid(tf.matmul(y1,W2) + b2)  # first layer neurons with sigmoid activation
  y = tf.matmul(y1,Wo) + bo
  
  return y



### PREPARE TRAINING DATA AS NUMPY ARRAYS ###

print("Reading training data...")
df = pd.read_csv('../'+output_dir+'/train_forest.csv', engine='c')
my_data = df.values[(df.values[:,ID_region] == r1) | (df.values[:,ID_region] == r2)]  # Region selection
print("DONE")
np.set_printoptions(precision=6, suppress=True)
print("--------------")
print("Input: (head) | y");
print(my_data[0:5, X_ids+[Y_id]].astype(float))
print("--------------")
print("Output: (head)");
print(tf.Session().run(tf.one_hot(my_data[0:5,Y_id], depth=n_classes, dtype=tf.int32)))
print("--------------")




Xmeans = np.mean(my_data[:,X_ids].astype(float), axis=0)
Xstd = np.std(my_data[:,X_ids].astype(float), axis=0)+1e-8
print("Input means/sd:");
print(Xmeans)
print("--------------")

#Xmeans = np.concatenate((Xmeans, [0]*len(ID_ft)), axis=0)
#Xstd   = np.concatenate((Xstd,   [1]*len(ID_ft)), axis=0)
#X_ids   = X_ids   + ID_ft

print("Reading evalutation data...")
dfe = pd.read_csv('../'+output_dir+'/eval_forest.csv', engine='c')
eval_data = dfe.values[(dfe.values[:,ID_region] == r1) | (dfe.values[:,ID_region] == r2)]  # Region selection
print(eval_data[0:5, X_ids+[Y_id]].astype(float))

print("Reading test data...")
dfte = pd.read_csv('../'+output_dir+'/test_forest.csv', engine='c')
test_data = dfte.values[(dfte.values[:,ID_region] == r1) | (dfte.values[:,ID_region] == r2)]  # Region selection
print(test_data[0:5, X_ids+[Y_id]].astype(float))
print("--------------")


# Layer 1 : with a single neuron
Wi = np.array(np.diag(1/Xstd), dtype='float32') # tf.constant([1/3/ts_std, 0,0, 0,1/3/dxl_std, 0, 0, 0, 1/3/lmois_std], shape=[3,3], dtype=tf.float32)
bi = np.matmul(Xmeans, Wi) #tf.constant([ts_mean, dxl_mean, lmois_mean], shape=[1,3], dtype=tf.float32)

print(Wi)	
print(bi)	

#### TENSORFLOW GRAPH BUILDING STARTS HERE ### 

xin = tf.matmul(tf.cast(my_data[:,X_ids].astype(float), tf.float32), Wi) - bi
yin = tf.one_hot(my_data[:,Y_id].astype(float), depth=n_classes, dtype=tf.int32)

print("Scaled inputs:")
print(tf.Session().run(xin[0:5,:]))
print(tf.Session().run(yin[0:5,:]))
print("--------------")

xeval = tf.matmul(tf.cast(eval_data[:,X_ids].astype(float), tf.float32), Wi) - bi
yeval = tf.one_hot(eval_data[:,Y_id].astype(float), depth=n_classes, dtype=tf.int32)

xtest = tf.matmul(tf.cast(test_data[:,X_ids].astype(float), tf.float32), Wi) - bi
ytest = tf.one_hot(test_data[:,Y_id].astype(float), depth=n_classes, dtype=tf.int32)

print("Scaled test data:")
print(tf.Session().run(xtest[0:5,:]))
print("--------------")
	

dat_train = tf.data.Dataset.from_tensor_slices((xin,yin))
dat_train = dat_train.repeat(10000)
dat_train = dat_train.shuffle(1000000)
dat_train = dat_train.batch(__batch_size)

#dat_valid = tf.data.Dataset.from_tensor_slices((xeval,yeval))

## prepare data iterator
it_handle = tf.placeholder(tf.string, shape=[])

iterator = tf.data.Iterator.from_string_handle(it_handle, dat_train.output_types, dat_train.output_shapes)

next_batch = iterator.get_next()

training_iterator = dat_train.make_one_shot_iterator()
#validation_iterator = dat_valid.make_one_shot_iterator()

n_inputs = len(X_ids)


x  = tf.reshape(next_batch[0], shape=[-1,n_inputs])
y_ = tf.reshape(next_batch[1], shape=[-1,n_classes])


## NN DEFINITION ##

# Layer 1 : with a single neuron
W1 = weight_variable([n_inputs, nn_h1])	
b1 = bias_variable([nn_h1])

## Layer 2 : with a single neuron
#W2 = weight_variable([3,5])	
#b2 = bias_variable([5])

# output layer
Wo = weight_variable([nn_h1, n_classes])	
bo = bias_variable([n_classes])

# forward prop
y = denseNet(x, W1,b1,Wo,bo)

### TRAINING ###

# training operation
cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y))
train_op = tf.train.AdamOptimizer(__learn_rate).minimize(cross_entropy)

# calculate accuracy on the training sample
y_soft = tf.nn.softmax(denseNet(x, W1,b1,Wo,bo))
correct_prediction = tf.equal(tf.argmax(y_soft,1), tf.argmax(y_,1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))


# calculate accuracy on the validation dataset
xv  = tf.reshape(xeval, shape=[-1, n_inputs])
yv_ = tf.reshape(yeval, shape=[-1, n_classes])
yv_soft = tf.nn.softmax(denseNet(xv, W1,b1,Wo,bo))
correct_prediction_v = tf.equal(tf.argmax(yv_soft,1), tf.argmax(yv_,1))
accuracy_v = tf.reduce_mean(tf.cast(correct_prediction_v, tf.float32))

# calculate accuracy on the test dataset
xt  = tf.reshape(xtest, shape=[-1, n_inputs])
yt_ = tf.reshape(ytest, shape=[-1, n_classes])
yt_soft = tf.nn.softmax(denseNet(xt, W1,b1,Wo,bo))
correct_prediction_t = tf.equal(tf.argmax(yt_soft,1), tf.argmax(yt_,1))
accuracy_t = tf.reduce_mean(tf.cast(correct_prediction_t, tf.float32))


acc_avg = 0
accv_avg = 0
acct_avg = 0
count = 0
ce_avg = 0
with tf.Session() as sess:

  training_handle = sess.run(training_iterator.string_handle())
#  validation_handle = sess.run(validation_iterator.string_handle())

  sess.run(tf.global_variables_initializer())

  for i in range(__n_steps):
  #  iterator = dat.make_one_shot_iterator()
    try:
      _,acc,ce,accv,acct = sess.run([train_op, accuracy, cross_entropy, accuracy_v, accuracy_t], feed_dict={it_handle: training_handle})
#       print(sess.run(x,feed_dict={it_handle: training_handle}))

      if (i % 100 == 0):
        print("train step ",i,": ",ce, acc, "|", accv, "|", acct)
        
        if (i > 0.8*__n_steps):
          acc_avg = acc_avg +  acc
          accv_avg = accv_avg +  accv
          acct_avg = acct_avg +  acct
          ce_avg = ce_avg + ce
          count = count +1	
#        print(sess.run([accuracy, cross_entropy], feed_dict={it_handle: validation_handle}))
        
    except tf.errors.OutOfRangeError:
      print("End of file")
      break


  print("---------------------------------------------------------")
  print("Average    ",0,": ",ce_avg/count, "\t", acc_avg/count, "\t", accv_avg/count, "\t", acct_avg/count)

  # final class probabilities for train, eval and test datasets
  y_tr = sess.run(tf.nn.softmax(denseNet(tf.reshape(xin,   [-1,n_inputs]),W1,b1,Wo,bo)))
  y_ev = sess.run(tf.nn.softmax(denseNet(tf.reshape(xeval, [-1,n_inputs]),W1,b1,Wo,bo)))
  y_ts = sess.run(tf.nn.softmax(denseNet(tf.reshape(xtest, [-1,n_inputs]),W1,b1,Wo,bo)))


  np.savetxt("../"+output_dir+"/"+model_dir+"/y_predic_ba_train.txt",y_tr,delimiter=" ")
  np.savetxt("../"+output_dir+"/"+model_dir+"/y_predic_ba_eval.txt",y_ev,delimiter=" ")
  np.savetxt("../"+output_dir+"/"+model_dir+"/y_predic_ba_test.txt",y_ts,delimiter=" ")


#  def rmb(s):	# small function to remove square brackets from printed arrays 
#    dic = {'[':'', ']':''}
#    return "".join(dic.get(x,x) for x in str(s))
#    
#  class myFloat(float):
#    def __str__(self):
#      return "%.12f"%self
           
  np.set_printoptions(precision=10)
  
  orig_stdout = sys.stdout
  f = open("../"+output_dir+"/"+model_dir+"/weights_ba.txt",'w')
  sys.stdout = f

  print(1)
  print([n_inputs,nn_h1,n_classes])
  print("\n")
  
  print(Wi)	
  print("\n")
  print(-bi)
  print("\n\n")

  print(sess.run(W1))
  print("\n")
  print(sess.run(b1))
  print("\n\n")

  print(sess.run(Wo))
  print("\n")
  print(sess.run(bo))
  print("\n\n")

  sys.stdout=orig_stdout;  
  f.close()
  
  np.savetxt("../"+output_dir+"/"+model_dir+"/ce_and_accuracy.txt", [ce_avg/count, acc_avg/count,  accv_avg/count, acct_avg/count])  
  
  

