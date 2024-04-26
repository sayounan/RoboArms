Lab_04262024_notes

lda2.mlx
posture feautrure matrices (train)
8 EMGs, 4 feaures
32 x 32 matrix - covariance matrix (of feautures)
add up, devide by numPostures

within-class scatter matrix
diag = standard dev. 

10 eigenvals
10 eigenvecs

W matrix will be 32 x n (transformation matrix: control)
n=10
filter 
remove off
extract all feautures
project features on eigenvec for each posture


array of postures
each transformed pooint 
euclidean/cartesian (pythagrian math)
compute distance of each point to center point from each data cluster
classification of each EMG bin


confusion matrix

FInal Exam: presentation - story about process
filtering and data processing
off removal
feauture extraction 
W matrix
classification strategy


run everything on both training and validation data
then onece you get W, feed validation data in
create a function for 1 posture