sourceCpp("test.cpp")
mat <- matrix(1:12, ncol=4)
mat
res <- cppf(mat)
mat