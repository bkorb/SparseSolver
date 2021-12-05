import numpy as np
import scipy as sp
from scipy import sparse
import matplotlib.pyplot as plt
import math as mt

# Thomas Alogrithm

def Thomas(A, b):

    n = len(A)
   
    r1_prime = A[0, :] / A[0,0]

    b1_prime = b[0] / A[0, 0]

    A[0, :] = r1_prime
    b[0] = b1_prime


    for i in range(1,n):
        # A[i, :] = A[i, i-1] - A[i, i-1]*A[i-1, :]
        #b[i] = A[i, i-1] - A[i, i-1]*b[i-1]

        



        
        
        A[i, :] = A[i, :] - A[i,i-1] * A[i-1, :] #R_i = R_i - b_i*R_i-1

        A[i, :] = A[i, :] /( A[i, i] - A[i-1, i] * A[i, i-1] / A[i, i]  ) #R_i = R_i / (a_i - (c_{i-2} * b_i / a_i))


        b[i] = b[i] - A[i, i-1] * b[i-1] #Performs same row operations to the b vector

        b[i] = b[i]/ ( A[i, i] - A[i-1, i] * A[i, i-1] / A[i, i]  )
        


   

    return(A,b)


    


def tridiagonal(d, u, l):

    diagonals = [d, u, l]
    B = sp.sparse.diags(diagonals, [0, 1, -1]).toarray() #Creates a tridiagonal matrix with vector d on the main diagonal, vector u on the upper diagonal, and l on the lower diagonal

    return B




# Cubic Splines

N = 5
x = np.linspace(-5, 5, N)


d = np.array([((x[i+2] - x[i+1]) + (x[i+1] - x[i])) / 3 for i in range(0, N-2)]) 
u = np.array([(x[i+2] - x[i+1]) / 6 for i in range(0, N-2)])
l = np.array([(x[i+1] - x[i])/6 for i in range(0, N-2)])

A = tridiagonal(d, u, l)


#f = lambda x: 1 / (1+x**2) #Function we are interpolating with cubic splines
f = lambda x: x**3 - 8*x

#b = np.array([((f(x[i+2]) - f(x[i+1]))/ (x[i+2] - x[i+1])) - (((f(x[i+1]) - f(x[i]))/ (x[i+1] - x[i]))) for i in range(0,N-2)]) 
b = np.zeros((N-2, 1))
for i in range(0, N-2):
    b[i] = ((f(x[i+2]) - f(x[i+1]))/ (x[i+2] - x[i+1])) - ((f(x[i+1]) - f(x[i]))/ (x[i+1] - x[i]))


A, b = Thomas(A,b)

print(A)
print(b)

#n = len(A)


M = np.dot(np.linalg.inv(A), b)
#for i in range(n-1, 0, -1):


M1 = np.append(M, 0)
M2 = np.insert(M1, 0, 0)
print(M2)

# Building the splines
for i in range(0, N-1):
    x_fine = np.linspace(x[i],x[i+1],101)
    
    C = (f(x[i])/(x[i+1] - x[i])) -  (x[i+1] - x[i])*M2[i] / 6
    D = (f(x[i+1])/(x[i+1] - x[i])) -  (x[i+1] - x[i])*M2[i+1] / 6
        
    S = ((x[i+1] - x_fine)**3) *M2[i] / (6 * (x[i+1] - x[i])) + ((x_fine - x[i])**3)*M2[i+1] / (6 * (x[i+1] - x[i])) + C*(x[i+1] - x_fine) + D*(x_fine - x[i])
    
    plt.plot(x_fine, S)
    
   

plt.show()
