import sys
import scipy as scipy
import math

# Counts critical points using the multinomial coefficient.
# The input is a class of critical points, the output is the number of critical
# points in that class along with the energy value associated to that class in
# terms of alpha, beta, gamma, and delta.
# The greek variables correspond to the energy values of the critcal points of
# internal bond type 2-2 in increasing order. 

def cps1321(n,k,i1,i2,j1,j2):
    cp = (math.factorial(n)/(math.factorial(n-k-j1-j2)*math.factorial(j1)
        *math.factorial(j2)*math.factorial(i1)*math.factorial(i2)
        *math.factorial(k-i1-i2)))

    print('There are ',cp,' cps corresponding to energy value ',  j1, '\u03B1 +'
            , j2, '\u03B1 hat+', n-k-j1-j2, '\u03B2+', k-i1-i2, '\u03B3+'
            , i1, '\u03B4+', i2, '\u03B4 hat.')

if __name__ == "__main__":
    if len(sys.argv) > 1:
        cps1321(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),
                int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]))
    else:
        raise SystemExit("usage:  python cps1321.py n k i1 i2 j1 j2 ")
