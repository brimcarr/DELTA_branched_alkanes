import sys
import scipy as scipy
import math

# Counts critical points using the multinomial coefficient

# n1 is the number of internal bonds of type 2-2
# n2 is the number of internal bonds of type 3-2

# The code iterates over all possible classes and counts the number of
# critical points in each class.
# Then, it prints the class (n1,n2,k1,k2,i11,j11,i21,i22,j21,j22),
# the number of points in each class, and the energy value associated to the
# class in terms of alphas, betas, gammas, and deltas.
# The greek variables correspond to the energy values of critical points in the
# internal bond energy landscapes. A subscript of 1 corresponds to the energy
# value of a critical point from the energy landscape of internal bond type 2-2
# and a subscript of 2 corresponds to the same for internal bond type 3-2.
def cps_1221_1321(n1,n2):
    avals=[]
    for k1 in range(0,n1+1):
        for i11 in range(0,k1+1):
            for j11 in range(0,n1-k1+1):
                for k2 in range(0,n2+1):
                    for i21 in range(0,k2+1):
                        for i22 in range(0,k2+1):
                            for j21 in range(0,n2-k2+1):
                                for j22 in range(0,n2-k2+1):
                                    if i21+i22<=k2 and j21+j22<=n2-k2 and i11<=k1 and j11<=n1-k1:
                                        cp = (2**(n1-i11-j11)*(math.factorial(n1)*math.factorial(n2))/(math.factorial(n1-k1-j11)
                                                *math.factorial(j11)*math.factorial(i11)*math.factorial(k1-i11)*math.factorial(n2-k2-j21-j22)*math.factorial(j21)
                                                *math.factorial(j22)*math.factorial(i21)*math.factorial(i22)*math.factorial(k2-i21-i22)))
                                        avals.append(cp)
                                        print('There are',cp,'cps of type',n1,n2,k1,k2,i11,j11,i21,i22,j21,j22,  'corresponding to energy value '
                                                , j11, '\u03B1_1 +', n1-k1-j11, '\u03B2_1+', k1-i11, '\u03B3_1+', i11, '\u03B4_1+'
                                                , j21, '\u03B1_2 +', j22, '\u03B1_2 hat+', n2-k2-j21-j22, '\u03B2_2+'
                                                , k2-i21-i22, '\u03B3_2+', i21, '\u03B4_2+', i22, '\u03B4_2 hat.')
    print(len(avals),sum(avals))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        cps_1221_1321(int(sys.argv[1]),int(sys.argv[2]))
    else:
        raise SystemExit("usage:  python cps1321.py n1 n2")
