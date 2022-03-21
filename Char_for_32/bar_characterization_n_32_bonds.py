import sys
import scipy as scipy
import math

# Counts number of bars for any number of copies of the bond 1-3-2-1.
# Input: n -- Number of bonds of type 1-3-2-1
# Output: Prints out the class, length of bars born, number of bars born,
# and the time of birth

def pc1321(n):
# Iterates through all classes of critical points that satisfy the following conditions:
# k<=n, i1+i2<=k, j1+j2<=n-k
    for k in range(0,n+1):
        for i1 in range(0,k+1):
            for i2 in range(0,k+1):
                for j1 in range(0,n-k+1):
                    for j2 in range(0,n-k+1):
                        if i1+i2<=k and j1+j2<=n-k:
                            avals=[]
                            bvals=[]
# Using other restrictions, the classes get sorted into four groups by bar length:
# gamma-alpha, delta-beta, semi infinite, and no bars born.
                            if j2 != 0: # gamma-alpha
                                for l in range(0,k-i1-i2+1):
                                    a=((-1)**(l)*math.factorial(n)
                                      /(math.factorial(n-k-j1-j2)*math.factorial(k-i1-i2-l)
                                      *math.factorial(i1)*math.factorial(i2)
                                      *math.factorial(j1)*math.factorial(j2+l)))
                                    avals.append(a)
                                print('Class',(n,k,i1,i2,j1,j2), u'\u03B3-\u03B1 bars born: ',
                                        (sum(avals)), ' at ', j1, '\u03B1 +'
                                        , j2, '\u03B1 hat +', n-k-j1-j2, '\u03B2 +'
                                        , k-i1-i2, '\u03B3 +', i1, '\u03B4 +', i2, '\u03B4 hat.')
                            elif i1+i2==k and j2==0 and j1<n-k: # delta-beta
                                for l in range(0,i1+1):
                                    b=((-1)**(l)*math.factorial(n)
                                      /(math.factorial(n-k-j1-j2+l)*math.factorial(k-i1-i2)
                                      *math.factorial(i1-l)*math.factorial(i2)
                                      *math.factorial(j2)*math.factorial(j1)))
                                    bvals.append(b)
                                print('Class',(n,k,i1,i2,j1,j2), u'\u03B4-\u03B2 bars born: '
                                        , (sum(bvals)), ' at ', j1, '\u03B1 +'
                                        , j2, '\u03B1 hat +', n-k-j1-j2, '\u03B2 +'
                                        , k-i1-i2, '\u03B3 +', i1, '\u03B4 +', i2, '\u03B4 hat.')
                            elif i1==0 and i2==k and j2==0 and j1==n-k: # semi-infinite
                                c=(math.factorial(n)
                                  /(math.factorial(n-k-j1-j2)*math.factorial(k-i1-i2)
                                  *math.factorial(i1)*math.factorial(i2)
                                  *math.factorial(j2)*math.factorial(j1)))
                                print('Class',(n,k,i1,i2,j1,j2), 'semi-\u221E bars born: '
                                        , c, ' at ', j1, '\u03B1 +', j2, '\u03B1 hat +'
                                        , n-k-j1-j2, '\u03B2 +', k-i1-i2, '\u03B3 +'
                                        , i1, '\u03B4 +', i2, '\u03B4 hat.')
                            else: # no bars born
                                print('Class',(n,k,i1,i2,j1,j2), '0 bars born at ', j1, '\u03B1 +'
                                        , j2, '\u03B1 hat +', n-k-j1-j2, '\u03B2 +', k-i1-i2, '\u03B3 +'
                                        , i1, '\u03B4 +', i2, '\u03B4 hat.') # no bars born
                        else:
                            continue

if __name__ == "__main__":
    if len(sys.argv) > 1:
        pc1321(int(sys.argv[1]))
    else:
        raise SystemExit("usage:  python pc1321.py <n>")
