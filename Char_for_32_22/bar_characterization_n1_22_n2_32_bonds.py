import sys
import scipy as scipy
import math

# Counts number of bars for any number of copies of the bonds 1-3-2-1 and
# 1-2-2-1.
# Input: n1 and n2 -- number of bonds of types 22 and 32 respectively
# Output: Prints out the class, length of bars born, number of bars born, and
# the class energy value.

# Define an abriviation for the factorial function.
def fac(z):
    return math.factorial(z)
# Defines a function that takes in a class and prints the number of critical
# points of a given type and the classes' corresponding energy value.
def eng(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22):
     return (j11, '\u03B1_1 +', n1-k1-j11, '\u03B2_1 +', k1-i11, '\u03B3_1 +', i11, '\u03B4_1 +'
            , j21, '\u03B1_2 +', j22, '\u03B1_2 hat +', n2-k2-j21-j22, '\u03B2_2 +'
            , k2-i21-i22, '\u03B3_2 +', i21, '\u03B4_2 +', i22, '\u03B4_2 hat.')

# Iterates through all classes of critical points given by the following
# conditions:
# k1<=n1, i11<=k1, j11<=n1-k1
# k2<=n2, i21+i22<=k2, j21+j22<=n2-k2
# k1+k2<=n1+n2, i11+i21+i22<=k1+k2, j11+j21+j22<=(n1+n2)-(k1+k2).

# Then it uses the formulas given in the paper to calculate the number of bars
# in each class.
# Finally, it prints the class, length of bars born, corresponding energy value,
# and number of bars born.

def pc_1221_1321(n1,n2):
    for k1 in range(0,n1+1):
        for k2 in range(0,n2+1):
            for i11 in range(0,k1+1):
                for i21 in range(0,k2+1):
                    for i22 in range(0,k2+1):
                        for j11 in range(0,n1-k1+1):
                            for j21 in range(0,n2-k2+1):
                                for j22 in range(0,n2-k2+1):
                                    if i11+i21+i22<=k1+k2 and i11<=k1 and i21+i22<=k2 and j11+j21+j22<=n1+n2-k1-k2 and j21+j22<=n2-k2 and j11<=n1-k1:
                                        avals=[]
                                        bvals=[]
                                        dvals=[]
# Using other restrictions, the classes get sorted into five groups by bar length:
# gamma_1-beta_1, gamma_2-alpha_2, delta_2-beta_2, semi infinite, and no bars born
                                        if n1-k1-j11>0 and n2-k2-j21-j22>=0:
                                            for l in range(0,k1-i11+1):
                                                d=(2**(n1-i11-j11))*((-1)**(l)*fac(n1)*fac(n2)
                                                  /(fac(n1-k1-j11+l)*fac(k1-i11-l)*fac(i11)*fac(j11)
                                                  *fac(j21)*fac(j21)*fac(n2-k2-j21-j22)
                                                  *fac(i21)*fac(i22)*fac(k2-i21-i22)))
                                                dvals.append(d)
                                            print('Class',(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22)
                                                    , u'\u03B3_1-\u03B2_1 bars born: ', (sum(dvals))
                                                    , ' at ', eng(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22))
                                                    #gamma_1-beta_1 length bars
                                        elif j22 > 0 and n1-(i11+j11)==0:
                                            for l in range(0,k2-i21-i22+1):
                                                a=((-1)**(l)*fac(n1)*fac(n2)
                                                  /(fac(i11)*fac(j11)*fac(n1-k1-j11)*fac(k1-i11)
                                                  *fac(n2-k2-j21-j22)*fac(k2-i21-i22-l)*fac(i21)
                                                  *fac(i22)*fac(j21)*fac(j22+l)))
                                                avals.append(a)
                                            print('Class',(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22)
                                                    , u'\u03B3_2-\u03B1_2 bars born: ', (sum(avals))
                                                    , ' at ', eng(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22))
                                                    #gamma_2-alpha_2 length bars
                                        elif i11+i22==k1+k2 and j22==0 and i21==0 and j11+j21<n1+n2-k1-k2:
                                            for l in range(0,i21+1):
                                                b=((-1)**(l)*fac(n1)*fac(n2)
                                                  /(fac(i11)*fac(j11)*fac(n1-k1-j11)*fac(k1-i11)
                                                  *fac(n2-k2-j21-j22+l)*fac(k2-i21-i22)*fac(i21-l)
                                                  *fac(i22)*fac(j21)*fac(j22)))
                                                bvals.append(b)
                                            print('Class',(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22)
                                                    , u'\u03B4_2-\u03B2_2 bars born: ', (sum(bvals))
                                                    , ' at ', eng(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22))
                                                    #delta_2-beta_2 length bars
                                        elif i11+i22==k1+k2 and i21==0 and j22==0 and j11+j21==n1+n2-k1-k2:
                                            c=2**(n1-i11-j11)*(fac(n1)*fac(n2)
                                              /(fac(i11)*fac(j11)*fac(n1-k1-j11)*fac(k1-i11)
                                              *fac(n2-k2-j21-j22)*fac(k2-i21-i22)*fac(i21)
                                              *fac(i22)*fac(j21)*fac(j22)))
                                            print('Class',(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22)
                                                    , 'semi-\u221E bars born:  ', c
                                                    , ' at ', eng(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22))
                                                    #semi-infinite bars
                                        else:
                                            print('Class',(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22)
                                                    , '0 bars born at              ', eng(n1,n2,k1,k2,i11,j11,i21,i22,j21,j22))
                                    else:
                                        continue

if __name__ == "__main__":
    if len(sys.argv) > 1:
        pc_1221_1321(int(sys.argv[1]), int(sys.argv[2]))
    else:
        raise SystemExit("usage:  python pc_1221_1321(n1,n2) <n1,n2>")
