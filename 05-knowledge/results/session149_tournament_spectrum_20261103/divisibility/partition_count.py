from collections import Counter
from math import gcd, factorial
from fractions import Fraction

def oddparts(n, lo=1):
    if n == 0:
        yield ()
    for x in range(lo,n+1,2):
        for p in oddparts(n-x,x):
            yield (x,)+p

def count(n):
    s=Fraction(0); count=0; mv=10**9
    for part in oddparts(n):
        m=Counter(part)
        z=1
        for i,k in m.items(): z*=i**k*factorial(k)
        E=(sum(a*b*gcd(i,j) for i,a in m.items() for j,b in m.items())-sum(m.values()))//2
        v=(z&-z).bit_length()-1
        mv=min(mv,E-v)
        s += Fraction(1<<E,z); count+=1
    if s.denominator != 1:
        raise ArithmeticError("Nonintegral Burnside sum")
    return s.numerator,count,mv
if __name__ == "__main__":
    for n in [3,4,5,10,20,40,77]:
        a,k,v=count(n)
        actual=(a & -a).bit_length()-1
        print(n,k,v,actual,len(str(a)),str(a)[-30:])
