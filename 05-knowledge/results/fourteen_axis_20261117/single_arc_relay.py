from replay import adj,menu,strong,extend,prime
from math import gcd
A=(8,17,35,54,101,195,143,31)
B=list(A);B[4]^=1<<7;B[7]^=1<<4
V=[menu(X,3) for X in (A,B)]
if V!=[(315,198,235,37),(329,198,235,37)]: raise AssertionError(V)
for word,expected in (((1,0,0,0,1),[(7343,6595,6632,37),(7357,6595,6632,37)]),
                      ((0,1,0,1),[(4571,3353,3390,37),(4585,3353,3390,37)])):
    Y=[A,B];Z=V
    for ch in word:
        Y=[adj(X,3,ch) for X in Y];Z=[extend(T,ch) for T in Z]
    if [menu(T,3) for T in Y]!=expected or Z!=expected or not all(map(strong,Y)):
        raise AssertionError(Z)
    if any((Y[0][i]^Y[1][i])!=mask for i,mask in enumerate([0,0,0,0,1<<7,0,0,1<<4]+[0]*len(word))):
        raise AssertionError
    print(word,expected)
if gcd(13190,37)!=1 or not (prime(1049) and prime(1051)):
    raise AssertionError
if pow(2,42,147)!=1:raise AssertionError
