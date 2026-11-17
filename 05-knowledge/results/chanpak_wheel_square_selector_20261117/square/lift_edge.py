from math import isqrt
from itertools import pairwise
base=[1,8,28,21,15,34,2,7,18,31,33,16,9,27,22,14,35,29,20,5,11,25,24,12,4,32,17,19,30,6,10,26,23,13,3]
# Consecutive block-end/block-start values from the audited Gerbicz table.
J=[(1,24),(74,26),(76,68),(18,31),(81,19),(69,75),(25,11),(11,70),(20,5),(5,4),(4,12),(12,13),(63,37),(87,82),(32,17),(67,77),(27,22),(72,9),(9,7),(7,29),(79,21),(71,10),(10,6),(6,30),(80,64),(14,2),(2,23),(73,8),(8,28),(78,66),(16,84),(34,15),(65,35),(85,36),(86,83),(33,3)]
def ck(x):
    if not x:raise ArithmeticError('control failed')
def square(x):return isqrt(x)**2==x
def check(a):ck(a[0]==1 and a[-1]==3 and len(a)%2 and set(a)==set(range(1,len(a)+1)) and all(square(x+y) for x,y in pairwise(a)))
def lift(a):
    ck(len(a)%2 and a[0]==1 and a[-1]==3)
    result=[1];params=[];single=[1,3]
    for (_,u),(v,_) in pairwise(J):
        if u==v:result.append(u);single.append(u)
        else:
            rev=u>v;c=min(u,v)-25
            ck(abs(u-v)==50 and -12<=c<=12)
            params.append(c)
            b=[25*x+(1-2*(i%2))*c for i,x in enumerate(a)]
            if rev:b=b[::-1]
            ck((b[0],b[-1])==(u,v));result.extend(b)
    result.append(3)
    ck(sorted(params)==list(range(-12,13)) and sorted(single)==list(range(1,13)))
    ck(all(square(x+y) for x,y in J))
    return result
a=base
for m in range(4):
    check(a);ck(len(a)==(71*25**m-1)//2)
    ck(any({x,y}=={19*25**m,30*25**m} for x,y in pairwise(a)))
    a=lift(a) if m<3 else a
print('PASS protected edge at 35, 887, 22187, 554687')
