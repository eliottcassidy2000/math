from math import isqrt
g=sum(1<<j for j in (0,2,4,5,6,10,11))
C={0}
for j in range(12):
    b=(g<<j)|(1<<23)
    C|={w^b for w in C}
trio=[(0,2,4,5,6,10,11,23),(1,3,7,8,12,13,15,18),(9,14,16,17,19,20,21,22)]
m=[sum(1<<i for i in p) for p in trio]
assert all(u in C and u.bit_count()==8 for u in m)
order=sum((list(p) for p in trio),[])
assert sorted(order)==list(range(24))
raw=4*2**50+4*2**25-1
values={}
for word in C:
    x=[(word>>i)&1 for i in order]
    D=1
    Q=2**25
    for i in range(24):
        for j in range(i+1,24):
            if x[i] and not x[j]:
                d=j-i-1
                D+=2**d+(d==0)
                Q+=2**(23-d)*(3**d+1)
    values[word]=4*Q*2**25+3*Q+D*2**25-D
assert {w for w,h in values.items() if h==raw}=={0,m[2],m[1]|m[2],sum(m)}
assert {h for h in values.values() if isqrt(h+2)**2==h+2}=={raw}
assert len(set(values.values()))==4067
print("single-apex trio: four square-preserving codewords and 4067 exact heights")
