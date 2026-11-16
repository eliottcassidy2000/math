from itertools import product
D=set(product(range(3),range(2),range(2)))-{(0,0,0),(2,1,1)}
C=lambda z:tuple(x-y for x,y in zip((2,1,1),z))
Q=lambda z:(z[0],z[2],z[1])
f3=lambda z:z[1]-z[2]
def f5(z):
    x=f3(z)
    return (x,2*abs(z[0]-1)-1) if x else (0,0)
def f7(z):
    x=f3(z)
    return (x,z[0]-1) if x else (0,0)
for f,n in [(f3,3),(f5,5),(f7,7)]:
    cc={f(z):f(C(z)) for z in D}
    qq={f(z):f(Q(z)) for z in D}
    for z in D:assert cc[f(z)]==f(C(z)) and qq[f(z)]==f(Q(z))
    assert len(cc)==len(qq)==n
    assert sum(k==v for k,v in qq.items())==1
    assert (cc==qq)==(n<7)
C_orbit=lambda z:frozenset([z,C(z)])
co={C_orbit(z) for z in D}
assert len(co)==5
assert sum(s==frozenset(map(Q,s)) for s in co)==3
print('3/5/7 equivariant quotient ladder, plus hostile 5-point complement quotient: PASS')
