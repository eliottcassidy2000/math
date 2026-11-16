from check_bridge_reset import count,join
from math import isqrt
def rows(code):
    return tuple(code>>(8*i)&255 for i in range(8))
V0,V1=map(rows,(4578838542729556224,4578836343708259584))
for V,ep in [(V0,(59,72)),(V1,(60,71))]:
    assert count(V,6)[2:]==(158,158,2848)
    assert count(V,6,True)[0]==109 and count(V,6,True)[2:4]==ep
def test(U,r,F):
    h,g,e,o,_=count(U,r)
    x,d=e+o,e-o
    assert d==-1
    vals=[count(join(U,r,V,6),r)[:4] for V in (V0,V1)]
    assert [v[0] for v in vals]==list(F)
    assert vals[0][1:]==vals[1][1:]==(109*g,1424*x-158*g,1424*x-158*g)
def seed(m):
    n=m+2;u=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    u[0]^=1<<(n-1);u[-1]^=1
    return tuple(u)
for m in range(1,5):
    h=1+2**m
    test(seed(m),m+1,(240*h-13,240*h-11))
test((0,),0,(227,229))
assert all(x%d for x in (227,229) for d in range(2,isqrt(x)+1))
print('Matched end-parity readers and direct arbitrary-height samples pass, including 17→4067/4069.')
def strong(A):
    F=(1<<len(A))-1
    for k in range(len(A)):
        s=1<<k
        while True:
            z=s
            for i,r in enumerate(A):
                if s>>i&1:z|=r
            if z==s:break
            s=z
        if s!=F:return False
    return True
U=seed(4)
for c,et,out in [(8871896820243039496,(111,156),7843),(8007222596720985376,(112,155),7845)]:
    V=rows(c);J=join(U,5,V,6)
    assert strong(V) and strong(J)
    assert count(V,6)==(197,count(V,6)[1],318,333,3392)
    assert count(V,6,True)[2:4]==et
    assert count(J,5)[:4]==(out,197,1696*17-333,1696*17-318)
print('Both strong 197-reader joins replayed.')
