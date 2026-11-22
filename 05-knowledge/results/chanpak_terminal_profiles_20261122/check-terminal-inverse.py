from itertools import product
def Z(t):
    d,u,v=t
    return d+2*u,v,d+u
def P(t,m):
    d,u,v=t
    return (1<<m)*(d+2*u+v),0,d+u+v
x,y=(1,0,0),(0,0,1)
dets=[]
for j in range(50):
    a,b,c,d=sum(x),sum(y),x[1],y[1]
    dets.append(a*d-b*c)
    if j==13:
        assert (a,b,c,d)==(47*56,47*47,9*56,9*47)
    x,y=Z(x),Z(y)
assert dets[:14]==[0,2,-3,1,4,-8,5,7,-20,18,9,-47,56,0]
assert all(d3==-d2-d1+d0 for d0,d1,d2,d3 in zip(dets,dets[1:],dets[2:],dets[3:]))
for (D,G),j,m in product(((2,1),(4,1),(26,5),(56,7),(100,13)),range(30),range(1,9)):
    t=(D,0,G)
    for _ in range(j):t=Z(t)
    d,_,g=P(t,m)
    m0=(d//g).bit_length()-1
    assert m0==m and (d>>m0)-g==t[1] and (t[1]==0)==(j==0)
    if j==13:
        K=56*D+47*G
        assert (d+g,g)==((47+56*(1<<m))*K,47*K)
print("unique width, determinant recurrence and rank-one reset PASS")
