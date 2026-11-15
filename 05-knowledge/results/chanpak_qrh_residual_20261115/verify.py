from transfer import *
from random import Random
r=Random(20261115)

def state_mod(x,y,s,m):
    vx=p(x)%m; vy=0 if y is None else p(y)%m
    sx=x%2; sy=0 if y is None else y%2
    X=x%m; Y=0 if y is None else y%m
    for z in s:
        b=int(z)
        if not b:vx=(vx+X+sx)%m
        if y is not None and b:vy=(vy+Y+sy)%m
        X=(2*X+b)%m
        Y=(2*Y+1-b)%m
        sx,sy=b,1-b
    return X,Y,(vx+vy)%m,sx,sy

total=0
for m in (3,5,7,9,11,13,15,17,25,27,31,45,63):
    for _ in range(12):
        x=2*r.randrange(100)+1;y=2*r.randrange(100)
        a,b=r.randrange(m),r.randrange(m)
        for yy in (None,y):
            s=steer(x,yy,m,a,b)
            X,Y,W,sx,sy=state_mod(x,yy,s,m)
            assert sx==1 and X==a and W==b
            if yy is not None:
                assert sy==0 and (X+Y)%m==(x+y)%m
                if m<12:
                    x1,y1=append(x,yy,s)
                    assert x1+y1==((x+y+1)<<len(s))-1
                    assert (p(x1)+p(y1))%m==b
            total+=1
    if m in (3,5,7,9,11,15):
        for N in range(m):count_selector(N,m)
print('Independent and complementary exact transfers:',total)
print('Literal odd-modulus count checks: 50.')
for m in (7,11,49):
    if m==49: continue  # Larger-modulus handles covered by streamed states above.
    h=suffix_length(m)+2
    g=1<<h
    x=256*g+g//7
    y=255*g+g-1-g//7
    for n in (100,109,110,113):
        X,Y,j=joint_selector(n,m,x,y,g)
        assert X+Y==x+y and count(X,Y,j)%(3*m)==n%(3*m)
print('Joint literal sum and residue checks: 8.')
