import math,json

def require(x):
    if not x:raise RuntimeError('check failed')
def hasse(c,k):return sum(math.comb(i,k)*v for i,v in enumerate(c) if i>=k)
def ev(c,t):return sum(v*t**i for i,v in enumerate(c))
def mul(c,d):
    a=[0]*(len(c)+len(d)-1)
    for i,x in enumerate(c):
        for j,y in enumerate(d):a[i+j]+=x*y
    return a
seeds=[[0,1],[0,0,1],[0,0,4,0,1],[0,0,22,0,14,0,1],[0,0,24,0,12,0,1]]
records=[]
for c in seeds:
    H,D,E,F=[hasse(c,k) for k in range(4)];n=len(c)-1
    def state(r):
        q=[0]*r+c
        h,d,e,f=[hasse(q,k) for k in range(4)]
        return (d+9*f)%27,e%9,q
    require(len({state(r)[0] for r in range(27)})==27)
    for par in (0,1):
        rs=[r for r in range(4,58) if (n+r)%2==par and state(r)[0]==0]
        require(len(rs)==1)
        r=rs[0];_,J,q=state(r)
        for t in range(81):
            if t%3==0:w=0
            else:
                s=1 if t%3==1 else -1;u=(s*t-1)//3
                w=s**par*(H+9*u*u*J)
            require(ev(q,t)%81==w%81)
        K=(2*H*E-D*D+H*D)%3
        require(K==2*H*J%3)
        records.append({'Q':c,'parity':par,'r':r,'J':J,'H':H})
z=[-1,0,1];cube=mul(mul(z,z),z);cube[0]+=9;cube[2]-=9
P=[0]*4+cube
require(all(ev(P,t)%81==0 for t in range(81)))
require([hasse(P,k) for k in range(4)]==[0,-18,-81,-136])
print(json.dumps({'normalizations':records,'alias_coefficients':P,'alias_Hasse':[hasse(P,k) for k in range(4)],'all_checks_passed':True},indent=2,sort_keys=True))
