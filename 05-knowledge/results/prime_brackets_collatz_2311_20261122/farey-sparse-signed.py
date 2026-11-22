"""No signed shortened Collatz script with <=2 odd steps on the E61 Farey charts."""
from itertools import combinations,product
from math import isqrt
N=(1<<60)*((1<<61)-1)-1
C4=2465,2293760,1344270,712944,1740,986,423385
C8=37265,144834560,84839910,57602544,24360,14906,33838015

def f(x):
    s=0
    while x:
        if not x&1:s+=(x>>1)+((x>>1)&1)
        x>>=1
    return s

def candidate(chart,target,A,B,C,D):
    h,q,c,t,e,r,k=chart
    # x=A*z+B, y=C*z+D, with both positive.
    z0=max(0,-((B-1)//A),-((D-1)//C))
    aa=q*A*C;bb=q*(A*D+B*C)+t*A+c*C
    cc=k+t*B+c*D+q*B*D-target
    if (aa*z0+bb)*z0+cc>0:return (0,None)
    z=(isqrt(bb*bb-4*aa*cc)-bb)//(2*aa)
    x=A*z+B;y=C*z+D
    R=-((aa*z+bb)*z+cc)-r*(x&1)-e*(y&1)
    if not 0<=R<=h*(x+y):return (0,None)
    if R%(2*h):return (1,None)
    if R!=h*(f(x)+f(y)):return (2,None)
    return (2,(x,y))

def scan(chart,target):
    h,q,c,t,e,r,k=chart
    bb=pp=counter=top=0
    hits=set()
    for sw in (False,True):
        cc,tt=(t,c) if sw else (c,t)
        xhi=(target-k-cc)//(q+tt)
        if xhi<=0:continue
        for m in range(3):
            dmax=xhi.bit_length()-1+2*m
            top=max(top,dmax)
            three=3**m
            for d in range(max(1,m),dmax+1):
                B=1<<d
                inv=pow(three,-1,B)
                for pos in combinations(range(d),m):
                    weights=[3**(m-1-j)<<t0 for j,t0 in enumerate(pos)]
                    for s in product((-1,1),repeat=m):
                        counter+=1
                        carry=sum(a*b for a,b in zip(s,weights))
                        x0=-carry*inv%B; y0=(three*x0+carry)//B
                        assert y0>=0 and (three*x0+carry)%B==0
                        C=(three,y0,B,x0) if sw else (B,x0,three,y0)
                        a,z=candidate(chart,target,*C)
                        bb+=a>=1;pp+=a>=2
                        if z:hits.add(z)
    return counter,top,bb,pp,sorted(hits)

if __name__=="__main__":
    for n,c in ((4,C4),(8,C8)):
        h,q,cy,tx,e,r,k=c
        x,y=29,8 # + at 29, then /2 twice, - at 11, then /2: 29 ->44->22->11->16->8.
        def height(x,y):return k+q*x*y+tx*x+cy*y+r*(x&1)+e*(y&1)+h*(f(x)+f(y))
        assert (x,y) in scan(c,height(x,y))[-1]
        assert (y,x) in scan(c,height(y,x))[-1]
        print(n,scan(c,N),flush=True)
