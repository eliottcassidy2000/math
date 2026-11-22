"""Independent integer-only exhaustive inverse for the canonical spine-8 chart.

CLI limit is on the smaller normalized coordinate; omit limits for the entire
exact search, which is much too large for an interactive run.
"""
from argparse import ArgumentParser
from fractions import Fraction
from math import gcd,isqrt

E=(1<<60)*((1<<61)-1)
H,q,c,t,e,r,k=(37265,144834560,84839910,57602544,24360,14906,33838015)
P1,P2=29,257
SQ={p:{(s*s)%p:tuple({s,(-s)%p}) for s in range(p)} for p in (P1,P2)}

def f(y):
    z=0
    while y:
        y,j=divmod(y,2)
        if not j:z+=y+(y&1)
    return z

def height(x,y):
    return k+q*x*y+t*x+c*y+r*(x&1)+e*(y&1)+H*(f(x)+f(y))

def ceil(a,b):
    return -((-a)//b)

def crt(a,m,b,n):
    d=gcd(m,n)
    if (b-a)%d:return None
    p=n//d
    return a+m*(((b-a)//d)*pow(m//d,-1,p)%p if p!=1 else 0),m*p

def roots(a,b,c,p):
    a,b,c=a%p,b%p,c%p
    if not a:
        if b:return (-(c*pow(b,-1,p))%p,)
        return () if c else None  # None means unrestricted.
    return tuple(((-b+d)*pow(2*a,-1,p))%p for d in SQ[p].get((b*b-4*a*c)%p,()))

def circ(n):
    if n<=0:return 0
    x=1<<((n.bit_length()+2)//3)
    while True:
        y=(2*x+n//(x*x))//3
        if y>=x:return x
        x=y

def bounds(R,Q,T,C,X0,W,a,b,y0):
    """Min and max of a*d+b*((R-T*x)/(Q*x+C)-y0) on integer 0<=d<W."""
    J=Q*R+T*C
    def val(d):
        x=X0+d;D=Q*x+C
        return Fraction((a*d-b*y0)*D+b*(R-T*x),D)
    ends=[val(0),val(W-1)]
    if W>1 and J:
        # Consecutive difference is a-b*J/(D*(D+Q)). It is
        # monotone, increasing for J>0 and decreasing for J<0.
        l,h=0,W-1
        while l<h:
            i=(l+h)//2
            D=Q*(X0+i)+C
            dif=a*D*(D+Q)-b*J
            side=(dif>=0 if J>0 else dif<=0)
            if side:h=i
            else:l=i+1
        ends.append(val(l))
    return min(ends),max(ends)

class Branch:
    def __init__(self,N,i,j,sw):
        self.N,self.sw=N,sw
        self.q=q*20
        self.o,self.p=(j,i) if sw else (i,j)
        self.s,self.u=(2,10) if sw else (10,2)
        self.x,self.y=(c,t) if sw else (t,c)
        self.T=self.s*(self.x+q*self.p)
        self.C=self.u*(self.y+q*self.o)
        self.R=N-height_affine(i,j)
        self.axis_bit=(i if sw else j)&1
        self.axis_parity=r if sw else e
        self.anchor_parity=e if sw else r
        self.height_candidates=self.blocks=0
        self.cut=max(1,circ(max(0,self.R)//self.q))
    def unscale(self,X,Y):
        x=self.s*X+self.o;y=self.u*Y+self.p
        return (y,x) if self.sw else (x,y)
    def test(self,X,Y):
        self.height_candidates+=1
        x,y=self.unscale(X,Y)
        if height(x,y)==self.N:return x,y
    def row(self,X):
        """Group all low suffixes with the necessary slope congruence."""
        x=self.s*X+self.o
        slope=self.y+q*x
        rem=self.N-k-self.x*x-self.anchor_parity*(x&1)-H*f(x)-self.axis_parity*self.axis_bit
        low=max(self.u*(X+int(self.sw))+self.p,ceil(rem,slope+H))
        high=rem//slope
        if low>high:return
        # A short exact progression is cheaper than constructing low suffixes.
        if high-low<16:
            for y in range(low+(self.p-low)%self.u,high+1,self.u):
                hit=self.test(X,(y-self.p)//self.u)
                if hit:yield hit
            return
        scale=max(2,isqrt(rem*H)//slope)
        B=1<<(scale.bit_length()-1)
        Nbits=B.bit_length()-1
        g=gcd(slope*B,H);mod=H//g
        inv=pow(slope*B//g,-1,mod) if mod>1 else 0
        for lowbit in range(self.axis_bit,B,2):
            rr=rem-slope*lowbit-H*f(lowbit)
            if rr<0:break
            if rr%g:continue
            reverse=int(f"{lowbit:0{Nbits}b}"[::-1],2)
            D=slope*B+H*(B-1-reverse)
            lower=max(0,ceil(rr-H,D+H),ceil(low-lowbit,B))
            upper=min(rr//D,(high-lowbit)//B)
            res=(rr//g)*inv%mod
            d=gcd(B,self.u)
            if (self.p-lowbit)%d:continue
            sub=self.u//d
            bound=((self.p-lowbit)//d)*pow(B//d,-1,sub)%sub if sub>1 else 0
            check=crt(res,mod,bound,sub)
            if check is None:continue
            res,step=check
            for z in range(lower+(res-lower)%step,upper+1,step):
                y=B*z+lowbit
                if y%self.u!=self.p:continue
                hit=self.test(X,(y-self.p)//self.u)
                if hit:yield hit
    def slab(self,X0,W):
        Q,T,C,R=self.q,self.T,self.C,self.R
        der=Fraction(Q*R+T*C,(Q*X0+C)**2).limit_denominator(W)
        a,b=der.numerator,der.denominator
        y0=(R-T*X0)//(Q*X0+C)
        lo,_=bounds(R-H*(self.o+self.p),Q,T+H*self.s,C+H*self.u,X0,W,a,b,y0)
        _,hi=bounds(R,Q,T,C,X0,W,a,b,y0)
        self.blocks+=1
        inv=pow(a,-1,b) if b>1 else 0
        for z in range(ceil(lo.numerator,lo.denominator),hi.numerator//hi.denominator+1):
            D=Q*X0+C
            J=b*y0+z
            rr1=roots(-Q*a,b*T+Q*J-a*D,b*T*X0+D*J-b*R,P1)
            rr2=roots(-Q*a,b*T+Q*J-a*D,b*T*X0+D*J-b*R,P2)
            if rr1==() or rr2==():continue
            base=(inv*z%b,b)
            for u in rr1 if rr1 is not None else (None,):
                branch1=base if u is None else crt(*base,u,P1)
                if branch1 is None:continue
                for v in rr2 if rr2 is not None else (None,):
                    branch2=branch1 if v is None else crt(*branch1,v,P2)
                    if branch2 is None:continue
                    d,step=branch2
                    for j in range(d,W,step):
                        X=X0+j;Y=y0+(z-a*j)//b
                        if Y>=X+int(self.sw):
                            hit=self.test(X,Y)
                            if hit:yield hit
    def scan(self,L,U):
        if self.R<0:return
        U=min(U,isqrt(self.R//self.q))
        X=L
        while X<=U:
            W=min(U-X+1,max(1,X//self.cut))
            if X<self.cut*2:
                yield from self.row(X)
            else:
                yield from self.slab(X,W)
            X+=W

def height_affine(x,y):
    return k+q*x*y+t*x+c*y+r*(x&1)+e*(y&1)

def scan(N,L=0,U=1<<100):
    for i in range(10):
        if (N-k-t*i-r*(i&1))%5:continue
        for j in range(2):
            for sw in (False,True):
                branch=Branch(N,i,j,sw)
                for hit in branch.scan(L,U):yield hit

if __name__=='__main__':
    cli=ArgumentParser(description=__doc__)
    cli.add_argument('--target',type=int,default=E-1)
    cli.add_argument('--low',type=int,default=0)
    cli.add_argument('--high',type=int,default=1<<100)
    cli.add_argument('--check',action='store_true')
    cli.add_argument('--budget',action='store_true')
    a=cli.parse_args()
    if a.budget:
        for i in range(10):
            if (a.target-k-t*i-r*(i&1))%5:continue
            for j in range(2):
                for sw in (False,True):
                    B=Branch(a.target,i,j,sw)
                    if B.R<0:continue
                    u=isqrt(B.R//B.q)
                    X=B.cut*2
                    slabs=0
                    while X<=u:
                        w=X//B.cut
                        n=ceil(min(u+1,(w+1)*B.cut)-X,w)
                        slabs+=n
                        X+=n*w
                    print(f"left={i} right={j} swapped={int(sw)} xmax={u} "
                          f"rows={min(u+1,B.cut*2)} slabs={slabs}")
    elif a.check:
        for x,y in ((0,0),(10,2100),(7,1184),(4442,0),(1234,620),(0,45678),
                    (76543,391),(771,156784)):
            N=height(x,y)
            u=min(x//10,y//2)
            assert (x,y) in set(scan(N,u,u))
        # Force the exact tangent stage on a planted large point.
        x,y=140_000_000_000_000,280_000_000_000_000
        N=height(x,y)
        br=Branch(N,x%10,y%2,False)
        assert (x,y) in set(br.slab(x//10-10,100))
        for x,y in ((0,0),(0,5),(7,0),(1,1),(7,5),(11,2)):
            N=height(x,y)
            brute=set()
            for xx in range((N-k)//t+1):
                lim=(N-k-t*xx)//(q*xx+c)
                for yy in range(lim+1):
                    if height(xx,yy)==N:brute.add((xx,yy))
            found=list(scan(N))
            assert set(found)==brute and len(found)==len(brute)
        for z in range(1,10000):
            a0=circ(z)
            assert a0**3<=z<(a0+1)**3
        print('exact axis, progression, tangent and full inverse controls PASS')
    else:
        found=list(scan(a.target,a.low,a.high))
        print(len(found),'solutions:', found)
