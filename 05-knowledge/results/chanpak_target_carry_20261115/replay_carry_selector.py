# Load the published predecessor packet (or its local working name).
from pathlib import Path
import sys
base=Path(__file__).resolve().parent.parent
for directory in ("chanpak_absolute_steering_20261115","chanpak_absolute_steering"):
    if (base/directory).is_dir():
        sys.path.insert(0,str(base/directory))
from fractions import Fraction
from random import Random
from replay_steering import B,Q,T,S,p,output,reverse

def check(x):
    if not x:
        raise AssertionError('selector certificate failed')

def cmp(x,y):
    return (x>y)-(x<y)

def directed(s,L,U,target,side):
    r=s.bit_length()
    if (side==1 and target<0) or (side==-1 and target>=(1<<r)):
        return None
    target=max(0,min((1<<r)-1,target))
    states={(0,0,0,0):(0,0)}
    for i in range(r):
        bL,bU,bs,bt=[int(z) for z in ((L>>i)&1,(U>>i)&1,(s>>i)&1,(target>>(r-1-i))&1)]
        nxt={}
        for (c,d,e,k),(cval,x) in states.items():
            for bx in (0,1):
                for by in (0,1):
                    z=bx+by+c-bs
                    if z not in (0,2):
                        continue
                    cp=z//2
                    key=(cp,cmp(bx,bL) or d,cmp(bx,bU) or e,k or cmp(cp,bt))
                    cand=(2*cval+cp,x+(bx<<i))
                    if key not in nxt or side*(cand[0]-nxt[key][0])>0:
                        nxt[key]=cand
        states=nxt
    eligible=[ans for (c,d,e,k),ans in states.items() if c==0 and d>=0 and e<=0 and k*side<=0]
    return max(eligible,key=lambda a: side*a[0]) if eligible else None

def extract(x,y,r):
    c=0
    val=0
    for i in range(r):
        c=(((x>>i)&1)+((y>>i)&1)+c)//2
        val=2*val+c
    return val

def steer(s,L,U,N):
    r=s.bit_length()
    K=1<<(r+1)
    t0=Fraction(L+U-s,2)
    w=Fraction(U-L,2)
    rho=Fraction(reverse(s,r),1<<r)
    center=B+Q*(Fraction(s*s,4)-t0*t0)+Fraction((T+S)*s,2)+(T-S)*t0+45*s*(1-rho/2)
    E=Q*(2*abs(t0)*w+w*w)+abs(T-S)*w+45*(abs(t0)+w)+90*(r+1)
    target=K*(center-N)/Fraction(135*s,2)
    n,d=target.numerator,target.denominator
    a=directed(s,L,U,n//d,1)
    b=directed(s,L,U,-(-n//d),-1)
    live=[z for z in (a,b) if z is not None]
    c,x=min(live,key=lambda z:(abs(z[0]-target),z[0]))
    model=center-Fraction(135*s*c,2*K)
    check(abs(output(x,s-x)-model)<=E)
    return c,x,target,E

if __name__=='__main__':
    rng=Random(170311)
    k=0
    for s in range(1,512):
        r=s.bit_length()
        for _ in range(5):
            L=rng.randrange(s+1); U=rng.randrange(L,s+1)
            N=output(L,s-L)+rng.randint(-10_000,10_000)
            c,x,t,E=steer(s,L,U,N)
            cs={extract(v,s-v,r) for v in range(L,U+1)}
            check(c==extract(x,s-x,r) and L<=x<=U)
            check(abs(c-t)==min(abs(z-t) for z in cs))
            for v in range(L,U+1):
                z=extract(v,s-v,r)
                check(2*(reverse(v,r)+reverse(s-v,r)-reverse(s,r))==3*z)
            k+=1
    print('PASS: target-routed nearest carries against exhaustive address intervals:',k)
