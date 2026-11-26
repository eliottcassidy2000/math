"""Exact all-spoke inverse for one rooted right bank and a binary left block."""
from exact_bit_inverse import F, rev
from math import isqrt
from pathlib import Path
from sys import argv

def count(v, h, e, R):
    return h*F(v)-e*(v%2==0)+R*(v+1)

def inverse(E,h,e,R):
    assert 0<=e<=h<=R
    w0=isqrt(E*h//2)
    lo,hi=max(1,w0-12*h),w0+12*h+1
    i=(hi//R).bit_length()-1
    pairs=[]; checks=0
    while i>=0:
        ell=1<<i
        u=R*ell
        if u+h*(ell-1)<lo:break
        z0=max(0,-((u-lo)//h));z1=min(ell-1,(hi-u)//h)
        gamma=e if not i else 0
        signs=[s for s in (1,-1) if (u+gamma-2*s)%h==0]
        for Z in range(z0,z1+1):
            if not signs:break
            r=(ell-1)^rev(Z,i); eps=(Z&1) if i else 0
            w=u+h*Z
            for s in signs:
                beta=(w+gamma-2*s)//h+eps
                if beta%2:continue
                for q in (beta,beta-1):
                    if q<0:continue
                    checks+=1
                    v=2*ell*q+r
                    M=2*w*q+R*(r+1)+h*(F(q)+beta+F(r)-1)-e*(v%2==0)+s
                    if M==E:
                        assert count(v,h,e,R)==E-s and count(v+ell,h,e,R)==E+s
                        pairs.append((v,i,s))
        i-=1
    return pairs,checks
if __name__=='__main__':
    bank=[tuple(map(int,l.split())) for l in Path(argv[1]).read_text().splitlines() if l.strip()]
    for h,e,R in bank[:70]:
        for v in range(48):
            for i in range(6):
                ell=1<<i
                if v&ell:continue
                a,b=count(v,h,e,R),count(v+ell,h,e,R)
                if abs(a-b)==2:
                    assert (v,i,(b-a)//2) in inverse((a+b)//2,h,e,R)[0]
    print("small rooted-bank controls passed")
    for p in map(int,argv[2:]):
        E=(1<<(p-1))*((1<<p)-1)
        out=[inverse(E,*row) for row in bank]
        print("exponent",p,"bank size",len(bank),"exact reads",sum(t[1] for t in out),
              "witnesses",[(B,t[0]) for B,t in zip(bank,out) if t[0]])
