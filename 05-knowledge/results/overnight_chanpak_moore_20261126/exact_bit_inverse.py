"""All-spoke F1 inverse on even b with bounded odd part; exact integer arithmetic."""
from math import isqrt
import argparse

def F(u):
    v=1
    while u:
        q,bit=divmod(u,2)
        if not bit:v+=q+(q&1)
        u=q
    return v

def H(b,c):
    return b*c+F(b-1)*F(c-1)-(b&1)*(c&1)

def rev(u,i):
    v=0
    for _ in range(i):v=(v<<1)|(u&1);u>>=1
    return v

def hits(E, dmax):
    eligible=rows=checks=0
    answer=[]
    for d in range(1,dmax+1,2):
        h=F(d-1)
        w0=isqrt(E*h//2)
        lo,hi=max(1,w0-9*h),w0+9*h+1
        k=(hi//d).bit_length()-1
        while k>=1:
            u=d<<k
            if 3*u<2*lo:break
            k-=1
            signs=[sgn for sgn in (1,-1) if (u-2*sgn)%h==0]
            if not signs:continue
            eligible+=1
            z0=max(0,-((u-lo)//h)); z1=(hi-u)//h
            if z1<z0:continue
            for s in range(1,k+2):
                i=k+1-s;ell=1<<i;b=d<<s
                if z0>=ell:break
                rows+=1
                for Z in range(z0,min(z1,ell-1)+1):
                    r=(ell-1)^rev(Z,i)
                    eps=(Z&1) if i else 0
                    w=u+h*Z
                    for sgn in signs:
                        beta=(w-2*sgn)//h+eps
                        assert beta%2==0
                        for q in (beta,beta-1):
                            if q<0:continue
                            checks+=1
                            M=2*w*q+b*(r+1)+h*(F(q)+beta+F(r)-1)+sgn
                            if M!=E:continue
                            c=2*ell*q+r+1
                            assert H(b,c)==E-sgn and H(b,c+ell)==E+sgn
                            answer.append((b,c,i,sgn))
    return answer,(eligible,rows,checks)

def controls():
    # Exhaust a small right-address window and compare to the exact inverse
    # on a collection of its first midpoint encounters.
    observed={}
    for b in range(2,38,2):
        for v in range(64):
            for i in range(6):
                ell=1<<i
                if v&ell:continue
                a,c=H(b,v+1),H(b,v+ell+1)
                if abs(a-c)==2:
                    observed.setdefault((a+c)//2,[]).append((b,v+1,i,(c-a)//2))
    for E,W in observed.items():
        A,_=hits(E,37)
        assert set(W).issubset(A)
    print("small controls:",len(observed),"midpoints",sum(map(len,observed.values())),"literal arcs")

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument("exponent",type=int)
    p.add_argument("odd_part_max",type=int)
    args=p.parse_args()
    controls()
    E=(1<<(args.exponent-1))*((1<<args.exponent)-1)
    A,stats=hits(E,args.odd_part_max)
    print("exponent",args.exponent,"odd-part maximum",args.odd_part_max,
          "eligible leading lanes, dyadic splits, exact reads",stats)
    print("witnesses",A)
