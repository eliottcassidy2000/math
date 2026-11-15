from fractions import Fraction
from itertools import product
from math import isqrt
from verify_scaled import h
from verify_discriminant import inverse

LETTERS=((4,7,8),(7,4,8),(7,11,13))
ROOT=(1,3)
def move(s,i):
    x,y=s;a,b,c=LETTERS[i]
    return a*a*(x+y)-x,b*b*(x+y)-y
def address(s):
    w=[]
    while s!=ROOT:
        x,y=s;S=x+y
        if min(x,y)<1:return None
        r=Fraction(x,S)
        cand=[(i,a,c) for i,(a,b,c) in enumerate(LETTERS)
              if Fraction(a*a-1,c*c)<r<Fraction(a*a,c*c)]
        if len(cand)!=1:return None
        i,a,c=cand[0]
        if S%(c*c):return None
        T=S//(c*c)
        xx=a*a*T-x
        if not 0<xx<T:return None
        s=xx,T-xx; w.append(i)
    return tuple(w[::-1])

def need(b):
    if not b:raise RuntimeError("mixed decoder mismatch")

def decode(N):
    if not N&1:return []
    out=[]
    j=2*N%3
    if N==h(*ROOT,j):out.append((ROOT,(),j))
    lo=isqrt(N//700)
    if 700*lo*lo<N:lo+=1
    hi=isqrt(N//35)
    v=4
    while v<=hi:
        S=v
        while S<lo:
            S*=169
        while S<=hi:
            if S>4:
                for ori in (0,1):
                    ans=inverse(S,N,ori)
                    if ans:
                        pair=ans[:2]
                        w=address(pair)
                        if w is not None:out.append((pair,w,ans[2]))
            S*=169
        v*=64
    return out

if __name__=="__main__":
    seen={ROOT}
    for length in range(1,7):
        for w in product(range(3),repeat=length):
            s=ROOT
            for i in w:s=move(s,i)
            need(s not in seen);seen.add(s)
            need(address(s)==w)
            N=h(*s,0)
            need(inverse(sum(s),N,int(s[0]<s[1]))==s+(0,))
            need((s,w,0) in decode(N))
    print(f"{len(seen)} exact mixed-scale host addresses and scalar inversions")
