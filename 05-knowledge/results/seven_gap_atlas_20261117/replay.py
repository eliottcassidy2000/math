"""Independent certificates for dyadic gap witnesses and the modular spine."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import cache
from math import gcd, isqrt
from pathlib import Path
import json


def need(ok, reason="certificate failed"):
    if not ok:
        raise AssertionError(reason)


def prime(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))


def arcs(rows):
    n=len(rows)
    need(all(0<=r<1<<n and not(r>>i&1) for i,r in enumerate(rows)))
    need(all((rows[i]>>j&1)+(rows[j]>>i&1)==1
             for i in range(n) for j in range(i)))


def subset(rows):
    arcs(rows); n=len(rows); f=(1<<n)-1
    d=[[0]*n for _ in range(f+1)]
    h=[1]+[0]*f
    for s in range(1,f+1):
        for i,r in enumerate(rows):
            if not s>>i&1:continue
            z=s^(1<<i)
            d[s][i]=sum(d[z][j] for j in range(n) if r>>j&1) if z else 1
        h[s]=sum(d[s])
    return h,d


def profile(rows,r):
    h,d=subset(rows); f=len(h)-1
    E=sum(v[r]*h[f^s] for s,v in enumerate(d) if not s.bit_count()&1)
    O=sum(v[r]*h[f^s] for s,v in enumerate(d) if s.bit_count()&1)
    return h[-1],h[f^(1<<r)],E,O


def build_root(rows,r,z):
    T=list(rows)
    for j in range(z.bit_length()):
        bit=z>>j&1;n=len(T)
        T=[a+((i!=r or not bit)<<n) for i,a in enumerate(T)]
        T.append(bit<<r)
    return T


def transfer(H,g,E,O,z):
    for j in range(z.bit_length()):
        if z>>j&1:
            H,E,O=H+2*E,E+O-g,E+O
        else:
            E,O=E+O,E+O
    return H


def strong_core(rows):
    """For these nonunit counts exactly one SCC has more than one vertex."""
    f=(1<<len(rows))-1
    reach=[]
    for i in range(len(rows)):
        seen=1<<i
        for _ in rows:
            z=seen
            for j,r in enumerate(rows):
                if seen>>j&1:z|=r
            if z==seen:break
            seen=z
        reach.append(seen)
    groups={tuple(j for j in range(len(rows)) if v>>j&1 and reach[j]>>i&1)
            for i,v in enumerate(reach)}
    big=[g for g in groups if len(g)>1]
    need(len(big)==1)
    I=big[0]
    return tuple(sum((rows[i]>>j&1)<<z for z,j in enumerate(I)) for i in I)


def modules(rows):
    f=(1<<len(rows))-1; out=[]
    for s in range(1,f):
        if s.bit_count()<2:continue
        if all(r&s in (0,s) for i,r in enumerate(rows) if not s>>i&1):
            out.append(s)
    return out


def pack_poly(rows):
    """Directed odd cycles anchored at their least vertex, then disjoint packs."""
    n=len(rows);f=(1<<n)-1
    bins=defaultdict(list)
    for u in range(n):
        F=f^((1<<(u+1))-1); B=1<<u
        @cache
        def end(s,v):
            if not s:return int(v==u)
            if not s>>v&1:return 0
            z=s^(1<<v)
            if not z:return rows[u]>>v&1
            return sum(end(z,w) for w in range(u+1,n) if z>>w&1 and rows[w]>>v&1)
        k=F
        while k:
            if k.bit_count()>=2 and k.bit_count()%2==0:
                c=sum(end(k,w) for w in range(u+1,n) if k>>w&1 and rows[w]>>u&1)
                if c:bins[B].append((k|B,c))
            k=(k-1)&F
    @cache
    def P(s):
        if not s:return {0:1}
        B=s&-s
        L=Counter({i+1:c for i,c in P(s^B).items()})
        for j,c in bins[B]:
            if s&j==j:
                for i,q in P(s^j).items():L[i+1]+=2*c*q
        return dict(L)
    return P(f)


def one(w):
    n=len(w)
    return 1+sum((1<<(j-i-1))+(j==i+1) for i in range(n) if w[i]=='1'
                 for j in range(i+1,n) if w[j]=='0')


def apex(w):
    n=len(w)
    return tuple([sum((b=='1')<<(i+1) for i,b in enumerate(w))]
                 +[sum(1<<(j+1) for j in range(i+1,n))+(m=='0') for i,m in enumerate(w)])


def two(rows):
    h,_=subset(rows)
    return h[-1],sum(v*h[(len(h)-1)^i] for i,v in enumerate(h))


def ordinal(A,B):
    n=len(A); m=len(B)
    return tuple(a|(((1<<m)-1)<<n) for a in A)+tuple(b<<n for b in B)


def c3(A,B):
    # A -> B -> c -> A
    n=len(A);m=len(B)
    return tuple(a|(((1<<m)-1)<<n) for a in A)+tuple(
        (b<<n)|(1<<(m+n)) for b in B)+((1<<n)-1,)


def block_pair(M):
    a=(M&-M).bit_length()-1; m=M>>a
    d=max(a,1)
    while pow(2,d,m)!=1%m:d+=1
    L=2*a+d+3
    def word(j):
        return ''.join('1' if i in (a,j) else '0' for i in range(L))
    return word(a+1),word(a+1+d)


def make_word(s,M,r,k=7):
    # This anchor puts H/2^s in (18,28) for every sufficiently large s.
    prefix='10011'+'0'*k
    suffix='0'*(5+k)
    A,B=block_pair(M)
    g=gcd(2,M); t=M//g-1
    n=s+5; z=n-len(prefix+suffix)-t*len(A)
    need(z>=0)
    x=prefix+'0'*z
    h=one(x+A*t+suffix)
    delta=(one(B)-one(A))%M
    need(gcd(delta,M)==g and (r-h)%g==0)
    y=((r-h)//g)*pow(delta//g,-1,M//g)%(M//g) if M//g>1 else 0
    out=x+B*y+A*(t-y)+suffix
    count=one(out)
    need(count%M==r%M and (18<<s)<count<(28<<s))
    return out,count


def replay():
    f=Path(__file__).with_name('witnesses.json')
    atlas={int(k):v for k,v in json.loads(f.read_text()).items()}
    expected={7*p for s in range(2,7) for p in range(((18<<s)+6)//7,(28<<s)//7+1)
              if prime(p)}
    need(set(atlas)==expected)
    W=defaultdict(list)
    for n in range(2,(max(atlas)-1).bit_length()+2):
        for q in range(1<<(n-2)):
            w='1'+(f'{q:0{n-2}b}' if n>2 else '')+'0'
            H=one(w)
            if H in atlas:W[H].append(w)
    for N,item in atlas.items():
        P=item['prime'];s=item['scale']; need(prime(P) and N==7*P and 18<<s<=N<=28<<s)
        rows,r,z=item['root_construction']
        H,g,E,O=profile(rows,r)
        B=build_root(rows,r,z); poly=pack_poly(B)
        core=strong_core(B);hc,_=subset(core); m=modules(core)
        need(hc[-1]==N and all(hc[s]%P for s in m))
        need(O==E+g and transfer(H,g,E,O,z)==N and subset(B)[0][-1]==N and sum(poly.values())==N,N)
        need(len(W[N])==item['apex_word_count'])
        if w:=item['apex_example']:
            need(w in W[N] and subset(apex(w))[0][-1]==N)
        if A:=item.get('small_rows'):
            need(subset(A)[0][-1]==N)
        if C:=item.get('cyclic_modules'):
            (u,R,z,A),(v,S,D)=C
            need(two(A)==(u,R) and two(D)==(v,S))
            for _ in range(z):A=ordinal(A,(0,))
            need((S-v)*(R<<z)+(3*v-S)*u==N and subset(c3(A,D))[0][-1]==N)
        mode,s0,data=item['family']
        if mode=='root':
            row,r,z=data; T=build_root(row,r,z)
            h,g,e,o=profile(T,r)
            need(h==N and o==e+g)
            alpha=2*e
        elif mode=='apex':
            need(data in W[N])
            alpha=int(data,2)
        elif mode=='cyclic':
            need(bool(C))
            alpha=(S-v)*(R<<z)
        else:raise ValueError(mode)
        delta=N-alpha
        need((18<<s0)<=alpha<=(28<<s0))
        k0=0 if s0==s else 6
        need((18<<(s0+k0))<=(alpha<<k0)+delta<=(28<<(s0+k0)))
        need((alpha*(pow(2,6,21)-1))%21==0)
        print(N, len(W[N]),len(B),len(core),len(m),mode,alpha,delta,s0,k0,
              ' + '.join(f'{c}q^{i}' for i,c in sorted(poly.items(),reverse=True)))
    need(profile((2,20,17,7,105,15,47),1)==(49,9,24,33))
    need(two((126,60,120,112,96,64,2))==(17,448))
    need(two((24,1,3,6,14))==(11,100))
    X=(6,12,96,21,7,27,59); Y=(126,124,88,112,96,68,0)
    for k in range(4):
        A=build_root(X,2,(1<<(3+k))-1)
        B=c3(Y,(2,0)); expect=(768<<k)-5
        need(subset(A)[0][-1]==subset(B)[0][-1]==expect)
        need(sum(v.bit_count()==k+7 for v in A)==1)
        need(sum(v.bit_count()==k+7 for v in B)==2)
        Y=ordinal(Y,(0,))
    need(pack_poly(apex('100110010'))=={10:1,8:24,6:94,4:84,2:14})
    need(prime(877))
    N=6139
    need(all(one('1'+(f'{q:0{n-2}b}' if n>2 else '')+'0')!=N
             for n in range(2,(N-1).bit_length()+2)
             for q in range(1<<(n-2))))
    for A,D,K in ((768,-5,(0,3,15)),(812,49,(1,5,9,21)),(1788,-31,(0,3))):
        for k in K:
            h=(A<<k)+D
            need(h%7==0 and prime(h//7))
    for t in (5,11,17,35,101):
        H=one('10011'+'0'*t)
        need(H==19*(1<<t)-13 and H%7==0 and H%3 and (18<<t)<H<(28<<t))
    need(prime((19*(1<<11)-13)//7) and prime((19*(1<<35)-13)//7))
    for M in (2,3,4,5,7,8,9,12,21,32,49,62,63,72,105):
        A,B=block_pair(M);s=3+2*len('10011'+'0'*7)+(M//gcd(2,M)-1)*len(A)+5
        L=lambda w:int(w,2)
        R=lambda w:sum((b=='0')<<i for i,b in enumerate(w))
        need((L(A)-L(B))%M==(R(A)-R(B))%M==0)
        for r in (1,7,13):
            make_word(s,M,r)
    print('35 prime targets, exact one-apex census, cyclic and rooted rows, modular construction: PASS')


if __name__=='__main__':
    replay()
