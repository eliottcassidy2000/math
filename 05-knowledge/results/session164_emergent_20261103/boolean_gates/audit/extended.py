#!/usr/bin/env python3
"""Independent complementary-ear, substitution, and full-observer audit."""
from independent import *
from math import gcd

def left(a,q):return apex(a,((1<<len(a))-1)^(1<<q))
def right(a,s):return apex(a,1<<s)
def starts(a):
    n=len(a)
    @lru_cache(None)
    def go(p,rem):
        if not rem:return 1
        return sum(go(v,rem^(1<<v)) for v in range(n) if rem>>v&1 and a[p][v])
    return tuple(go(v,((1<<n)-1)^(1<<v)) for v in range(n))
def profile(a,q):
    n=len(a);r=[v for v in range(n) if a[v][q]];need(len(r)==1,'guard')
    vs=[v for v in range(n) if v!=q];st=starts(a)
    return sum(st),sum(st)-st[q],starts(sub(a,vs))[vs.index(r[0])]
def dual(a):return tuple(zip(*a))
def obs(a,p):
    h=H(a);b=apex(a,1<<p);j=H(b);k=H(apex(b,1<<len(a)))
    need((j-h)%2==0 and (k-j)%2==0,'observer parity')
    return h,(j-h)//2,(k-j)//2

def protected(a,q,s):
    n=len(a)
    if q==s or sum(a[v][q] for v in range(n))!=1 or sum(a[s])!=1:return False
    bs=[b for b in range(n) if b!=q and all(v in (q,b) or a[b][v] for v in range(n))]
    cs=[c for c in range(n) if c!=s and all(v in (s,c) or a[v][c] for v in range(n))]
    return any(len({q,b,s,c})==4 for b in bs for c in cs)

def inflate(q,blocks):
    nodes=[(i,j) for i,b in enumerate(blocks) for j in range(len(b))]
    return tuple(tuple(blocks[i][u][v] if i==j else q[i][j] for j,v in nodes) for i,u in nodes)

def gates():return {(x,y):apex(fromrows((8,1,3,6)),2+x+4*y) for x,y in product(range(2),repeat=2)}

def main():
    count=0;paircases=[]
    for n in range(3,6):
        for code in range(1<<comb(n,2)):
            a=tour(n,code)
            if len(sccs(a))!=1:continue
            for q in range(n):
                if sum(a[v][q] for v in range(n))!=1:continue
                h,u,d=profile(a,q);b=left(a,q)
                need(profile(b,n)==(h+u+d,h+d,h-u),'C transfer')
                need(profile(b,q)==(h+u+d,2*u,0),'D transfer');count+=1
                for s in range(n):
                    if protected(a,q,s) and len(paircases)<40:paircases.append((a,q,s))
    print('Independent C/D transfer cases, strong n=3..5:',count,flush=True)
    for a,q,s in paircases:
        h=H(a);U=H(left(a,q))-h;V=H(right(a,s))-h;W=H(right(left(a,q),s))-h-U-V
        lr=right(left(a,q),s);rl=left(right(a,s),q)
        need(lr==sub(rl,list(range(len(a)))+[len(a)+1,len(a)]),'commutation')
        for k,l in [(0,2),(2,0),(1,2),(2,1),(2,2)]:
            b=a
            for _ in range(k):b=left(b,q)
            for _ in range(l):b=right(b,s)
            need(protected(b,q,s) and len(sccs(b))==1,'protected strong')
            need(H(b)==h+((1<<k)-1)*U+((1<<l)-1)*V+((1<<k)-1)*((1<<l)-1)*W,'bilinear')
    print('Independent protected-pair bilinear cases:',len(paircases),'with five lattice points each',flush=True)
    for rows,q,p,defect in [((8,1,3,6),2,3,1),((20,1,2,7,14),3,4,2)]:
        a=fromrows(rows);h,u,d=profile(a,q);need(d==defect,'initial defect');b=a;g=q;steps=[]
        for i,ch in enumerate('CDDC'):
            b=left(b,g)
            if ch=='C':g=len(b)-1
            if i%2:steps.append(H(b))
        need(steps==[3*h+2*d,9*h+8*d],'two equal offsets')
        need(profile(b,g)==(9*h+8*d,5*h+4*d,h),'CDDC final')
        old=obs(a,p);need(obs(b,p)==(9*old[0]+8*d,9*old[1],9*old[2]),'named observer')
    seeds=[((8,1,19,6,43,15),4,0,(17,12,12,8)),((20,1,2,39,14,87,31),5,1,(31,22,22,16))]
    for rows,q,s,(h,U,V,W) in seeds:
        a=fromrows(rows)
        for k,l in product(range(4),repeat=2):
            b=a
            for _ in range(k):b=left(b,q)
            for _ in range(l):b=right(b,s)
            want=h+((1<<k)-1)*U+((1<<l)-1)*V+((1<<k)-1)*((1<<l)-1)*W
            need(H(b)==want and protected(b,q,s) and len(sccs(b))==1,'named family')
        print('Independent named family: seed',h,'all 0<=k,l<=3 PASS',flush=True)
    gs=gates()
    for s in [1,2]:
        for bits in product(range(2),repeat=2*s):
            blocks=[gs[bits[2*i],bits[2*i+1]] for i in range(s)]
            if s%2:blocks.append(((False,),))
            Y=ordinal(blocks);Z=wrap(Y);k=sum(bits[2*i+1]*(1-bits[2*i]) for i in range(s));r=s-k
            P=9**r*11**k;R=2**(s%2)*84**r*100**k
            for p in [len(Y),len(Y)+1]:
                v=obs(Z,p);need(v3(gcd(*v))==r,'full observer valuation')
                if p==len(Y):need(v==(P+R,R//2,P+R//2),'exact a observer')
        print('Independent full-observer valuation, both wrapper ports: s=',s,'PASS',flush=True)
    q=fromrows((8,1,3,6));base=obs(q,3)
    for bits in product(range(2),repeat=4):
        z=inflate(q,[gs[bits[0],bits[1]],gs[bits[2],bits[3]],((False,),),((False,),)])
        v=obs(z,11);factor=4*bits[1]*(1-bits[0])*bits[3]*(1-bits[2])
        need(all((a-factor*b)%3==0 for a,b in zip(v,base)),'uniform substitution')
        need(len(sccs(z))==1,'inflation strong')
    print('Independent uniform-substitution full observers: all 16 two-gate settings PASS',flush=True)
    print('PASS: extended independent controls',flush=True)

if __name__=='__main__':main()
