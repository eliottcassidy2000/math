#!/usr/bin/env python3
"""Weighted fixed-set moments and exact joint orbit masses; standard library."""
from fractions import Fraction
from itertools import product
from math import comb,prod
from collections import Counter
import json

CHECKS=0
def ck(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(msg)

def compose(a,b):return tuple(a[x] for x in b)
def neg(k):return tuple(sum((- (v//3**i)%3)*3**i for i in range(k)) for v in range(3**k))
def aut(k):
    if not k:return [(0,)]
    child=aut(k-1);s=3**(k-1)
    return [tuple(((b+r)%3)*s+abc[b][j] for b in range(3) for j in range(s))
            for r in range(3) for abc in product(child,repeat=3)]
def cent(k):
    if not k:return [(0,)]
    s=3**(k-1);t=neg(k-1)
    return [tuple(c)+tuple(s+x for x in a)+tuple(2*s+x for x in compose(t,compose(a,t)))
            for c in cent(k-1) for a in aut(k-1)]

def moments(kind,k,w,rmax):
    if not k:return [Fraction(w[0]**r) for r in range(rmax+1)]
    s=3**(k-1)
    if kind=='A':
        m=[moments('A',k-1,w[b*s:(b+1)*s],rmax) for b in range(3)]
        return [Fraction(1)]+[sum(comb(r,i)*comb(r-i,j)*m[0][i]*m[1][j]*m[2][r-i-j]
                                  for i in range(r+1) for j in range(r-i+1))/3
                               for r in range(1,rmax+1)]
    t=neg(k-1)
    left=moments('C',k-1,w[:s],rmax)
    paired=moments('A',k-1,[w[s+j]+w[2*s+t[j]] for j in range(s)],rmax)
    return [sum(comb(r,i)*left[i]*paired[r-i] for i in range(r+1)) for r in range(rmax+1)]

def orbit_mass(kind,k,ports,weights):
    """Sum product_i weights[i][y_i] over DISTINCT tuples in the joint orbit."""
    if not ports:return 1
    if not k:return prod(w[0] for w in weights)
    s=3**(k-1);blocks=[[i for i,x in enumerate(ports) if x//s==b] for b in range(3)]
    if kind=='A':
        total=0
        for shift in range(3):
            term=1
            for b,inds in enumerate(blocks):
                dest=(b+shift)%3
                term*=orbit_mass('A',k-1,[ports[i]%s for i in inds],
                                 [weights[i][dest*s:(dest+1)*s] for i in inds])
            total+=term
        return total
    t=neg(k-1);zero=blocks[0];fold=[i for i in range(len(ports)) if ports[i]//s]
    first=orbit_mass('C',k-1,[ports[i]%s for i in zero],[weights[i][:s] for i in zero])
    p=[];w=[]
    for i in fold:
        b=ports[i]//s
        p.append(ports[i]%s if b==1 else t[ports[i]%s])
        w.append(weights[i][s:2*s] if b==1 else [weights[i][2*s+t[j]] for j in range(s)])
    return first*orbit_mass('A',k-1,p,w)

def orbit_size(kind,k,ports):
    return orbit_mass(kind,k,ports,[[1]*3**k for _ in ports])

def berggren(k,root=(1,2)):
    row=[root]
    for _ in range(k):
        # digit0=M, digit1=L, digit2=R: negation fixes M and exchanges L,R.
        row=[child for m,n in row for child in [(n,2*n+m),(n,2*n-m),(m,2*m+n)]]
    return [(n*n-m*m,2*m*n,n*n+m*m) for m,n in row]

def word_parameters(word,root):
    m,n=root
    for digit in word:
        m,n=[(n,2*n+m),(n,2*n-m),(m,2*m+n)][digit]
    return m,n

def main():
    rows=[]
    for k in range(1,4):
        n=3**k;G=cent(k);triples=berggren(k)
        for a,b,c in triples:
            ck(a*a+b*b==c*c,'Pythagorean equation')
            import math
            ck(math.gcd(a,b)==1 and a%2==1 and b%2==0,'primitive canonical triple')
        weights=[c for a,b,c in triples]
        m=moments('C',k,weights,4)
        for r in range(5):
            literal=Fraction(sum(sum(weights[i] for i in range(n) if g[i]==i)**r for g in G),len(G))
            ck(m[r]==literal,'weighted fixed-set moment')
        # Every joint orbit, using actual source-product weights, not quotient averages.
        orbit_counts=[]
        for r in range(1,4):
            seen=set();total=0;count=0
            ws=[[triples[v][i] for v in range(n)] for i in range(r)]
            for ports in product(range(n),repeat=r):
                if ports in seen:continue
                orb={tuple(g[v] for v in ports) for g in G};seen.update(orb);count+=1
                actual=sum(prod(ws[i][v] for i,v in enumerate(t)) for t in orb)
                predicted=orbit_mass('C',k,ports,ws)
                ck(predicted==actual,'source-preserving joint orbit mass')
                ck(orbit_size('C',k,ports)==len(orb),'joint orbit size')
                total+=predicted
            ck(total==prod(sum(w) for w in ws),'total named source mass')
            orbit_counts.append(count)
        rows.append(dict(k=k,vertices=n,centralizer=len(G),joint_orbit_counts=orbit_counts,
                         hypotenuse_fixed_moments=list(map(str,m))))
    # Two equal marginals can place all mass in disjoint joint-charge orbits.
    laws=[[(j,(charge-j)%3) for j in range(3)] for charge in [0,1]]
    for side in [0,1]:ck(Counter(p[side] for p in laws[0])==Counter(p[side] for p in laws[1]),'equal marginals')
    ck({(j+l)%3 for j,l in laws[0]}!={(j+l)%3 for j,l in laws[1]},'different joint charge')
    # The outer mirror is an exact arithmetic intertwiner only with its root frame.
    for k in range(6):
        for word in product(range(3),repeat=k):
            m,n=word_parameters(word,(1,2));mirror=tuple((-d)%3 for d in word)
            u,v=word_parameters(mirror,(1,3))
            ck((u,v)==(n-m,n+m),'Berggren root-frame intertwiner')
            ck((v*v-u*u,2*u*v,v*v+u*u)==(4*m*n,2*(n*n-m*m),2*(n*n+m*m)), 'leg swap and factor two')
    print(json.dumps(dict(scope='uniform group moments; separate product-source orbit masses',
                         rows=rows,active_checks=CHECKS),indent=2,sort_keys=True))

if __name__=='__main__':main()
