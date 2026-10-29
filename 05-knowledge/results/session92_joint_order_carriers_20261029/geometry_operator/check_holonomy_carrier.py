#!/usr/bin/env python3
"""Exact signed boundary incidence: independent construction, rational ranks."""
from itertools import combinations, product
from fractions import Fraction
import json
CHECKS=0
def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok: raise RuntimeError(message)
def rank(a):
    a=[[Fraction(x) for x in row] for row in a]
    if not a:return 0
    r=0
    for c in range(len(a[0])):
        p=next((p for p in range(r,len(a)) if a[p][c]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r]
        v=a[r][c];a[r]=[x/v for x in a[r]]
        for i in range(r+1,len(a)):
            v=a[i][c]
            if v:a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

def incidence(m, hyperedges):
    # State digits: -1 empty, 0 first port, 1 second port.
    states=list(product((-1,0,1),repeat=m));index={s:i for i,s in enumerate(states)}
    rows=[]
    for i in range(m):
        for rest in product((-1,0,1),repeat=m-1):
            s=list(rest);s.insert(i,0);t=s[:];t[i]=1
            exponent=sum(all(s[j]==1 for j in F if j!=i) for F in hyperedges if max(F)==i)
            row=[0]*len(states);row[index[tuple(s)]]=1;row[index[tuple(t)]]=(-1)**exponent
            rows.append(row)
    return states,rows

def faces(m,hyperedges):
    return sum(not any(all(S>>j&1 for j in F) for F in hyperedges) for S in range(1<<m))

def verify(m,hyperedges,gram=False):
    states,R=incidence(m,hyperedges)
    expected=3**m-faces(m,hyperedges)
    check(rank(R)==expected,('rank',m,hyperedges))
    if gram:
        D=[[sum(row[i]*row[j] for row in R) for j in range(len(states))] for i in range(len(states))]
        check(rank(D)==expected,('gram',m,hyperedges))
    # Disjoint sectors; local cube rank checked without invoking global rank formula.
    for S in range(1<<m):
        cols=[k for k,s in enumerate(states) if sum((x!=-1)<<j for j,x in enumerate(s))==S]
        sub=[[row[k] for k in cols] for row in R if any(row[k] for k in cols)]
        kernel=len(cols)-rank(sub)
        check(kernel==int(not any(all(S>>j&1 for j in F) for F in hyperedges)),('sector',S))
    return expected

def full_operator_control(m,hyperedges):
    states,R=incidence(m,hyperedges)
    masks=[sum(1<<(2*j+x) for j,x in enumerate(s) if x!=-1) for s in states]
    feasible=set(masks);n=1<<(2*m)
    B=[[int((i^j).bit_count()==1) for j in range(n)] for i in range(n)]
    # Independently translate row construction into physical cube edges.
    for s in states:
        mask=sum(1<<(2*j+x) for j,x in enumerate(s) if x!=-1)
        for i in range(m):
            if s[i]==1:
                y=mask|(1<<(2*i))
                exponent=sum(all(s[j]==1 for j in F if j!=i) for F in hyperedges if max(F)==i)
                B[mask][y]=B[y][mask]=(-1)**exponent
    forbidden=[j for j in range(n) if j not in feasible]
    for i,s in enumerate(masks):
        for j,t in enumerate(masks):
            leakage=sum(B[s][y]*B[y][t] for y in forbidden)
            check(leakage==sum(row[i]*row[j] for row in R),'physical incidence')
            check(sum(B[s][y]*B[y][t] for y in range(n))-sum(B[s][y]*B[y][t] for y in masks)==leakage,'block identity')
    return B

def main():
    cases=0
    for m in range(4):
        subsets=[F for r in range(2,m+1) for F in combinations(range(m),r)]
        for mask in range(1<<len(subsets)):
            H=[F for i,F in enumerate(subsets) if mask>>i&1]
            verify(m,H,gram=m<=2);cases+=1
    # All graph-valued twists on four modules, and higher-order controls.
    edges=list(combinations(range(4),2))
    for mask in range(64):
        verify(4,[e for i,e in enumerate(edges) if mask>>i&1]);cases+=1
    for m in (4,5):
        verify(m,[tuple(range(m))]);cases+=1
    for m in (2,3):
        base=full_operator_control(m,[]);twisted=full_operator_control(m,[tuple(range(m))])
        check(all(abs(a)==abs(b) for x,y in zip(base,twisted) for a,b in zip(x,y)),'magnitudes')
        check(sum(a!=b for x,y in zip(base,twisted) for a,b in zip(x,y))==2,'one Hermitian edge')
        masks=[s for s in range(1<<(2*m)) if all(((s>>(2*j))&3)!=3 for j in range(m))]
        check(all(base[s][t]==twisted[s][t] for s in masks for t in masks),'same legal operator')
        difference=[[sum(twisted[s][y]*twisted[y][t]-base[s][y]*base[y][t] for y in range(len(base))) for t in masks] for s in masks]
        check(rank(difference)==2,'rank gap is not difference rank')
    two=[verify(2,[]),verify(2,[(0,1)])]
    three=[verify(3,[]),verify(3,[(0,1,2)])]
    print(json.dumps({'checks':CHECKS,'families':cases,'two_modules_ranks':two, 'three_modules_ranks':three},sort_keys=True))
if __name__=='__main__':main()
