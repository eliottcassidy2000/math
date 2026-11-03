#!/usr/bin/env python3
"""Exact controls for the C/D reset language; standard library only."""
from itertools import combinations, product
from pathlib import Path
import json

I=((1,0,0),(0,1,0),(0,0,1))
C=((0,1,0),(1,1,1),(1,0,0))
D=((1,0,1),(0,2,0),(0,0,0))

def need(c,msg):
    if not c: raise RuntimeError(msg)
def mul(a,b): return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))) for i in range(len(a)))
def power(a,n):
    p=I
    for _ in range(n): p=mul(a,p)
    return p
def transfer(w):
    m=I
    for c in w: m=mul(C if c=='C' else D,m)
    return m
def apply(m,v): return tuple(sum(a*b for a,b in zip(row,v)) for row in m)
def det(a):
    x,y,z=a
    return x[0]*(y[1]*z[2]-y[2]*z[1])-x[1]*(y[0]*z[2]-y[2]*z[0])+x[2]*(y[0]*z[1]-y[1]*z[0])
def rank(a):
    if det(a): return 3
    for i,j in combinations(range(3),2):
        for k,l in combinations(range(3),2):
            if a[i][k]*a[j][l]!=a[i][l]*a[j][k]: return 2
    return int(any(any(row) for row in a))
def outer(v,r): return tuple(tuple(a*b for b in r) for a in v)
def add(rows,out):
    n=len(rows)
    return tuple(r|((0 if out>>i&1 else 1)<<n) for i,r in enumerate(rows))+(out,)
def build(rows,q,w):
    for c in w:
        rows=add(rows,((1<<len(rows))-1)^(1<<q))
        if c=='C': q=len(rows)-1
    return rows,q
def subset(rows,vs): return tuple(sum(1<<j for j,v in enumerate(vs) if rows[u]>>v&1) for u in vs)
def starts(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if mask>>i&1:
                rem=mask^(1<<i)
                if rem: dp[mask][i]=sum(dp[rem][j] for j in range(n) if rows[i]>>j&1)
    return tuple(dp[-1])
def state(rows,q):
    v=starts(rows); incoming=[i for i in range(len(rows)) if rows[i]>>q&1]
    need(len(incoming)==1,'guard is not unique-in-neighbor')
    vs=[i for i in range(len(rows)) if i!=q]
    d=starts(subset(rows,vs))[vs.index(incoming[0])]
    return (v[q],sum(v)-v[q],d)
def query_matrix(rows,q,p):
    cols=[]
    for i in range(3):
        cols.append(state(rows,q))
        rows=add(rows,1<<(p if i==0 else len(rows)-1))
    return tuple(zip(*cols))
def observer(v):
    h=tuple(v[0][i]+v[1][i] for i in range(3))
    need((h[1]-h[0])%2==0 and (h[2]-h[1])%2==0,'observer parity')
    return (h[0],(h[1]-h[0])//2,(h[2]-h[1])//2)
def strong(rows):
    for i in range(len(rows)):
        seen=1<<i
        while True:
            nxt=seen
            for j in range(len(rows)):
                if seen>>j&1:nxt|=rows[j]
            if seen==nxt:break
            seen=nxt
        if seen!=(1<<len(rows))-1:return False
    return True

def polyadd(a,b): return tuple(a[i]+b[i] for i in range(3))
def polysub(a,b): return tuple(a[i]-b[i] for i in range(3))
def polymul(a,b): return tuple(sum(a[i]*b[k-i] for i in range(k+1)) for k in range(3))
def jetmul(a,b):
    return tuple(tuple(tuple(sum(a[i][k][j]*b[k][l][degree-j] for k in range(3) for j in range(degree+1)) for degree in range(3)) for l in range(3)) for i in range(3))
def slot_controls():
    z=(0,0,0);one=(1,0,0);r=(1,1,0);t=(0,1,0)
    cj=((t,r,z),(one,one,one),(one,z,z));dj=((r,z,one),(z,(2,1,0),z),(z,z,t))
    ij=tuple(tuple(one if i==j else z for j in range(3)) for i in range(3))
    for k in range(41):
        m=dj
        for _ in range(k):m=jetmul(cj,m)
        m=jetmul(dj,m)
        minors=[polysub(polymul(m[i][a],m[j][b]),polymul(m[i][b],m[j][a])) for i,j in combinations(range(3),2) for a,b in combinations(range(3),2)]
        order=next((h for h in range(3) if any(v[h] for v in minors)),3)
        need(order==(1 if k in (1,14) else 0),'two-minor jet order')
        determinant=z
        for j,sign in ((0,1),(1,-1),(2,1)):
            cols=[a for a in range(3) if a!=j]
            minor=polysub(polymul(m[1][cols[0]],m[2][cols[1]]),polymul(m[1][cols[1]],m[2][cols[0]]))
            term=polymul(m[0][j],minor)
            determinant=polyadd(determinant,tuple(sign*x for x in term))
        need(determinant==(0,0,4),'gap determinant jet')
    tests=0
    for length in range(6):
        for letters in product('CD',repeat=length):
            for r in range(6):
                cr=((r-1,r,0),(1,1,1),(1,0,0));dr=((r,0,1),(0,r+1,0),(0,0,r-1));m=I
                for ch in letters:m=mul(cr if ch=='C' else dr,m)
                need(det(m)==r**length*(r*r-1)**letters.count('D'),'full slot determinant')
                tests+=1
    return {'full_slot_determinant_controls':tests,'gap_local_jet_controls':41}

def main():
    need(det(C)==1 and rank(D)==2,'initial ranks')
    pairs=list(combinations(range(3),2))
    E=tuple(tuple(C[a][c]*C[b][d]-C[a][d]*C[b][c] for c,d in pairs) for a,b in pairs)
    need(E==((-1,0,1),(-1,0,0),(-1,-1,0)),'exterior matrix')
    E2=mul(E,E);E3=mul(E2,E)
    need(all(E3[i][j]+E2[i][j]+E[i][j]-I[i][j]==0 for i in range(3) for j in range(3)),'exterior polynomial')
    T={0:0,1:1,2:1}
    for n in range(-1,-205,-1):T[n]=T[n+3]-T[n+2]-T[n+1]
    deltas=[]
    for k in range(201):
        a=mul(D,power(C,k));delta=a[0][0]*a[1][1]-a[0][1]*a[1][0]
        need(delta==-2*T[-k-3],f'reflection identity {k}')
        deltas.append(delta)
    need([k for k,a in enumerate(deltas) if a==0]==[1,14],'finite zero controls')
    w0=(1,2,0);w1=(47,112,0);l0=(1,2,1);l1=(47,112,47)
    R0=transfer('DCD');R1=transfer('D'+'C'*14+'D')
    need(R0==outer(w0,l0) and R1==outer(w1,l1),'reset outer products')
    need(mul(R1,R1)==tuple(tuple(14753*x for x in row) for row in R1),'new-ray multiplier')
    # Independent exact minor tests, not the rank-language recurrence.
    census=[];front=[('',I)];seen=0
    for depth in range(19):
        hist=[0,0,0,0];nxt=[]
        for w,a in front:
            actual=rank(a)
            expect=3 if 'D' not in w else (1 if ('DCD' in w or 'D'+'C'*14+'D' in w) else 2)
            need(actual==expect,'rank language '+w)
            hist[actual]+=1;seen+=1
            if depth<18:
                # Sparse left multiplication to keep the exhaustive test compact.
                c=(a[1],tuple(a[0][j]+a[1][j]+a[2][j] for j in range(3)),a[0])
                d=(tuple(a[0][j]+a[2][j] for j in range(3)),tuple(2*x for x in a[1]),(0,0,0))
                nxt.extend(((w+'C',c),(w+'D',d)))
        census.append(hist[1:]);front=nxt
    # Coefficients from R2(z)=z/((1-z)(1-2z+z²-z³+z^15-z^16)).
    q=[1,-3,3,-2,1]+[0]*10+[1,-2,1]
    a=[]
    for n in range(19):
        a.append((1 if n==1 else 0)-sum(q[j]*a[n-j] for j in range(1,min(len(q)-1,n)+1)))
        need(a[n]==census[n][1],'rank-two generating series')
    seed=(2,4,1);q=0;p=1;V=query_matrix(seed,q,p)
    need(V==((1,2,3),(2,3,6),(0,1,2)) and det(V)==-2,'physical full-rank seed')
    controls=0;examples=[]
    for n in range(6):
        for letters in product('CD',repeat=n):
            w=''.join(letters);rows,r=build(seed,q,w)
            need(strong(rows),'strongness '+w)
            Y=query_matrix(rows,r,p)
            need(Y==mul(transfer(w),V),'physical query transfer '+w)
            controls+=1
            if w in ('DCD','CDD'):examples.append({'word':w,'rows':rows,'guard':r,'query_matrix':Y,'observer':observer(Y)})
    R1V=mul(R1,V)
    need(observer(R1V)==(43089,16377,34185),'long reset observer')
    longrows,longq=build(seed,0,'D'+'C'*14+'D')
    need(strong(longrows),'long reset strongness')
    # Locking: rank-one query state and all 127 suffixes through length six.
    locked=mul(R0,V);v=observer(locked)
    for n in range(7):
        for letters in product('CD',repeat=n):
            obs=observer(mul(transfer(''.join(letters)),locked))
            need(all(obs[i]*v[j]==obs[j]*v[i] for i,j in combinations(range(3),2)),'locked observer ray')
    out={**slot_controls(),'rank_words_checked':seen,'rank_census_by_depth_0_through_18':census,'gap_determinants_0_through_20':deltas[:21],
         'finite_gap_zero_check_through':200,'gap_zeros_in_control':[1,14],'all_depth_zero_source':'Mignotte–Tzanakis (1991), as restated by Bilu et al. (2024), p.1333; not established by finite search',
         'R0':R0,'R14':R1,'physical_query_controls':controls,'physical_seed_query_matrix':V,'examples':examples,
         'long_reset':{'word':'D'+'C'*14+'D','rows':longrows,'guard':longq,'order':len(longrows),'formula_query_matrix':R1V,'formula_observer':observer(R1V)},
         'new_ray_return_multiplier':14753,'status':'PASS'}
    Path('results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
