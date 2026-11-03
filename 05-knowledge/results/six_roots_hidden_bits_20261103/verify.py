#!/usr/bin/env python3
"""Exact finite controls for BOUNDARIES and PORT_MEMORY; not universal coverage."""
import itertools,json,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
REV='2af5e7549da1f7a6cba03e42e15500eb060774c2'
def require(ok,label):
    if not ok: raise RuntimeError(label)
def profile(rows):
    n=len(rows); full=(1<<n)-1
    d=[[0]*n for _ in range(1<<n)]
    for v in range(n):d[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            if not d[s][v]:continue
            todo=rows[v]&~s
            while todo:
                bit=todo&-todo;todo-=bit
                d[s|bit][bit.bit_length()-1]+=d[s][v]
    h=[sum(x) for x in d];h[0]=1
    return h[-1],[sum(d[s][q]*h[full^s] for s in range(1,1<<n)) for q in range(n)]
def en(n):
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    rows[0]^=1<<(n-1);rows[-1]|=1
    return rows
def ear(rows,q):
    n=len(rows);a=n;b=n+1;r=rows+[0,0]
    r[q]|=1<<a;r[a]=(1<<b)|(((1<<n)-1)^(1<<q));r[b]=(1<<n)-1
    return r
def mul(x,y):
    a,b=x;c,d=y
    return a*c+b*d,a*d+b*c+b*d
def conj(x):a,b=x;return a+b,-b
def norm(x):a,b=x;return a*a+a*b-b*b
def gp(k):
    z=(1,0);p=(0,1) if k>=0 else (-1,1)
    for _ in range(abs(k)):z=mul(z,p)
    return z
def main():
    # Independent brute permutation/cut audit lives in audit.py.
    profiles=[]
    for n in range(3,15):
        h,k=profile(en(n));f=[h]+[2**q+3**(q-1)*2**(n-1-q)-(q==n-2) for q in range(1,n-1)]+[2**(n-1)-1]
        require(h==2**(n-2)+1 and k==f,('E_n formula',n,k,f))
        require(n==3 or len(set(k))==n,('port separation',n))
        if n<=9:
            for q in range(n):require(profile(ear(en(n),q))[0]==h+2*k[q],('literal ear',n,q))
        profiles.append({'n':n,'H':h,'K':k,'distinct':len(set(k))})
    seeds=json.loads(subprocess.check_output(['git','show',REV+':05-knowledge/results/triple_reentry_interfaces_20261103/scalar_seeds.json']))
    actual=[]
    for h,s in sorted(seeds.items(),key=lambda p:int(p[0])):
        got,k=profile(s['rows']);require(got==int(h),'seed H')
        outputs=[profile(ear(s['rows'],q))[0] for q in range(len(k))]
        require(outputs==[got+2*v for v in k],'seed ear')
        actual.append({'H':got,'rows':s['rows'],'K':k,'ear_counts':outputs,'bits_needed_for_one_ear':(len(set(k))-1).bit_length()})
    # Full b-ary heap, root0. E excludes root, children bu+1,...,bu+b.
    checks=0
    for b in range(2,6):
        for mask in range(1<<10):
            e={i+1 for i in range(10) if mask>>i&1}
            roots={0}|{b*u+j for u in e for j in range(1,b+1) if b*u+j not in e}
            internal=sum((u-1)//b in e for u in e)
            require(len(roots)==1+(b-1)*len(e)+(len(e)-internal),'forest boundary');checks+=1
    forbidden={7,21};roots={1}|{3*u+2*j for u in forbidden for j in range(3) if 3*u+2*j not in forbidden}
    require(roots=={1,23,25,63,65,67},'six roots')
    mersenne=[(a,b) for a in range(1,41) for b in range(1,41) if 3**a*(2**b-1)==2**(2*b)-1]
    require(mersenne==[(1,1),(2,3)],'Mersenne audit')
    units={gp(k) for k in range(-30,31)};units|={(-a,-b) for a,b in list(units)}
    for a in range(-100,101):
        for b in range(-100,101):
            if abs(norm((a,b)))==1:require((a,b) in units,'unit box')
    for k in range(-50,51):
        z=gp(k);w=gp(-k);sgn=(-1)**k
        require(conj(z)==tuple(int(sgn*x) for x in w),'golden conjugation')
        require(norm(z)==(-1)**abs(k),'golden norm')
        require(mul(z,w)==(1,0),'golden inverse')
        for l in range(-15,16):require(mul(z,gp(l))==gp(k+l),'golden clock')
    require(norm((1,2))== -1,'norm positive control')
    require(norm((-1,2))== -5,'addition escape sqrt5')
    A=[n for n in range(1,10000,2) if n not in (7,21)]
    B=[n for n in range(1,10000) if n in (1,15,16,17,23) or n>=25]
    rank_map=dict(zip(A,B))
    require([rank_map[h] for h in sorted(roots)]==[1,29,30,49,50,51],'rank six roots')
    for h in A:
        if h>=23:require(rank_map[h]==(h+35)//2,'cofinal rank formula')
        if h>=23 and 3*h+4 in rank_map:
            for j in range(3):require(rank_map[3*h+2*j]==3*rank_map[h]-35+j,'conjugated ternary')
    out={'status':'FINITE-EXACT controls; see hand proofs; no odd-spectrum coverage claim','E_profiles':profiles,'actual_seeds':actual,'deleted_sets_checked':checks,'Mersenne_box':[40,40],'Mersenne_solutions':mersenne,'golden_unit_box':100,'golden_clock_range':[-50,50],'rank_root_images':[rank_map[h] for h in sorted(roots)]}
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: exact ports, literal ears, 4096 deleted sets, golden units, Mersenne exception, rank transport')
if __name__=='__main__':main()
