#!/usr/bin/env python3
"""Independent scalar-operation audit; never imports candidate verifier."""
import hashlib,json,itertools
from pathlib import Path
from functools import lru_cache
ROOT=Path(__file__).resolve().parent

def check(ok,label):
    if not ok: raise RuntimeError(label)
def physical(x,c):
    D,E,U,d=x
    return {'O':(D,0,U+U,d),'Z':(D-E+U,d,U+U-d,d),'J':(D-E+U,D-E,D+U,D)}[c]
def lift(y):
    a,f,e,p=y
    return (a+2*e+p,e,a+f+2*e+2*p,e+p)
def coords(x):
    D,E,U,d=x
    return (D-E-d,U-D+E-d,E,d-E)
def action(y,c):
    a,f,e,p=y
    return {'O':(a+e,a+2*f+e+2*p,0,e+p),'Z':(2*a+f+e+p,f,e+p,0),'J':(f+p,0,a+e+p,e)}[c]
def pairing(s,t):
    D,E,U,d=s; A,B,V,b=t
    return D*b+E*B+U*V+d*A-E*b-U*b-d*B-d*V+d*b

def apply(s,w):
    for c in w:s=physical(s,c)
    return s
@lru_cache(None)
def g(n):
    if n==1:return 1
    if n%2==0:return g(n//2)
    h=(n+1)//2
    return g(h)+h-(h%2)
def tau(n):return (g(n),n%2,n,1)
def unit_word(n):
    if n==2:return ''
    if n%2==0:return unit_word(n//2)+'O'
    return unit_word((n+1)//2)+'Z'
B=[tuple(int(i==j) for i in range(4)) for j in range(4)]
menu=json.loads((ROOT/'input/INPUT-MENU.json').read_text())
check(len(menu)==12 and len({e['word'] for e in menu})==12,'twelve distinct words')
check(sum(e['word']=='' for e in menu)==1,'empty present once')
expected={'O':[[1,1,2,2],[1,1,2,2],[1,1,2,2],[1,1,3,2]],'Z':[[1,1,1,1],[1,1,1,1],[2,2,2,3],[2,2,2,2]],'J':[[1,0,2,1],[0,0,1,1],[2,1,3,3],[1,1,3,2]]}
residual={}; matrices={}
for c in 'OZJ':
    for b in B:
        check(coords(lift(b))==b,'inverse on basis')
        check(physical(lift(b),c)==lift(action(b,c)),'conjugate physical action')
        check(min(action(b,c))>=0,'cone coefficients')
        check(sum(action(b,c)[2:])>=sum(b[2:]),'d coefficient monotonicity')
    R=[[pairing(lift(x),physical(lift(y),c))-pairing(lift(x),lift(y)) for y in B] for x in B]
    check(R==expected[c],'literal residual '+c)
    check(all(R[i][j]>=2*int(i>=2 and j>=2) for i in range(4) for j in range(4)),'residual strict lower bound')
    residual[c]=R
for e in menu:
    w=e['word']; A=[[pairing(x,apply(y,w)) for y in B] for x in B]
    check(A==e['matrix'],'word matrix '+w)
    matrices[w]=A
sources=[apply((1,0,2,1),w) for k in range(7) for w in itertools.product('OZJ',repeat=k)]
check(len(sources)==1093,'source universe')
for s in sources:
    y=coords(s)
    check(min(y)>=0 and s[3]>=1 and lift(y)==s,'source boundary')
inc=0
for s in sources:
    for t in sources[::17]:
        for c in 'OZJ':
            check(pairing(s,physical(t,c))-pairing(s,t)>=2*s[3]*t[3],'strict increment finite')
            inc+=1
for n in range(2,4097):
    check(apply((1,0,2,1),unit_word(n))==tau(n),'all unit sources constructible')
    check(min(coords(tau(n)))>=0,'unit source boundary')
mins=0
for u in range(2,98):
    for v in range(2,98):
        baseline=pairing(tau(u),tau(v))
        check(baseline>=1,'output natural')
        for e in menu:
            val=pairing(tau(u),apply(tau(v),e['word']))
            check(val-baseline>=2*len(e['word']),'full menu minimum')
            check((val==baseline)==(e['word']==''),'unique minimum template')
        if u==2:check(baseline==g(v)+v-v%2,'fixed row')
        mins+=1
formula=0
for a in range(2,26):
    for b in range(1,26):
        n=2**(a+b)+2**b+1
        check(g(n)==n-2**(a-1),'sparse recurrence identity')
        formula+=1
separators=0
for k in range(1,41):
    for m in range(1,41):
        letters=[(1,1)]+[(0,0)]*k+[(1,0),(0,1)]+[(0,0)]*(m-1)+[(1,1)]
        u=int(''.join(str(x[0]) for x in letters),2)
        y=int(''.join(str(x[1]) for x in letters),2)
        check(u==2**(k+m+2)+2**(m+1)+1 and y==2**(k+m+2)+2**m+1,'word decoding')
        check((g(u)==y)==(k==m),'distinguishing suffix')
        separators+=1
print(json.dumps({'status':'PASS','method':'scalar operations and bilinear coefficient extraction, independent of candidate imports','menu_sha256':hashlib.sha256((ROOT/'input/INPUT-MENU.json').read_bytes()).hexdigest(),'physical_states':len(sources),'increment_checks':inc,'unit_sources':4095,'union_fibers':mins,'sparse_formula_checks':formula,'separator_pairs':separators,'residuals':residual,'regenerated_menu':matrices},indent=2,sort_keys=True))
