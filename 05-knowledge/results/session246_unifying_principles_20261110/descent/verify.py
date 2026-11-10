#!/usr/bin/env python3
"""Standalone exact structural checks. No external packages or repo imports."""
from fractions import Fraction as Q
from itertools import product
import json, sys

def check(v, msg):
    if not v: raise RuntimeError(msg)

def step(s, c):
    D,E,U,d=s
    if c=='O': return D,0,2*U,d
    if c=='Z': return D-E+U,d,2*U-d,d
    if c=='J': return D-E+U,D-E,D+U,D
    raise ValueError(c)

def state(w):
    s=(1,0,2,1)
    for c in w: s=step(s,c)
    return s

def graph(w):
    rows=[2,0]; root=0
    for c in w:
        n=len(rows); rows.append(0)
        for i in range(n):
            if c!='O' and i==root: rows[n]|=1<<i
            else: rows[i]|=1<<n
        if c=='J': root=n
    return rows,root

def count(rows, omit=None):
    vs=[v for v in range(len(rows)) if v!=omit]; n=len(vs)
    if not n: return 1,{}
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1<<n):
        for i in range(n):
            if not dp[mask][i]: continue
            for j in range(n):
                if not mask>>j&1 and rows[vs[i]]>>vs[j]&1:
                    dp[mask|1<<j][j]+=dp[mask][i]
    ends={v:dp[-1][i] for i,v in enumerate(vs)}
    return sum(ends.values()),ends

def words(alphabet, cap):
    for l in range(cap+1):
        for w in product(alphabet,repeat=l): yield ''.join(w)

def tail_check(s,w):
    D,E,U,d=s; a=U-d; y=Q(D,U)
    for n,c in enumerate(w,1):
        s=step(s,c); D,E,U,dd=s; y=(y+(c=='Z'))/2
        check(dd==d and U>=d+2**n*a,'growth')
        check(abs(Q(D,U)-y)<=Q(n*d,2**n*a),'shadow')
    return s,y

def main():
    if '--fail' in sys.argv: check(False,'intentional failing sentinel')
    nw=0
    for w in words('OZJ',7):
        s=(1,0,2,1); ds=[]
        for c in w:
            before=s; s=step(s,c)
            check(s[0]>=before[0] and s[3]>=before[3],'monotonicity')
            if c=='J':
                ds.append(s[3]); check(s[0]>=before[0]+before[3],'J increment')
                if len(ds)>=3: check(ds[-1]>=ds[-2]+ds[-3],'Fibonacci growth')
        nw+=1
    ng=0
    for w in words('OZJ',4):
        s=state(w); rows,r=graph(w); D,ends=count(rows); d,_=count(rows,r)
        child,_=graph(w+'Z'); Dz,_=count(child)
        check((D,ends[r],Dz-D+ends[r],d)==s,'independent full profile')
        ng+=1
    roots=['','J','ZJ','OJJ','JOZJ','ZZOJ','JJOJ']
    nt=0
    for root in roots:
        for w in words('OZ',8): tail_check(state(root),w); nt+=1
    controllers=0
    for root in roots:
        initial=state(root); d=initial[3]; a=initial[2]-d; r=8
        for j in range(65):
            t=Q(j,64); A=min(2**r-1,int(2**r*t))
            w=''.join('Z' if A>>i&1 else 'O' for i in range(r))
            out,_=tail_check(initial,w)
            check(abs(Q(out[0],out[2])-t)<=Q(1+r*Q(d,a),2**r),'finite controller')
            controllers+=1
    cycles=[]
    for block in ['OZ','ZO','OOZ','OZZ','OZOZ','ZZOOO']:
        s,y=tail_check(state('ZJ'),block*30)
        A=sum((c=='Z')*2**i for i,c in enumerate(block)); limit=Q(A,2**len(block)-1)
        n=len(block)*30; a=state('ZJ')[2]-state('ZJ')[3]; d=state('ZJ')[3]
        check(abs(Q(s[0],s[2])-limit)<=Q(n*d,2**n*a)+Q(1,2**n),'cycle limit')
        cycles.append([block,str(limit)])
    # Dense-path construction checked via increasingly fine rational nets, not float plots.
    dense=''.join(w for m in range(1,7) for w in words('OZ',m) if w)
    seen={}; s=(1,0,2,1); y=Q(1,2); n=0
    for c in dense:
        n+=1; s=step(s,c); y=(y+(c=='Z'))/2
        if n>=32:
            bucket=min(31,int(y*32)); seen[bucket]=n
        check(abs(Q(s[0],s[2])-y)<=Q(n,2**n),'dense shadow')
    check(len(seen)==32,'dense rational net')
    print(json.dumps({'status':'PASS','ternary_words':nw,'independent_graph_profiles':ng,
      'binary_tails':nt,'finite_controllers':controllers,'periodic_limits':cycles,'dense_prefix_length':len(dense),
      'dyadic_bins_hit':len(seen),'scope':'Structural controls, not scalar coverage census'},sort_keys=True,indent=2))
if __name__=='__main__': main()
