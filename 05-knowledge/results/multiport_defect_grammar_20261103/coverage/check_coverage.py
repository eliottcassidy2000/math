"""Exact arithmetic grammar controls and literal terminal-C3 transport.

The arithmetic universality claim is proved in COVERAGE.md. Finite arithmetic
and graph universes below are controls, not evidence for spectrum completeness.
"""
from pathlib import Path
from itertools import combinations,permutations
from functools import lru_cache
from math import isqrt
import json,hashlib
HERE=Path(__file__).resolve().parent
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def allowed(n):return n>0 and n%2==1 and n not in (7,21)
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
SEEDS={'triple':[1,23,25,63,65,67],'multiply':[1,3,23,63,67],'both':[1,23,63,67]}
def predecessors(n,mode):
    out=[]
    for c in (2,4):
        if (n-c)%3==0 and allowed((n-c)//3):out.append(('B'+str(c),(n-c)//3))
    if mode in ('triple','both') and n%3==0 and allowed(n//3):out.append(('M3',n//3))
    if mode in ('multiply','both'):
        for a in range(3,isqrt(n)+1,2):
            if n%a==0 and allowed(a) and allowed(n//a):out.append(('P',a,n//a))
    return out
def arithmetic(limit=10001):
    out={}
    for mode,seeds in SEEDS.items():
        heights={n:0 for n in seeds};parents={}
        for n in range(1,limit+1,2):
            if not allowed(n) or n in heights:continue
            choices=[p for p in predecessors(n,mode) if all(a in heights for a in p[1:])]
            need(choices,('finite inverse control',mode,n));p=min(choices,key=lambda p:1+max(heights[a] for a in p[1:]));parents[n]=p;heights[n]=1+max(heights[a] for a in p[1:])
        for seed in seeds:need(not predecessors(seed,mode),('irreducible mandatory seed',mode,seed))
        # Separate forward closure, increasing-value topological sweep.
        seen=set(seeds)
        for h in range(1,limit+1,2):
            if h not in seen:continue
            for c in (2,4):
                if allowed(3*h+c) and 3*h+c<=limit:seen.add(3*h+c)
            if mode in ('triple','both') and allowed(3*h) and 3*h<=limit:seen.add(3*h)
            if mode in ('multiply','both'):
                for a in sorted(seen):
                    if a>limit//h:break
                    if allowed(a*h):seen.add(a*h)
        target={n for n in range(1,limit+1,2) if allowed(n)}
        need(seen==target,('independent forward closure',mode))
        out[mode]=dict(minimal_seeds=seeds,checked_through=limit,covered=len(seen),max_inverse_height=max(heights.values()),exceptions=[7,21],sample_parents={str(n):parents.get(n,['seed']) for n in (3,5,23,25,49,63,65,67,343,170807) if n<=limit})
    # Unary-only words from the actual170807 source, exact ternary addresses.
    h0=170807;levels=[{h0}];rows=[]
    for k in range(11):
        vals=levels[k];need(len(vals)==2**k,'binary distinct descendants')
        lo=3**k*(h0+1)-1;hi=3**k*(h0+2)-2
        need(min(vals)==lo and max(vals)==hi,'exact endpoint envelope')
        for h in vals:
            x=(h+1)//2;prefix=(h0+1)//2
            for j in range(k):need(x%3 in (0,1),'appended ternary digit');x//=3
            need(x==prefix,'fixed inherited ternary prefix')
        rows.append(dict(depth=k,size=len(vals),minimum=lo,maximum=hi))
        if k<10:levels.append({3*h+c for h in vals for c in (2,4)})
    atoms=[]
    for n in range(3,limit+1,2):
        if not allowed(n):continue
        irreducible=not any(n%a==0 and allowed(a) and allowed(n//a) for a in range(3,isqrt(n)+1,2))
        predicted=prime(n) or (n%7==0 and prime(n//7) and n//7!=3) or n in (63,343)
        need(irreducible==predicted,('multiplicative atom classification',n))
        if irreducible:atoms.append(n)
    roots={}
    for n in range(1,limit+1,2):
        if not allowed(n):continue
        h=n
        while h not in SEEDS['triple']:
            d={0:0,1:4,2:2}[h%3];h=(h-d)//3;need(allowed(h),'canonical component ancestor')
        roots[n]=h
    def address(n):
        z=(n-1)//2;digits=[]
        while z:
            d=(z-1)%3+1;digits.append(str(d));z=(z-d)//3
        return ''.join(reversed(digits))
    return dict(grammars=out,unary170807=rows,multiplicative_atoms_through_limit=atoms,ternary_components=dict(roots=SEEDS['triple'],root_addresses={str(n):address(n) for n in SEEDS['triple']},forbidden_addresses={str(n):address(n) for n in (7,21)},counts_through_limit={str(n):sum(r==n for r in roots.values()) for n in SEEDS['triple']}))
def strong(T):
    f=(1<<len(T))-1
    for r in range(len(T)):
        s=1<<r
        while True:
            v=s
            for i in range(len(T)):
                if s>>i&1:v|=T[i]
            if v==s:break
            s=v
        if s!=f:return False
    return True
def join(T,U):
    n=len(T);return tuple(a|(((1<<len(U))-1)<<n) for a in T)+tuple(a<<n for a in U)
def birth(T,cut):
    n=len(T);return tuple(a|(((cut>>i)&1)<<n) for i,a in enumerate(T))+(((1<<n)-1)^cut,)
@lru_cache(None)
def profile(T):
    n=len(T);f=(1<<n)-1;R=tuple(sum(((T[j]>>i)&1)<<j for j in range(n)) for i in range(n));tables=[]
    for U in (T,R):
        e=[[0]*n for _ in range(f+1)]
        for i in range(n):e[1<<i][i]=1
        for m in range(1,f+1):
            for v,a in enumerate(e[m]):
                q=U[v]&~m
                while q:
                    b=q&-q;q-=b;e[m|b][b.bit_length()-1]+=a
        tables.append(e)
    e,s=tables;h=[sum(a) for a in e];h[0]=1
    E=[sum(e[m][r]*h[f^m] for m in range(f+1)) for r in range(n)]
    S=[sum(s[m][r]*h[f^m] for m in range(f+1)) for r in range(n)]
    W=[[sum(e[m][r]*s[f^m][t] for m in range(f+1)) for t in range(n)] for r in range(n)]
    return dict(H=h[f],P2=sum(h[m]*h[f^m] for m in range(f+1)),E1=e[f],S1=s[f],deleted=[h[f^(1<<r)] for r in range(n)],E2=E,S2=S,W2=W)
def brute_h(T):return sum(all(T[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(len(T))))
def graph_controls():
    C=(2,4,1);count=0
    for n in range(1,5):
        edges=list(combinations(range(n),2))
        for code in range(1<<len(edges)):
            T=[0]*n
            for k,(a,b) in enumerate(edges):
                if code>>k&1:T[a]|=1<<b
                else:T[b]|=1<<a
            T=tuple(T);U=join(T,C);q=profile(T);v=profile(U);count+=1
            need(v['H']==3*q['H'] and v['P2']==12*q['P2'],'whole join readouts')
            for r in range(n):
                need(v['E1'][r]==0 and v['S1'][r]==3*q['S1'][r],'one-path endpoints')
                need(v['deleted'][r]==3*q['deleted'][r] and v['E2'][r]==3*q['E2'][r] and v['S2'][r]==12*q['S2'][r],'rooted join coefficients')
                for s in range(n):need(v['W2'][r][s]==3*q['W2'][r][s],'joint retained old ports')
            need(not strong(U),'terminal SCC remains separate')
    T=join(C,C);q=profile(T);need(q['H']==9 and brute_h(T)==9,'literal zero-defect weak H9 parent')
    hits=[];spectrum=[]
    for cut in range(64):
        U=birth(T,cut);v=profile(U);spectrum.append(v['H'])
        if v['H'] in (29,31):
            need(brute_h(U)==v['H'],'independent permutation count for absorption witness')
            hits.append(dict(cut=cut,rows=U,H=v['H'],strong=strong(U),profile=v))
    A=(20,1,2,7,14);need(profile(A)['H']==brute_h(A)==9 and strong(A),'actual alternate strong H9 source');other=[]
    for cut in range(32):
        U=birth(A,cut);v=profile(U)
        if v['H'] in (29,31):
            need(brute_h(U)==v['H'],'independent strong-source positive branch')
            other.append(dict(cut=cut,rows=U,H=v['H'],strong=strong(U)))
    return dict(join_controls=count,universe='every labelled T of order1..4, all retained old ports; append terminal directed C3',source_H9=dict(rows=T,profile=q),all_64_one_birth_counts=spectrum,one_birth_29_31_hits=hits,alternate_strong_H9=dict(rows=A,profile=profile(A),all32cuts_29_31_hits=other),scope='Exact finite SCC-absorption probes, not an all-source tripling/return constructor.')
def main():
    result=dict(status='PROVED_ARITHMETIC_CRITERIA_WITH_FINITE_EXACT_CONTROLS',arithmetic=arithmetic(),graphs=graph_controls())
    (HERE/'coverage.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: three exact seed grammars through10001, all unary170807 depths0..10, 75 exhaustive terminal-C3 port identities and 64 actual H9 birth cuts')
    print('absorption29/31 hits',[(z['H'],z['cut'],z['strong']) for z in result['graphs']['one_birth_29_31_hits']])
if __name__=='__main__':main()
