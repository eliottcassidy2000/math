#!/usr/bin/env python3
"""Literal checks for the delta-three Hamilton-path orbit spectrum."""
import argparse, collections, functools, itertools, json, math, pathlib

checks=0

def check(ok, message):
    global checks
    checks+=1
    if not ok: raise RuntimeError(message)

def words(maxweight):
    yield ()
    for first in (1,3):
        if maxweight>=first:
            for rest in words(maxweight-first): yield (first,)+rest

def tournament(word):
    n=2+sum(word); out=[0]*n; blocks=[]; k=1
    for size in word:
        block=tuple(range(k,k+size)); blocks.append(block); k+=size
    for u in range(n):
        for v in range(u+1,n): out[u]|=1<<v
    out[0]&=~(1<<(n-1)); out[n-1]|=1
    for block in blocks:
        if len(block)==3:
            u,_,v=block; out[u]&=~(1<<v); out[v]|=1<<u
    return out,blocks

def hp_dp(out):
    n=len(out); full=(1<<n)-1
    @functools.lru_cache(None)
    def f(mask,last):
        if mask==full:return 1
        return sum(f(mask|1<<v,v) for v in range(n) if not mask>>v&1 and out[last]>>v&1)
    return sum(f(1<<v,v) for v in range(n))

def paths(out):
    n=len(out); full=(1<<n)-1
    def walk(prefix,mask):
        if mask==full:yield prefix;return
        avail=out[prefix[-1]]&~mask
        while avail:
            bit=avail&-avail;avail-=bit;v=bit.bit_length()-1
            yield from walk(prefix+(v,),mask|bit)
    for v in range(n): yield from walk((v,),1<<v)

def isomorphisms(out,anti=False):
    n=len(out);deg=[x.bit_count() for x in out]
    cand=[[v for v in range(n) if deg[v]==(n-1-deg[u] if anti else deg[u])] for u in range(n)]
    source=sorted(range(n), key=lambda u:len(cand[u]));p=[-1]*n
    def rec(k,used):
        if k==n:yield tuple(p);return
        u=source[k]
        for v in cand[u]:
            if used>>v&1:continue
            if all(((out[u]>>s)&1)==(((out[v]>>p[s])&1)^anti) for s in source[:k]):
                p[u]=v;yield from rec(k+1,used|1<<v);p[u]=-1
    return list(rec(0,0))

def feedback(out):
    n=len(out);f=[0]*(1<<n)
    for mask in range(1,1<<n):
        f[mask]=min(f[mask^(1<<v)]+(out[v]&(mask^(1<<v))).bit_count() for v in range(n) if mask>>v&1)
    return f[-1]

def canonical(path,autos):
    return min(tuple(p[v] for v in path) for p in autos)

def poly_e(e):
    return {e+t+2:math.comb(e-t,t) for t in range(e//2+1)}

def main():
    rows=[];total_paths=0
    for word in words(10):
        if sum(word)<2:continue
        out,blocks=tournament(word);n=len(out);a=word.count(1);t=word.count(3);e=a+2*t
        autos=isomorphisms(out);antis=isomorphisms(out,True)
        H=hp_dp(out); ps=list(paths(out));total_paths+=len(ps)
        check(H==len(ps),f'DP/DFS disagreement {word}')
        check(H==3**t+2**a*12**t,f'raw H {word}')
        check(len(autos)==3**t,f'actual aut {word}')
        support=sum(any(p[v]!=v for p in autos) for v in range(n))
        check(support==3*t,f'support {word}')
        check(3*feedback(out)-support==3,f'actual slack {word}')
        reps={canonical(p,autos) for p in ps}
        check(len(reps)==1+2**e,f'path orbit H/aut {word}')
        # Unique direct orbit, then independent letter occupation numbers.
        direct=set();decoded={}
        for p in reps:
            if p.index(n-1)+1==p.index(0):
                left=set(p[:p.index(n-1)])
                js=tuple(len(left.intersection(block)) for block in blocks)
                check(js not in decoded,f'collision in split orbit encoding {word},{js}')
                decoded[js]=p
            else: direct.add(p)
        check(len(direct)==1,f'direct orbit {word}')
        check(set(decoded)==set(itertools.product(*(range(s+1) for s in word))),f'complete local alphabets {word}')
        check(bool(antis)==(word==word[::-1]),f'converse criterion {word}')
        b=0
        if antis:
            anti=antis[0]
            for p in reps:
                if canonical(tuple(anti[v] for v in p[::-1]),autos)==p:b+=1
            expected=1 if n%2 else 1+2**(e//2)
            check(b==expected,f'blue count {word}: {b}!={expected}')
            check(all(canonical(tuple(anti[v] for v in p[::-1]),autos) in direct for p in direct),f'direct blue {word}')
            for js,p in decoded.items():
                rev=tuple(s-j for s,j in zip(word,js[::-1]))
                transformed=canonical(tuple(anti[v] for v in p[::-1]),autos)
                check(transformed==decoded[rev],f'blue action alphabet {word},{js}')
        rows.append({'word':list(word),'n':n,'e':e,'H':H,'aut':len(autos),'h':len(reps),'self_converse':bool(antis),'b':b})
    # Independent finite word coefficients; weights are deliberately separate.
    gf=[]
    for e in range(2,15):
        count=collections.Counter();pal=collections.Counter();even=collections.Counter();odd=collections.Counter()
        def ew(rem,prefix=()):
            if rem==0:yield prefix;return
            if rem>=1:yield from ew(rem-1,prefix+(1,))
            if rem>=2:yield from ew(rem-2,prefix+(3,))
        for w in ew(e):
            n=2+sum(w);count[n]+=1
            if w==w[::-1]:
                pal[n]+=1;(even if n%2==0 else odd)[n]+=1
        check(dict(count)==poly_e(e),f'fixed e polynomial {e}')
        check(sum(count.values())==sum(math.comb(e-t,t) for t in range(e//2+1)),f'Fibonacci total {e}')
        gf.append({'e':e,'strong_coefficients':dict(sorted(count.items())),'self_converse_even':dict(sorted(even.items())),'self_converse_odd':dict(sorted(odd.items()))})
    check(poly_e(2)=={4:1,5:1},'binom(e,t) hostile e=2')
    check(next(r for r in rows if r['word']==[1,1])['b']==3,'S4 blue control')
    check(next(r for r in rows if r['word']==[3])['b']==1,'nonregular H15 control')
    # Symbolic series identities, including the two parity channels and padding.
    import sympy as sp
    z,u=sp.symbols('z u');D=1-u*z-u*u*z**3;P=1-z-z**3;Q=1-u*u*z*z-u**4*z**6
    S=z*z*(1/D-1-u*z);E=z*z*(1/Q-1);O=z*z*((u*z+u*u*z**3)/Q-u*z)
    for item in gf:
        e=item['e']
        check(sp.expand(sp.series(S,u,0,e+1).removeO()).coeff(u,e)==sum(c*z**n for n,c in item['strong_coefficients'].items()),f'symbolic ordinary {e}')
        check(sp.expand(sp.series(E,u,0,e+1).removeO()).coeff(u,e)==sum(c*z**n for n,c in item['self_converse_even'].items()),f'symbolic even {e}')
        check(sp.expand(sp.series(O,u,0,e+1).removeO()).coeff(u,e)==sum(c*z**n for n,c in item['self_converse_odd'].items()),f'symbolic odd {e}')
    check(sp.simplify((S/P**2).subs(u,1)-z*z*(1/P-1-z)/P**2)==0,'inherited slack-three GF recovery')
    q=1-z*z-z**6
    check(sp.Poly(z**3+z+1,z,modulus=2).is_irreducible,'mod2 cubic irreducible')
    check(sp.Poly(z*z+1,z,modulus=3).is_irreducible,'mod3 quadratic irreducible')
    check(sp.Poly(z**4-z*z-1,z,modulus=3).is_irreducible,'mod3 quartic irreducible')
    check(sp.gcd(P,q)==1,'coprime p and q')
    for item in gf:
        Pe=sum(c*z**n for n,c in item['strong_coefficients'].items())
        Ee=sum(c*z**n for n,c in item['self_converse_even'].items())
        Oe=sum(c*z**n for n,c in item['self_converse_odd'].items())
        for numerator in (Ee,Oe):
            if numerator!=0:check(sp.gcd(numerator,q)==1,'exact selfconverse denominator')
        check(sp.gcd(Pe*q-(Ee+Oe)*P**2,P**2*q)==1,'exact nonselfconverse denominator')
    # Two literal class squares: midpoint cut, or one straddling triangle.
    @functools.lru_cache(None)
    def tilings(e):
        if e==0:return ((),)
        return tuple((s,)+rest for s in (1,2) if s<=e for rest in tilings(e-s))
    def atomize(w):return tuple(1 if s==1 else 3 for s in w)
    def bridge(A,B):
        level=sum(A)
        while level:
            if B[0]==1:return (level,atomize(A)+atomize(B[1:])[::-1])
            A,B=B[1:],A;level-=1
        return (0,())
    cumulative=[]
    for r in range(1,8):
        balanced={atomize(U)+atomize(V)[::-1] for U in tilings(r) for V in tilings(r)}
        straddle={atomize(U)+(3,)+atomize(V)[::-1] for U in tilings(r-1) for V in tilings(r-1)}
        whole={atomize(w) for w in tilings(2*r)}
        check(not balanced.intersection(straddle),'cut and straddling disjoint')
        check(balanced|straddle==whole,'cut or straddling exhaustive')
        check(len(balanced)==len(tilings(r))**2,'balanced class square injective')
        check(len(straddle)==len(tilings(r-1))**2,'straddling class square injective')
        for U in tilings(r):
            for V in tilings(r):
                W=atomize(U)+atomize(V)[::-1]
                check((W==W[::-1])==(U==V),'balanced class converse diagonal')
        rr=sum(z**sum(atomize(w)) for w in tilings(r))
        rm=sum(z**sum(atomize(w)) for w in tilings(r-1))
        check(sp.expand(z*z*rr**2+z**5*rm**2)==sum(c*z**n for n,c in poly_e(2*r).items()),'weighted midpoint identity')
    for r in range(0,11):
        image=[bridge(A,B) for A in tilings(r) for B in tilings(r+1)]
        target={(0,())}
        for e in range(1,r+1):
            target.update((e,atomize(U)+atomize(V)[::-1]) for U in tilings(e) for V in tilings(e))
        check(len(image)==len(set(image)),'old tiling pair bridge injective')
        check(set(image)==target,'old tiling pair bridge onto')
        cumulative.append(len(target))
    check(cumulative[:7]==[1,2,6,15,40,104,273],'inherited cumulative sequence')
    check((1,3,1) not in {atomize(U)+atomize(V)[::-1] for U in tilings(2) for V in tilings(2)},'131 missing midpoint hostile')
    data={'checks':checks,'literal_word_count':len(rows),'literal_paths':total_paths,'max_core_order':12,'rows':rows,'fixed_e':gf,'cumulative_balanced_cut_counts':cumulative,'status':'PASS'}
    parser=argparse.ArgumentParser();parser.add_argument('--receipt',default='receipt.json');args=parser.parse_args()
    pathlib.Path(args.receipt).write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ('rows','fixed_e')}))

if __name__=='__main__':main()
