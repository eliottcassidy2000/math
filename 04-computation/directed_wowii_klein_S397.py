#!/usr/bin/env python3
"""
klein-2026-07-21-S397 -- DIRECTED WOWII: formulate and EXHAUSTIVELY TEST directed-tournament
analogs of Written-on-the-Wall-II invariant inequalities.  Report candidate theorems (hold on
all iso classes n<=7) and explicit counterexamples (fail, with the smallest witness).

Owner: "work on the directed analogies of the WOWII inequalities."

DIRECTED INVARIANT ZOO (undirected -> directed):
  alpha (independence, no edge)      -> tr(T)  = largest TRANSITIVE subtournament (no 3-cycle)
  chromatic chi                      -> dichromatic dichr(T) = min acyclic-colouring
  clique omega                       -> tr(T) again
  domination gamma                   -> gamma(T) = min dominating set
  (no analog)                        -> fas(T) = min feedback arc set; #C3; H(T) (Redei); diam
"""
import itertools
from math import log2, ceil, floor, comb, log

def pairs_of(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def relabel(om,perm,n):
    new=[0]*n
    for v in range(n):
        mv,t=om[v],0
        while mv:
            b=mv&-mv; w=b.bit_length()-1; mv^=b; t|=1<<perm[w]
        new[perm[v]]=t
    return tuple(new)
def word(om,n):
    w=0
    for v in range(n): w=(w<<n)|om[v]
    return w
def refine(om,n):
    col=[bin(om[v]).count("1") for v in range(n)]
    while True:
        sig=[]
        for v in range(n):
            cnt={}; mv=om[v]
            while mv:
                b=mv&-mv; w=b.bit_length()-1; mv^=b; cnt[col[w]]=cnt.get(col[w],0)+1
            sig.append((col[v],tuple(sorted(cnt.items()))))
        order=sorted(set(sig)); nc=[order.index(sig[v]) for v in range(n)]
        if nc==col: break
        col=nc
    cells={}
    for v in range(n): cells.setdefault(col[v],[]).append(v)
    return [tuple(cells[k]) for k in sorted(cells)]
def canon(om,n):
    cells=refine(om,n); base=[]; pos=0
    for c in cells: base.append((c,pos)); pos+=len(c)
    best=None
    for ch in itertools.product(*[itertools.permutations(c) for (c,_) in base]):
        perm=[0]*n
        for (blk,(c,st)) in zip(ch,base):
            for k,v in enumerate(blk): perm[v]=st+k
        w=word(relabel(om,perm,n),n)
        if best is None or w<best: best=w
    return best
def classes(n):
    P=pairs_of(n); om0=tuple(sum(1<<j for j in range(i)) for i in range(n))
    seen={canon(om0,n):om0}; fr=[om0]
    while fr:
        nx=[]
        for om in fr:
            for (i,j) in P:
                nm=list(om)
                if om[i]>>j&1: nm[i]&=~(1<<j); nm[j]|=1<<i
                else: nm[j]&=~(1<<i); nm[i]|=1<<j
                nm=tuple(nm); cc=canon(nm,n)
                if cc not in seen: seen[cc]=nm; nx.append(nm)
        fr=nx
    return list(seen.values())

def beats(om,a,b): return (om[a]>>b)&1
def is_transitive_subset(om, verts):
    # acyclic iff repeatedly a source (out-degree = size-1 within subset) exists
    S=list(verts)
    while S:
        src=None
        for v in S:
            if all(v==u or beats(om,v,u) for u in S): src=v; break
        if src is None: return False
        S.remove(src)
    return True
def tr(om,n):
    for size in range(n,0,-1):
        for sub in itertools.combinations(range(n),size):
            if is_transitive_subset(om,sub): return size
    return 1
def num_c3(om,n):
    c=0
    for a,b,cc in itertools.combinations(range(n),3):
        # 3-cycle iff not transitive triple
        wins=[0,0,0]; tri=[a,b,cc]
        for i in range(3):
            for j in range(3):
                if i!=j and beats(om,tri[i],tri[j]): wins[i]+=1
        if sorted(wins)==[1,1,1]: c+=1
    return c
def gamma(om,n):
    for size in range(1,n+1):
        for S in itertools.combinations(range(n),size):
            Sset=set(S)
            if all(v in Sset or any(beats(om,s,v) for s in S) for v in range(n)): return size
    return n
def fas(om,n):
    best=comb(n,2)
    for perm in itertools.permutations(range(n)):
        back=0
        pos={v:i for i,v in enumerate(perm)}
        for a in range(n):
            for b in range(n):
                if a!=b and beats(om,a,b) and pos[a]>pos[b]: back+=1
        best=min(best,back)
        if best==0: break
    return best
def dichromatic(om,n):
    # min k: partition vertices into k acyclic (transitive) sets
    if is_transitive_subset(om,range(n)): return 1
    for k in (2,3,4):
        # try all colourings into k classes, each acyclic
        for col in itertools.product(range(k),repeat=n):
            ok=True
            for c in range(k):
                cls=[v for v in range(n) if col[v]==c]
                if not is_transitive_subset(om,cls): ok=False; break
            if ok: return k
    return 5
def strong(om,n):
    # strongly connected?
    for s in range(n):
        seen={s}; st=[s]
        while st:
            u=st.pop()
            for w in range(n):
                if beats(om,u,w) and w not in seen: seen.add(w); st.append(w)
        if len(seen)!=n: return False
    return True
def diam(om,n):
    if not strong(om,n): return None
    d=0
    for s in range(n):
        dist=[-1]*n; dist[s]=0; q=[s]
        while q:
            nq=[]
            for u in q:
                for w in range(n):
                    if beats(om,u,w) and dist[w]<0: dist[w]=dist[u]+1; nq.append(w)
            q=nq
        d=max(d,max(dist))
    return d
def H(om,n):
    c=0
    def go(l,u,dep):
        nonlocal c
        if dep==n: c+=1; return
        for w in range(n):
            if beats(om,l,w) and not (u>>w&1): go(w,u|(1<<w),dep+1)
    for s in range(n): go(s,1<<s,1)
    return c

# ---- collect invariants ----
data={}
for n in (4,5,6,7):
    rows=[]
    for om in classes(n):
        rows.append(dict(n=n, tr=tr(om,n), c3=num_c3(om,n), gam=gamma(om,n),
                         fas=fas(om,n), dic=dichromatic(om,n), H=H(om,n),
                         diam=diam(om,n), strong=strong(om,n), om=om))
    data[n]=rows
    print(f"n={n}: {len(rows)} iso classes; tr range {min(r['tr'] for r in rows)}..{max(r['tr'] for r in rows)},"
          f" gamma range {min(r['gam'] for r in rows)}..{max(r['gam'] for r in rows)},"
          f" dichr max {max(r['dic'] for r in rows)}")

# ---- candidate directed WOWII inequalities ----
INEQ = [
 ("A: tr >= n - fas (reverse min-FAS -> transitive)", lambda r: r['tr'] >= r['n']-r['fas']),
 ("B: tr * dichr >= n (transitive colour classes)",   lambda r: r['tr']*r['dic'] >= r['n']),
 ("C: dichr <= ceil(n/tr)",                            lambda r: r['dic'] <= ceil(r['n']/r['tr'])),
 ("D: gamma <= floor(log2 n)+1 (domination bound)",    lambda r: r['gam'] <= floor(log2(r['n']))+1),
 ("E: tr >= floor(log2 n)+1 (Erdos-Moser)",            lambda r: r['tr'] >= floor(log2(r['n']))+1),
 ("F: tr + gamma <= n+1",                              lambda r: r['tr']+r['gam'] <= r['n']+1),
 ("G: DIRECTED-103  tr <= floor(n - log(diam)) [strong]", lambda r: (r['diam'] is None) or (r['tr'] <= floor(r['n']-log(r['diam'])))),
 ("H: c3 <= comb(n,3) - comb(tr,3)",                   lambda r: r['c3'] <= comb(r['n'],3)-comb(r['tr'],3)),
 ("I: gamma <= n - tr + 1",                            lambda r: r['gam'] <= r['n']-r['tr']+1),
 ("J: H <= 2^{n - tr} (transitivity caps path count)", lambda r: r['H'] <= 2**(r['n']-r['tr'])),
]
print("\n"+"="*84)
print("DIRECTED WOWII INEQUALITIES -- exhaustive over all iso classes n=4..7")
print("="*84)
for name,f in INEQ:
    fails=[]
    for n in (4,5,6,7):
        for r in data[n]:
            if not f(r): fails.append(r); break
        if fails: break
    if not fails:
        print(f"  HOLDS  {name}")
    else:
        r=fails[0]
        print(f"  FAILS  {name}")
        print(f"         smallest counterexample n={r['n']}: tr={r['tr']} gamma={r['gam']} fas={r['fas']}"
              f" dichr={r['dic']} c3={r['c3']} H={r['H']} diam={r['diam']}")
