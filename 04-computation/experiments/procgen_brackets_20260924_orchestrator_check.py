#!/usr/bin/env python3
"""Orchestrator's independent audit of THM-4470 (pairing ladder) and the brackets tables, written without reading the lane's code.

Checks (each raises on failure):
  Q1  pair-sum identity to i = 10^6; sum_(n<=2M) T(n) = sum n; the (q, r, offset) census of sum-preserving pairings
      over odd q <= 39, |r| <= 39 is exactly {(3, 1, {2i-1,2i}), (3, -1, {2i,2i+1})}.
  Q2  single-flip fragile pairs: exactly the 24 listed for T and the 12 listed for 3n-1, searched to i = 20000.
  Q3  bracket-internal Collatz moves (both sheets) and the complete (k, p) list; all-integer (k, n) list.
  Q4  a square-sum chain of [1,25] (found by DFS) and the 49 alternating chains partitioning [25, 1249].
  Q5  falsifying chain c_(t+1) = c_t + ceil(c_t/2) from 3 to 10^30: 169 terms, 82 flipped pairs, the flipped map follows it.
Run: python3 04-computation/experiments/procgen_brackets_20260924_orchestrator_check.py   (about 1-2 minutes, < 300 MB)
"""
import math, sys
sys.setrecursionlimit(10000)
def T(n, q=3, r=1): return n // 2 if n % 2 == 0 else (q * n + r) // 2
# Q1
N = 10 ** 6
assert all(T(2*i-1) + T(2*i) == 4*i - 1 for i in range(1, N))
assert all(T(2*i+1, 3, -1) + T(2*i, 3, -1) == 4*i + 1 for i in range(0, N))
assert sum(T(n) for n in range(1, 2*N + 1)) == sum(range(1, 2*N + 1))
hits = []
for q in range(1, 40, 2):
    for r in range(-39, 40, 2):
        for off in (0, 1):
            if all(T(2*i-1+off, q, r) + T(2*i+off, q, r) == 4*i - 1 + 2*off for i in range(1, 200)):
                hits.append((q, r, off))
assert sorted(hits) == [(3, -1, 1), (3, 1, 0)], hits
print("Q1  pair sums preserved; only (3,+1) on {2i-1,2i} and (3,-1) on {2i,2i+1}: ok")
# ---------- 1. single-flip fragile pairs ----------
def make_F(sheet, flip_i):
    # sheet +: pairs {2i-1,2i}; length ceil(n/2); odd goes up.  sheet -: pairs {2i,2i+1}; length floor(n/2); odd goes up (3n-1)/2.
    def F(n):
        if sheet=='+':
            i=(n+1)//2; L=i; up=(n%2==1)
        else:
            i=n//2; L=i; up=(n%2==1)
        if i==flip_i: up=not up
        return n+L if up else n-L
    return F
def fragile(sheet, imax, maxsteps=100000):
    out=[]
    base=set([1,2]) if sheet=='+' else None
    for i in range(1,imax+1):
        F=make_F(sheet,i)
        starts=[2*i-1,2*i] if sheet=='+' else [2*i,2*i+1]
        newcyc=False
        for s in starts:
            x=s; seen={}
            for t in range(maxsteps):
                if x in seen: break
                seen[x]=t; x=F(x)
            else:
                raise SystemExit(f"no cycle found from {s} (i={i})")
            # x is on the cycle: collect it
            cyc=[x]; y=F(x)
            while y!=x: cyc.append(y); y=F(y)
            cs=set(cyc)
            if sheet=='+':
                trivial = cs<= {0,1,2} and (cs=={1,2} or cs=={0})
            else:
                known=[{0},{1},{5,7,10},{17,25,37,55,82,41,61,91,136,68,34}]
                trivial = any(cs==k for k in known)
            if not trivial: newcyc=True
        if newcyc: out.append(i)
    return out
fp=fragile("+",20000); assert fp==[1,4,5,10,11,13,20,22,40,61,84,122,126,167,189,217,244,325,334,433,445,577,1154,2308], fp; print("Q2  T fragile pairs (i<=20000):",fp)
fm=fragile("-",20000); assert fm==[1,2,20,30,41,43,45,86,205,365,410,12029], fm; print("Q2  3n-1 fragile pairs (i<=20000):",fm)
# ---------- 2. escape sets / microcosm ----------
def br(n): return (math.isqrt(n-1)+1)//2 if n>1 else 0   # bracket index m with n in ((2m-1)^2,(2m+1)^2]
def br_real(x):  # bracket of real x>=1
    m=0
    while (2*m+1)**2 < x: m+=1
    return m
assert all(((2*br(n)-1)**2 < n <= (2*br(n)+1)**2) for n in range(2,10**5))
odd_in=[n for n in range(1,10**5,2) if br((3*n+1)//2)==br(n)]
even_in=[n for n in range(2,10**5,2) if br(n//2)==br(n)]
minus_in=[n for n in range(1,10**5,2) if n>0 and (3*n-1)//2>=1 and br((3*n-1)//2)==br(n)]
print("(3n+1)/2 stays:",odd_in); print("n/2 stays:",even_in); print("(3n-1)/2 stays:",minus_in)
# (k,p) with kp in p's bracket
def primes(N):
    s=bytearray([1])*(N+1); s[0]=s[1]=0
    for i in range(2,int(N**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(N+1) if s[i]]
kp=[(k,p) for p in primes(10**5) for k in range(2,10) if br(k*p)==br(p)]
print("(k,p):",kp)
kn=[(k,n) for n in range(1,10**4) for k in range(2,10) if br(k*n)==br(n)]
print("(k,n) all integers:",kn)
# ---------- 3. square-sum: find a chain of 1..25, check the 49-chain partition of [25,1249] ----------
N=25; sq={s*s for s in range(2,80)}
adj={a:[b for b in range(1,N+1) if b!=a and a+b in sq] for a in range(1,N+1)}
def dfs(path,used):
    if len(path)==N: return path
    for b in adj[path[-1]]:
        if b not in used:
            used.add(b); r=dfs(path+[b],used)
            if r: return r
            used.discard(b)
    return None
chain=None
for s0 in range(1,N+1):
    chain=dfs([s0],{s0})
    if chain: break
assert chain and sorted(chain)==list(range(1,26)) and all(chain[k]+chain[k+1] in sq for k in range(24))
allvals=[]
for c in range(-24,25):
    ch=[49*a+(c if k%2==0 else -c) for k,a in enumerate(chain)]
    assert all(math.isqrt(ch[k]+ch[k+1])**2==ch[k]+ch[k+1] for k in range(24))
    allvals+=ch
assert sorted(allvals)==list(range(25,1250))
print("square-sum chain of [1,25]:",chain,"; 49 alternating chains partition [25,1249]: ok")
# ---------- 4. falsifying chain: flip pairs of even chain members; modified map climbs forever ----------
c=[3]
while c[-1]<10**30: c.append(c[-1]+(c[-1]+1)//2)
flipped={x//2 for x in c if x%2==0}   # pair index of even member 2i is i (sheet +)
oddpairs={(x+1)//2 for x in c if x%2==1}
assert not (flipped & oddpairs)        # no chain odd member shares a flipped pair
def Fmod(n):
    i=(n+1)//2; up=(n%2==1)
    if i in flipped: up=not up
    return n+i if up else n-i
x=3; ok=True
for t in range(len(c)-1):
    assert Fmod(c[t])==c[t+1]
print("falsifying chain: modified map follows",len(c),"chain terms to 1e30 with",len(flipped),"flips: ok")
assert odd_in == [3,5,11,13,15,27,29,31,51,53] and even_in == [4,6,8,20,22,24]
assert minus_in == [1,3,5,11,13,15,17,27,29,31,33,51,53]
assert kp == [(2,2),(3,2),(4,2),(2,3),(3,3),(2,11)]
assert kn == [(2,2),(3,2),(4,2),(2,3),(3,3),(2,4),(2,10),(2,11),(2,12)]
assert len(c) == 169 and len(flipped) == 82
print("Q3-Q5 asserted: ok")
print("ALL CHECKS PASSED")
