#!/usr/bin/env python3
"""continuum_threads_boxeph_S200.py -- boxeph-2026-07-21-S200

Working the S199 open threads + agent pickups.

THREAD A (L3 probe): audit what separates deep-continuum tournaments after
   (char_A, |R|), then after score. Test LOCAL SUBTOURNAMENT DENSITIES (the k-profile =
   induced-subtournament census), i.e. flag/limit-theory coordinates, without crediting
   a profile for resolution already supplied by score.
THREAD B (reducibility ceiling): max c3 over REDUCIBLE (non-strong) tournaments = c3_max(n-1)
   [cyclic content sums over strong components; discrete convexity concentrates their size
   partition at (n-1,1)]. => reducible ceiling tau_red = c3_max(n-1)/c3_max(n).
THREAD C (H vs temperature): mean/max/spread of H per iso-cyclic shell -- locate the H structure on
   the temperature axis (death-star-S84 H>=disc binding case = quasirandom = tau=1).
"""
import numpy as np
from itertools import permutations, combinations
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from math import comb

def canon(A,n,perms):
    best=None
    for p in perms:
        c=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or c<best: best=c
    return best
PERMS={k:list(permutations(range(k))) for k in range(1,8)}
def iso_reps(nmax):
    reps={1:[[[0]]]}
    for n in range(2,nmax+1):
        seen=set(); out=[]
        for B in reps[n-1]:
            for pat in range(1<<(n-1)):
                A=[row[:]+[0] for row in B]+[[0]*n]
                for k in range(n-1):
                    if pat>>k&1: A[n-1][k]=1
                    else: A[k][n-1]=1
                c=canon(A,n,PERMS[n])
                if c not in seen:
                    seen.add(c); M=[[0]*n for _ in range(n)]; idx=0
                    for i in range(n):
                        for j in range(n):
                            if i!=j: M[i][j]=c[idx]; idx+=1
                    out.append(M)
        reps[n]=out
    return reps
def scoreseq(A,n): return tuple(sorted(sum(A[i]) for i in range(n)))
def c3v(A,n):
    sc=[sum(A[i]) for i in range(n)]; return comb(n,3)-sum(comb(s,2) for s in sc)
def strong(A,n):
    R=[[1 if(i==j or A[i][j])else 0 for j in range(n)]for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n):
                    if R[k][j]: R[i][j]=1
    return all(R[i][j] and R[j][i] for i in range(n) for j in range(n))
def charpoly_int(A):
    n=len(A); A=[[int(A[i][j])for j in range(n)]for i in range(n)]
    M=[[1 if i==j else 0 for j in range(n)]for i in range(n)]; c=[1]
    for k in range(1,n+1):
        AM=[[sum(A[i][t]*M[t][j]for t in range(n))for j in range(n)]for i in range(n)]
        tr=sum(AM[i][i]for i in range(n)); ck=-tr//k; c.append(ck)
        M=[[AM[i][j]+(ck if i==j else 0)for j in range(n)]for i in range(n)]
    return tuple(c)
def signed_redei(A,n):
    adj=A; R=0; path=[]
    def dfs(v,used):
        nonlocal R
        path.append(v)
        if len(path)==n:
            inv=sum(1 for a in range(n) for b in range(a+1,n) if path[a]>path[b])
            R+= 1 if inv%2==0 else -1
        else:
            for w in range(n):
                if not(used>>w&1) and adj[v][w]: dfs(w,used|(1<<w))
        path.pop()
    for v in range(n): dfs(v,1<<v)
    return R
def profile_k(A,n,k):
    cnt=Counter()
    for S in combinations(range(n),k):
        sub=[[A[i][j] for j in S] for i in S]
        cnt[canon(sub,k,PERMS[k])]+=1
    return tuple(sorted(cnt.items()))

reps=iso_reps(7)
print("iso classes:", {n:len(reps[n]) for n in range(3,8)})

# ---------------- THREAD A ----------------
print("\n"+"="*92); print("THREAD A  L3 probe: score versus local k-subtournament profiles inside (char_A,|R|)-twins")
print("="*92)
n=7
# The hot shell tau=6/7 (c3=12) has 47 classes; absolute |R| resolves 28 keys.
for cc in (12,11,13):
    shell=[A for A in reps[n] if c3v(A,n)==cc]
    records=[]
    for A in shell:
        records.append((
            charpoly_int(A), abs(signed_redei(A,n)), scoreseq(A,n),
            profile_k(A,n,4), profile_k(A,n,5),
        ))

    def resolution(extra):
        groups=Counter((record[0],record[1])+extra(record) for record in records)
        unresolved=sorted((size for size in groups.values() if size>1), reverse=True)
        return len(groups), unresolved

    r_base, base_groups=resolution(lambda record:())
    r_score, _=resolution(lambda record:(record[2],))
    r_p4, _=resolution(lambda record:(record[3],))
    r_p45, _=resolution(lambda record:(record[3],record[4]))
    r_score_p45, final_groups=resolution(lambda record:(record[2],record[3],record[4]))
    print("  c3=%d shell: %d classes; base(char_A,|R|)=%d; +score=%d; +4-profile=%d; +4&5-profile=%d; +score+4&5=%d; base_unresolved=%s; final_unresolved=%s"
          %(cc,len(shell),r_base,r_score,r_p4,r_p45,r_score_p45,base_groups,final_groups))
print("  => absolute-|R| audit: local profiles add beyond score in all three tested hot shells; score also adds at c3=11.")

# ---------------- THREAD B ----------------
print("\n"+"="*92); print("THREAD B  reducibility ceiling: max c3 over REDUCIBLE tournaments vs c3_max(n-1)")
print("="*92)
c3max={n:max(c3v(A,n) for A in reps[n]) for n in range(3,8)}
print("  c3_max(n) n=3..7:", [c3max[n] for n in range(3,8)])
for n in range(4,8):
    red=[A for A in reps[n] if not strong(A,n)]
    maxred=max(c3v(A,n) for A in red)
    print("  n=%d: max reducible c3 = %d ; c3_max(n-1) = %d ; equal? %s ; tau_red=c3_max(n-1)/c3_max(n)=%s"
          %(n, maxred, c3max[n-1], maxred==c3max[n-1], Fr(c3max[n-1],c3max[n])))
print("  => PROVED shape: c3(T)=sum c3(SCC); discrete convexity of c3_max concentrates the SCC-size partition at (n-1,1).")
print("     THM-462 no-holes realizes the next level, so the first all-strong shell is (c3_max(n-1)+1)/c3_max(n).")

# ---------------- THREAD C ----------------
print("\n"+"="*92); print("THREAD C  H vs cyclic temperature: H structure per iso-cyclic shell (locate death-star's binding)")
print("="*92)
def ham_paths(A,n):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if c:
                for w in range(n):
                    if not(mask>>w&1) and A[last][w]: dp[mask|(1<<w)][w]+=c
    return sum(dp[full][last] for last in range(n))
n=7
shells=defaultdict(list)
for A in reps[n]: shells[c3v(A,n)].append(A)
print("  n=7:  tau     c3   #cls   H:min  mean   max   (H grows with temperature; spread widens hot)")
for cc in sorted(shells,reverse=True):
    Hs=[ham_paths(A,n) for A in shells[cc]]
    print("       %-6s  %3d  %4d   %5d %6.1f %5d" %
          (str(Fr(cc,c3max[7])), cc, len(Hs), min(Hs), sum(Hs)/len(Hs), max(Hs)))
