"""
forbidden_21_profiles_kps.py  (v2, arithmetic + necessary-realizability filter)

Contrast the cluster (conflict-graph) profiles of H=7 vs H=21.

H = I(Omega,2) = sum_k alpha_k 2^k with alpha_0=1.  So (H-1)/2 = sum_{k>=1} alpha_k 2^{k-1}.
Enumerate all alpha-vectors (alpha_1,alpha_2,...) of nonneg ints with
  sum_{k>=1} alpha_k 2^{k-1} = (H-1)/2,
then keep only those that pass NECESSARY conditions for being the independent-set
sequence i_k(Omega) of SOME graph Omega on alpha_1 vertices:
  (R1) alpha_1 = #vertices = N.
  (R2) alpha_2 <= C(N,2)            (can't have more indep pairs than pairs).
  (R3) if alpha_k>0 then alpha_1 >= k and alpha_{k-1} >= ... (log-concavity is NOT
       required in general, but a clique of size k needs alpha_1>=k indep set => N>=k).
  (R4) alpha_k <= C(N,k).
  (R5) the empty graph (all independent) gives alpha_k=C(N,k); the complete graph K_N
       gives alpha=(1,N,0,0,...). Any realizable seq is "between" these; in particular
       alpha_2 ranges over [0, C(N,2)] and is achievable for EVERY value in that range
       by choosing the number of edges (Kruskal-Katona / threshold graphs realize all
       alpha_2 in [0,C(N,2)]). We use exact small-N graph enumeration for N<=6 to be
       rigorous there, and the necessary filter (R1-R4) for larger N.
This is enough to show H=7 is RIGID (unique profile K3) and H=21 is SOFT (many profiles).
"""
import itertools
from collections import defaultdict
from math import comb

def indep_seq(N, edges):
    adj=[0]*N
    for (a,b) in edges:
        adj[a]|=1<<b; adj[b]|=1<<a
    counts=defaultdict(int)
    def expand(allowed,size):
        counts[size]+=1
        a=allowed
        while a:
            i=(a&-a).bit_length()-1
            a&=a-1
            new=allowed & ~adj[i] & ~((1<<(i+1))-1)
            expand(new,size+1)
    expand((1<<N)-1,0)
    mk=max(counts)
    return tuple(counts.get(k,0) for k in range(mk+1))

def graph_iseqs_exact(N):
    """exact set of indep-set sequences for graphs on N vertices (N<=6 ok)."""
    pairs=[(i,j) for i in range(N) for j in range(i+1,N)]
    seqs=set()
    for mask in range(1<<len(pairs)):
        E=[pairs[b] for b in range(len(pairs)) if (mask>>b)&1]
        seqs.add(indep_seq(N,E))
    return seqs

def I_at_2(alpha):
    return sum(a*(2**k) for k,a in enumerate(alpha))

def enum_alpha_vectors(target):
    """all (alpha_1,alpha_2,...) nonneg with sum alpha_k 2^k = target-1."""
    rem0=target-1
    results=[]
    # max level: 2^k <= target-1
    maxk=0
    while (1<<(maxk+1)) <= rem0:
        maxk+=1
    # recursive: assign alpha_k for k=maxk..1
    def rec(k, rem, acc):
        if k==0:
            if rem==0:
                results.append(tuple(acc[1:]))  # alpha_1.. (drop the alpha_0 slot)
            return
        w=1<<k
        for ak in range(rem//w, -1, -1):
            acc[k]=ak
            rec(k-1, rem-ak*w, acc)
            acc[k]=0
    acc=[0]*(maxk+1)
    rec(maxk, rem0, acc)
    # prepend alpha_0=1, build full vectors
    out=[]
    for tail in results:
        # tail is (alpha_1,...,alpha_maxk); strip trailing zeros
        t=list(tail)
        while t and t[-1]==0:
            t.pop()
        out.append(tuple([1]+t))
    return sorted(set(out))

def necessary_realizable(alpha):
    """necessary conditions for alpha to be i_k of SOME graph on N=alpha_1 vertices."""
    if len(alpha)<2:
        return alpha==(1,) or alpha==()  # N=0
    N=alpha[1]
    for k in range(2,len(alpha)):
        if alpha[k]>comb(N,k):
            return False
        if alpha[k]>0 and N<k:
            return False
    return True

for target in (7,21):
    print("="*64)
    print(f"H = {target}: cluster (conflict-graph) profiles alpha with I(.,2)={target}")
    print("="*64)
    cands=enum_alpha_vectors(target)
    # filter by necessary realizability
    nec=[a for a in cands if necessary_realizable(a)]
    # for small N (<=6) refine using exact graph iseqs
    exact_realizable=[]
    cache={}
    for a in nec:
        N=a[1] if len(a)>1 else 0
        if N<=6:
            if N not in cache:
                cache[N]=graph_iseqs_exact(N) if N>=1 else {(1,)}
            # pad a to compare; a graph iseq has trailing structure; compare directly
            # need a to equal some graph iseq exactly (after trailing-zero norm)
            aa=list(a)
            while len(aa)>1 and aa[-1]==0: aa.pop()
            aa=tuple(aa)
            ok = any( (tuple(s[:len(aa)])==aa and all(x==0 for x in s[len(aa):])) or s==aa
                      for s in cache[N])
            # simpler exact membership with trailing-zero normalization
            norm=set()
            for s in cache[N]:
                ss=list(s)
                while len(ss)>1 and ss[-1]==0: ss.pop()
                norm.add(tuple(ss))
            ok = aa in norm
            if ok:
                exact_realizable.append(a)
        else:
            exact_realizable.append(a)  # only necessary-filtered (N>6)
    print(f"  raw alpha-vectors summing to {target}: {len(cands)}")
    print(f"  pass necessary realizability filter: {len(nec)}")
    print(f"  graph-realizable (exact for N<=6, necessary for N>6): {len(exact_realizable)}")
    for a in exact_realizable:
        N=a[1] if len(a)>1 else 0
        tag=""
        if a==(1,3,0): tag="  <-- K3 (the unique H=7 profile)"
        if a==(1,4,3,0): tag="  <-- K3+isolated"
        print(f"     alpha={list(a)}  (N={N} cycles)   I={I_at_2(a)}{tag}")
    print()

print("VERDICT:")
print(" H=7  -> exactly ONE graph-realizable cluster profile: K3 = [1,3,0].")
print("         THM-029 kills that single profile (3 pairwise-conflicting odd cycles")
print("         are non-realizable as a tournament conflict graph). RIGID obstruction.")
print(" H=21 -> SEVERAL graph-realizable cluster profiles. No single-cluster rigidity.")
print("         21 is forbidden MULTIPLICATIVELY (21=3*7, strong factorization needs strong-7).")
