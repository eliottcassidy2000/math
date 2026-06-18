"""Fast decisive check: the M3b map depends ONLY on residues mod 14.
So the COMPLETE set of n=5 tournaments producible by M3b = run the map over
ALL residue-vectors mod 14 (with all reference choices and both tie-break orders).
If a forbidden class is not in that complete set, it is provably unreachable for n=5
(as a top-level tournament). We also report which of the 12 n=5 iso classes ARE reachable.
"""
from itertools import combinations, permutations, combinations_with_replacement
from math import gcd

UNITS=[1,3,5,9,11,13]
def d0(v,a):
    s=(v*a)%14
    return min(s,14-s)
def w_ref(v0):
    return {a:(1 if (v0*a)%14%2==0 else -1) for a in UNITS}

def build(res, ref_res, tie_order):
    """res = tuple of 5 residues mod 14 (the vertices' v mod 14).
    ref_res = residue of reference runner. tie_order = permutation giving speed order
    for tie-break (larger speed wins). We model 'speed order' abstractly by tie_order:
    tie_order is a permutation of indices 0..4 giving relative size (last = largest)."""
    n=len(res)
    w=w_ref(ref_res)
    # rank[i] = relative speed rank; larger rank = larger speed
    rank=[0]*n
    for pos,idx in enumerate(tie_order):
        rank[idx]=pos
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j: continue
            vote=0
            for a in UNITS:
                dx=d0(res[i],a); dy=d0(res[j],a)
                if dx>dy: vote+=w[a]
                elif dx<dy: vote-=w[a]
            if vote>0: adj[i][j]=1
            elif vote<0: adj[i][j]=0
            else:
                adj[i][j]=1 if rank[i]>rank[j] else 0
    for i in range(n):
        for j in range(i+1,n):
            assert adj[i][j]+adj[j][i]==1
    return adj

def scoreseq(adj):
    n=len(adj); return tuple(sorted(sum(adj[i]) for i in range(n)))
def c3(adj):
    n=len(adj); c=0
    for i in range(n):
        for j in range(n):
            if i==j or not adj[i][j]: continue
            for k in range(n):
                if k in (i,j): continue
                if adj[j][k] and adj[k][i]: c+=1
    return c//3
def H(adj):
    n=len(adj); full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            cur=dp[mask][last]
            if not cur: continue
            for nxt in range(n):
                if mask&(1<<nxt): continue
                if adj[last][nxt]: dp[mask|(1<<nxt)][nxt]+=cur
    return sum(dp[full][v] for v in range(n))
def sk(adj): return (H(adj),c3(adj),scoreseq(adj))

FORBIDDEN={(9,3,(1,1,2,3,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))}

def main():
    # complete residue exhaust for n=5 top-level
    produced=set()
    hits={f:[] for f in FORBIDDEN}
    # residues: any 5 from 0..13 (with replacement allowed). For loneliness usually
    # we want section!=0, but residue 0 (=mult of 14) is allowed (the parked runner).
    res_vecs=list(combinations_with_replacement(range(0,14),5))
    print(f"residue vectors (multisets of 5 from 0..13): {len(res_vecs)}")
    # to fully cover tie-break, try all reference residues present, and several tie orders.
    # tie-break only matters when vote ties; trying all 120 permutations is overkill but safe for a sample;
    # we use a sufficient set: identity and reverse and a few rotations -> but to be airtight, try ALL perms.
    perms=list(permutations(range(5)))
    cnt=0
    for res in res_vecs:
        refset=set(res)
        for ref in refset:
            for to in perms:
                adj=build(res,ref,to)
                k=sk(adj)
                produced.add(k)
                if k in FORBIDDEN:
                    hits[k].append((res,ref,to))
        cnt+=1
    print(f"distinct n=5 short-keys produced by M3b (residue-complete): {len(produced)}")
    # the 12 actual n=5 iso classes have short-keys:
    print("\nAll produced short-keys (H,c3,score):")
    for k in sorted(produced): print("  ",k)
    print("\nForbidden-target check:")
    for f in sorted(FORBIDDEN):
        if hits[f]:
            print(f"  REFUTED {f}: example {hits[f][0]}  (total {len(hits[f])})")
        else:
            print(f"  NOT produced (forbidden holds at residue-complete level): {f}")

if __name__=="__main__":
    main()
