"""
DECISIVE residue-complete enumeration for the M3b section-parity character vote, n=5.

KEY OBSERVATION (proved below by construction): the M3b map's arc direction x->y
depends ONLY on:
  (1) the residues v mod 14 of x and y (via d0(v,a)=min(v*a mod14, 14-v*a mod14)),
  (2) the reference runner's residue v0 mod 14 (via the parity weight w(a)),
  (3) ON TIE (vote sum == 0): the relative SPEED order of x and y (tie-break x>y).

Therefore the COMPLETE set of n=5 top-level tournaments producible by M3b over ALL
positive-integer speed sets = enumerate all residue multisets (5 residues from 0..13),
all reference residues, and -- to cover the tie-break fully -- all relative orderings
of the tie-involved vertices. We handle tie-break exactly by enumerating, for each
(residue-vector, ref), every total order of the vertices (a linear extension that the
actual speeds could induce). Since only tie pairs depend on the order, we collapse to
the set of tie pairs and try every acyclic orientation consistent with SOME speed order.

To stay airtight AND fast: for a given (res,ref), compute the vote matrix. Pairs with
vote!=0 are FIXED. Pairs with vote==0 form a "tie graph"; their orientation is given by
speed order, which is a STRICT TOTAL ORDER on the 5 vertices. So the tie edges must be
oriented consistently with a total order -> we enumerate all 5! relative orders but only
re-derive the tie edges (cheap). We collect every resulting tournament's iso short-key.

This is EXHAUSTIVE for n=5 top-level. No sampling.
"""
from itertools import combinations_with_replacement, permutations
import sys

UNITS=[1,3,5,9,11,13]
def d0(v,a):
    s=(v*a)%14
    return min(s,14-s)

def vote_matrix(res,ref):
    n=len(res)
    w={a:(1 if (ref*a)%14%2==0 else -1) for a in UNITS}
    # precompute d0 per vertex
    D=[[d0(res[i],a) for a in UNITS] for i in range(n)]
    V=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            vote=0
            for ai,a in enumerate(UNITS):
                dx=D[i][ai]; dy=D[j][ai]
                if dx>dy: vote+=w[a]
                elif dx<dy: vote-=w[a]
            V[i][j]=vote; V[j][i]=-vote
    return V

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
        row=dp[mask]
        for last in range(n):
            cur=row[last]
            if not cur: continue
            al=adj[last]
            for nxt in range(n):
                if mask&(1<<nxt): continue
                if al[nxt]: dp[mask|(1<<nxt)][nxt]+=cur
    return sum(dp[full][v] for v in range(n))
def sk(adj): return (H(adj),c3(adj),scoreseq(adj))

FORBIDDEN={(9,3,(1,1,2,3,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))}

def main():
    res_vecs=list(combinations_with_replacement(range(0,14),5))
    print(f"residue vectors (multisets of 5 from 0..13): {len(res_vecs)}", flush=True)
    produced=set()
    hits={f:[] for f in FORBIDDEN}
    perms=list(permutations(range(5)))  # 120 relative speed orders
    processed=0
    for res in res_vecs:
        refs=set(res)
        for ref in refs:
            Vm=vote_matrix(res,ref)
            # find tie pairs
            tie_pairs=[(i,j) for i in range(5) for j in range(i+1,5) if Vm[i][j]==0]
            if not tie_pairs:
                # deterministic tournament
                adj=[[0]*5 for _ in range(5)]
                for i in range(5):
                    for j in range(i+1,5):
                        if Vm[i][j]>0: adj[i][j]=1
                        else: adj[j][i]=1
                k=sk(adj); produced.add(k)
                if k in FORBIDDEN: hits[k].append((res,ref,"no-tie"))
            else:
                # orientation of tie pairs depends on speed order; enumerate all orders
                # but distinct orders give same tie-orientation if they agree on tie pairs.
                seen_orient=set()
                for to in perms:
                    rank=[0]*5
                    for pos,idx in enumerate(to): rank[idx]=pos
                    # orientation signature on tie pairs
                    sig=tuple(1 if rank[i]>rank[j] else 0 for (i,j) in tie_pairs)
                    if sig in seen_orient: continue
                    seen_orient.add(sig)
                    adj=[[0]*5 for _ in range(5)]
                    for i in range(5):
                        for j in range(i+1,5):
                            if Vm[i][j]>0: adj[i][j]=1
                            elif Vm[i][j]<0: adj[j][i]=1
                            else:
                                if rank[i]>rank[j]: adj[i][j]=1
                                else: adj[j][i]=1
                    k=sk(adj); produced.add(k)
                    if k in FORBIDDEN: hits[k].append((res,ref,to))
        processed+=1
        if processed%1000==0:
            print(f"  processed {processed}/{len(res_vecs)} residue vectors; produced {len(produced)} keys", flush=True)
    print(f"\nDISTINCT n=5 short-keys produced by M3b (residue-complete, all tie-breaks): {len(produced)}", flush=True)
    print("All produced short-keys (H,c3,score):", flush=True)
    for k in sorted(produced): print("   ",k, flush=True)
    print("\nForbidden-target check:", flush=True)
    refuted=False
    for f in sorted(FORBIDDEN):
        if hits[f]:
            refuted=True
            print(f"  *** REFUTED *** {f}: e.g. {hits[f][0]}  (total witnesses {len(hits[f])})", flush=True)
        else:
            print(f"  CONFIRMED unreachable (top-level n=5): {f}", flush=True)
    print(f"\nOverall n=5 top-level verdict: {'REFUTED' if refuted else 'forbidden classes hold'}", flush=True)

if __name__=="__main__":
    main()
