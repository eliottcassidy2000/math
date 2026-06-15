"""
ADVERSARIAL independent verification of the 'forbidden-7-21' worker claims.
Written from scratch — does NOT import or reuse the worker's scripts.

Claims to check:
 (A) OCF: H(T) = #directed Hamiltonian paths = I(Omega,2) = sum_k alpha_k 2^k.
 (B) I(K3,2) = 1+2*3+4*0 = 7.
 (C) 7 and 21 absent from achievable H spectrum at n<=6 (exhaustive).
 (D) H=21 has FOUR graph-realizable cluster profiles [1,4,3],[1,6,2],[1,8,1],[1,10].
 (E) H=7 has EXACTLY ONE graph-realizable profile [1,3,0]=K3.
 (F) strong-min(m) = 3,5,9,15 for m=3..6 (exhaustive over strong tournaments).
"""
import itertools
import sys

def all_tournaments(n):
    """Yield adjacency matrices (tuple of tuples) of all tournaments on n vertices.
    Edge orientation for pair (i<j): bit=1 means i->j, bit=0 means j->i."""
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(1<<m):
        A = [[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1:
                A[i][j]=1
            else:
                A[j][i]=1
        yield tuple(tuple(r) for r in A)

def ham_paths(A):
    """Count directed Hamiltonian paths in tournament A (n! / brute over perms,
    but use DP for speed)."""
    n=len(A)
    # DP over subsets: dp[mask][v] = number of directed paths covering 'mask' ending at v
    full=(1<<n)-1
    # dp as dict-of-arrays
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            cur=dp[mask][v]
            if cur==0: continue
            if not (mask>>v)&1: continue
            for w in range(n):
                if (mask>>w)&1: continue
                if A[v][w]==1:
                    dp[mask|(1<<w)][w]+=cur
    return sum(dp[full][v] for v in range(n))

def directed_odd_cycles(A):
    """Return list of directed odd cycles, each as a frozenset of vertices.
    A directed cycle = a cyclic sequence v0->v1->...->v_{k-1}->v0 all arcs present.
    We enumerate by vertex subset then count Hamiltonian directed cycles within;
    but for Omega we need the cycle's VERTEX SET (for the conflict graph each
    distinct directed cycle is a node). We must count cycles, not just vertex sets,
    because two distinct directed cycles can share a vertex set.
    Return list of (frozenset_of_vertices) one per distinct directed cycle, but we
    need them as distinct nodes -> represent each directed cycle by its canonical
    rotation-min tuple."""
    n=len(A)
    cycles=[]
    # enumerate over subsets of size>=3 odd
    verts=list(range(n))
    for size in range(3,n+1,2):  # odd sizes only
        for sub in itertools.combinations(verts,size):
            subset=set(sub)
            # find all directed Hamiltonian cycles of the sub-tournament induced on sub
            # fix smallest vertex as start to avoid rotation duplicates
            start=sub[0]
            rest=[x for x in sub if x!=start]
            for perm in itertools.permutations(rest):
                seq=(start,)+perm
                ok=True
                for k in range(len(seq)):
                    a=seq[k]; b=seq[(k+1)%len(seq)]
                    if A[a][b]!=1:
                        ok=False;break
                if ok:
                    # canonical: this directed cycle. To dedup direction-rotation:
                    # rotation fixed by start; reflection is a DIFFERENT directed cycle
                    cycles.append((seq, frozenset(sub)))
    return cycles

def build_omega(A):
    """Build conflict graph Omega: nodes = directed odd cycles, edge iff share>=1 vertex."""
    cyc=directed_odd_cycles(A)
    nodes=list(range(len(cyc)))
    vsets=[c[1] for c in cyc]
    adj=[set() for _ in nodes]
    for i in range(len(cyc)):
        for j in range(i+1,len(cyc)):
            if vsets[i]&vsets[j]:
                adj[i].add(j); adj[j].add(i)
    return adj

def independence_poly_at(adj, z):
    """Compute I(Omega,z)=sum_k alpha_k z^k via counting independent sets.
    For small Omega we can enumerate independent sets exactly using a recursive
    branch (alpha_k counts). We return the full alpha vector AND I(z)."""
    n=len(adj)
    # count independent sets by size via DP/branching
    alpha=[0]*(n+1)
    alpha[0]=1
    # use recursion: pick vertex order; standard independent set enumeration
    order=list(range(n))
    def rec(remaining, chosen):
        # remaining = sorted list of available vertices
        if not remaining:
            return
        # branch on first vertex v: include or exclude
        v=remaining[0]
        rest=remaining[1:]
        # exclude v
        # include v: remove neighbors
        # We count by enumerating all independent sets. To get alpha_k we increment.
        # include
        newrem=[u for u in rest if u not in adj[v]]
        alpha[len(chosen)+1]+=1
        rec(newrem, chosen+[v])
        # exclude
        rec(rest, chosen)
    rec(order, [])
    Iz=sum(alpha[k]*(z**k) for k in range(n+1))
    return alpha, Iz

def main():
    print("=== (A)/(B) OCF identity + spectra at n<=6 ===")
    for n in range(3,7):
        spectrum=set()
        cnt=0
        mismatch=0
        for A in all_tournaments(n):
            cnt+=1
            H=ham_paths(A)
            adj=build_omega(A)
            alpha,Iz=independence_poly_at(adj,2)
            if Iz!=H:
                mismatch+=1
                if mismatch<=3:
                    print(f"  MISMATCH n={n}: H={H} I(Omega,2)={Iz}")
            spectrum.add(H)
        print(f"n={n}: #tournaments={cnt}, OCF mismatches={mismatch}, "
              f"H-spectrum={sorted(spectrum)}, 7 in? {7 in spectrum}, 21 in? {21 in spectrum}")

if __name__=="__main__":
    main()
