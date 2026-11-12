"""The O/Z/J system as a linear semigroup on states (alpha, M, E, d) of a rooted tournament:
   s1 (tail letter 0: new vertex z, B -> z, z -> v):  (a + M + d - E, 2M, d, d)
   s0 (tail letter 1: v -> z):                        (a, 2M + d, 0, d)
   rho (tail letter 0, then re-root at z):            (a + M + d - E, M + d, a - E, a)
1) verify against literal counts; 2) search for common eigenvectors of words (invariant lines)."""
import itertools, random, sys
from fractions import Fraction
import sympy as sp
S1=sp.Matrix([[1,1,-1,1],[0,2,0,0],[0,0,0,1],[0,0,0,1]])
S0=sp.Matrix([[1,0,0,0],[0,2,0,1],[0,0,0,0],[0,0,0,1]])
RH=sp.Matrix([[1,1,-1,1],[0,1,0,1],[1,0,-1,0],[1,0,0,0]])
GEN={'1':S1,'0':S0,'r':RH}
def literal(word):
    # tournament as out-lists; root index
    out={0:set()}; root=0; n=1
    for ch in word:
        z=n; n+=1; out[z]=set()
        for b in range(z):
            if b!=root: out[b].add(z)
        if ch in '1r': out[z].add(root)       # z -> v
        else: out[root].add(z)                # v -> z   (letter '0' = s0)
        if ch=='r': root=z
    # count Hamiltonian paths
    masks=[sum(1<<w for w in out[v]) for v in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for Sx in range(1,1<<n):
        row=dp[Sx]
        for v in range(n):
            c=row[v]
            if c:
                t=masks[v]&~Sx
                while t:
                    b=t&-t; dp[Sx|b][b.bit_length()-1]+=c; t^=b
    return sum(dp[(1<<n)-1])
def state(word):
    v=sp.Matrix([1,0,1,1])
    for ch in word: v=GEN[ch]*v
    return v
if __name__=="__main__":
    rng=random.Random(1); ok=0; tot=0
    for _ in range(300):
        w=''.join(rng.choice('10r') for _ in range(rng.randint(1,10)))
        tot+=1; ok+=(state(w)[0]==literal(w))
    print("linear system vs literal count on random words (length <= 10): %d/%d" % (ok,tot))
    # value set of short words
    vals=set()
    for L in range(0,11):
        for w in itertools.product('10r',repeat=L): vals.add(int(state(''.join(w))[0]))
    odd=[v for v in sorted(vals) if v<=2000]
    print("words of length <= 10: %d values; odd numbers <= 400 missing: %s" % (len(vals), [x for x in range(1,400,2) if x not in vals]))
