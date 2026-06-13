"""'A loop of the input swaps the two outputs' = MONODROMY of the FTA coeff→root map; the branch
locus = the discriminant = the worry-set (roots collide, S699l). The LRC runners = a SPACETIME
BRAID; crossings (pair-clock zeros, S699) = the swaps. 'Rado graph as a tournament' = the countable
RANDOM tournament (Fraïssé limit), the GENERIC fiber with the extension property = 'always a
witness'. opus-2026-06-06-S699p."""
import cmath, math, random
def main():
    print("(1) MONODROMY: loop the input around the discriminant ⟹ the two roots SWAP.")
    print("    z² − t : roots ±√t. Loop t = e^{iθ}, θ: 0→2π, track the root starting at +1:")
    r=1.0+0j
    for k in range(0,13):
        th=2*math.pi*k/12; t=cmath.exp(1j*th)
        # continuously-tracked sqrt: pick the branch nearest previous r
        cand=[cmath.sqrt(t), -cmath.sqrt(t)]
        r=min(cand, key=lambda c: abs(c-r))
        if k in (0,3,6,9,12): print(f"      θ={th:.3f}: tracked root = {r:.3f}")
    print(f"    ⟹ started at +1, returned at {r:.3f} = −1: the two roots SWAPPED (a Galois transposition / braid half-twist).")
    print("    The branch point t=0 is the DISCRIMINANT (roots collide) = the worry-set analog (S699l: roots collide = multiplicities).\n")
    print("(2) LRC SPACETIME BRAID: runners v_i t mod 1 are strands; crossings (v_i−v_j)t∈ℤ are the swaps.")
    V=[1,2,3]
    crossings=[]
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            d=abs(V[i]-V[j])
            for k in range(1,d+1): crossings.append((k/d if d else 0,(V[i],V[j])))
    crossings.sort()
    print(f"    V={V}: crossings (t, pair) over one period = {[(round(t,3),p) for t,p in crossings]}")
    print(f"    each crossing = a transposition of two runners (a braid generator); the pair-clock=0 of S699.")
    print(f"    # crossings = Σ|v_i−v_j| = {sum(abs(V[i]-V[j]) for i in range(len(V)) for j in range(i+1,len(V)))} (the braid's length / the writhe).\n")
    print("(3) RADO = the countable RANDOM tournament (Fraïssé limit); EXTENSION PROPERTY = 'always a witness':")
    rng=random.Random(0)
    for n in (10,20,40,80):
        # random tournament; for random disjoint A,B check ∃ vertex beating all A, beaten by all B
        adj=[[0]*n for _ in range(n)]
        for a in range(n):
            for b in range(a+1,n):
                if rng.random()<.5: adj[a][b]=1
                else: adj[b][a]=1
        ok=0; trials=200
        for _ in range(trials):
            S=rng.sample(range(n), min(4,n)); A=set(S[:2]); B=set(S[2:])
            found=any(all(adj[w][a] for a in A) and all(adj[b][w] for b in B) for w in range(n) if w not in A|B)
            ok+=found
        print(f"    n={n}: P(∃ witness beating A, beaten by B) ≈ {ok/trials:.3f}  → 1 (the extension property, the Rado/random tournament limit)")
    print("\n=> SYNTHESIS: witness-finding is a COVERING MAP; 'loop input ⟹ swap outputs' is its MONODROMY")
    print("   (Galois/braid, transpositions at the discriminant = worry-set). The Rado/random tournament")
    print("   is the GENERIC FIBER (extension property = always a witness, full symmetry); the worry-set is")
    print("   the SPECIAL fiber (witnesses collide, branch points). LRC 'always lonely' = the extension property;")
    print("   the worry-set = where it degenerates (measure-0, the discriminant/branch).")
if __name__=='__main__': main()
