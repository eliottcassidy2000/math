"""SIGNED LRC, rigorous bottom-up. A sign vector ε∈{±1}^{n-1} = a 2-COLORING of the runners.
Pair (i,j) clock = ε_i v_i − ε_j v_j ∈ {±(v_i−v_j) [monochromatic], ±(v_i+v_j) [bichromatic=cut]}.
Work mod C=2n−1 (THM-401). ZERO clock = pair-sum ≡0 mod C = a SHELL-PARTNER. Observer M is
sign-invariant (theorem). Enumerate ALL sign patterns at each n; classify the sign-orbit, the
cut/sum structure, the zero-clock-exposing signs. opus-2026-06-06-S699c."""
from itertools import combinations, product
def clocks_modC(V, eps, C):
    # signed pairwise frequencies mod C, as a sorted multiset of residues in [0,C)
    n=len(V); out=[]
    for i in range(n):
        for j in range(i+1,n):
            f=(eps[i]*V[i]-eps[j]*V[j])%C
            out.append(min(f, C-f))   # |clock| as distance to 0 mod C (the resonance strength)
    return tuple(sorted(out))
def shell_partners(V,C):
    return [(V[i],V[j]) for i in range(len(V)) for j in range(i+1,len(V)) if (V[i]+V[j])%C==0]
def main():
    print("Each config: C=2n−1; shell-partners (v_i+v_j≡0 mod C); over ALL 2^(n-1) signs:")
    print("  #distinct pair-clock multisets (sign-orbit), max #zero-clocks exposable, maxcut sum-count")
    configs={
      3:[(1,2)],
      4:[(1,2,3)],
      5:[(1,2,3,4),(1,2,4,7)],
      6:[(1,2,3,4,5),(2,3,4,5,9)],          # AP (no shell-partner); (2..9) has 2+9=11 shell-partner
      7:[(1,2,3,4,5,6),(1,2,3,4,9,12)],     # AP; a shell-partner config (1+12=13, 4+9=13)
      8:[(1,2,3,4,5,6,7),(1,2,5,6,9,12,13)],
    }
    for n in sorted(configs):
        C=2*n-1
        for V in configs[n]:
            sp=shell_partners(V,C)
            orbit=set(); maxzero=0; maxcut=0
            for eps in product((1,-1),repeat=len(V)):
                cl=clocks_modC(V,eps,C); orbit.add(cl)
                z=sum(1 for x in cl if x==0); maxzero=max(maxzero,z)
                # cut size = # bichromatic pairs
                cut=sum(1 for i in range(len(V)) for j in range(i+1,len(V)) if eps[i]!=eps[j])
                maxcut=max(maxcut,cut)
            print(f"  n={n} V={V}: C={C}; shell-partners={sp if sp else 'none'}; "
                  f"sign-orbit={len(orbit)} (of 2^{len(V)}={2**len(V)}); max zero-clocks={maxzero}; maxcut={maxcut}")
    print("\nRIGOROUS: sign=2-coloring; sum-clock iff bichromatic (cut edge), diff-clock iff mono.")
    print("Zero clock mod C ⟺ shell-partner (v_i+v_j≡0 mod 2n-1, THM-401). AP {1..n-1}: max sum 2n-3<C,")
    print("so AP has NO shell-partner / NO zero clock. Configs with a shell-partner expose a zero clock")
    print("under the cut separating that pair — a sign-orbit invariant FINER than the observer M.")
if __name__=='__main__': main()
