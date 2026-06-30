"""Sanity: the 4 cusps of X_0(14) as a subset of Z/6, and check R(complement) is NOT [-1]."""
# from cuspidal2: base S=d14=0; T=d1 order6 (generator g), +=d2 order3, -=d7 order2.
# realize: g=1 => T=1, +=2 (order3), -=3 (order2), S=0. 4 cusps = {0,1,2,3} in Z/6.
cusps={'S':0,'T':1,'+':2,'-':3}
def order(x): 
    for k in range(1,7):
        if (k*x)%6==0: return k
print("4 cusps as Z/6 elements (base S=0, generator T=1):", cusps)
print("orders:",{k:order(v) for k,v in cusps.items()},"  (T=6 generator, +=3, -=2, S=1)")
print("sum of 4 cusps mod 6:",sum(cusps.values())%6,"(the 4 cusps sum to 0 -- a cuspidal relation)")
neg={k:( -v)%6 for k,v in cusps.items()}
print("\n[-1] (negation) sends: ",{k:f'{v}->{neg[k]}' for k,v in cusps.items()})
print("  [-1] fixes 0(S),3(-); swaps 1(T)<->5, 2(+)<->4.  So [-1] would FIX -, move T,+ OFF the cusp set.")
print("complement R sends:  T->T, S->S, + <-> -  i.e. fixes 1(T),0(S); swaps 2(+)<->3(-).")
print("  => R != [-1]: R swaps an ORDER-3 cusp(+) with an ORDER-2 cusp(-), impossible for a group automorphism.")
print("  HONEST: the complement is NOT the curve negation; SC=2-torsion/pair=3-torsion is FALSE.")
print("  R is a curve involution fixing the T-S (order-6) axis and flipping the +/- pair (an affine reflection).")
print()
print("SOLID: n=4 metagraph (4 classes) embeds in E_tors(14a)=Z/6 (Manin-Drinfeld). The TRANSITIVE class T")
print("       is a GENERATOR (order 6). The 4 cusps {0,1,2,3} sum to 0. 6=phi(14)=units=existence column.")
