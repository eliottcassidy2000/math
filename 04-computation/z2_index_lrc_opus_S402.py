# opus-2026-07-20-S402 -- HYP-8225: the Z/2-INDEX of the LRC configuration.
#
# THM-1380 left the route owing "one involution, free AND carrying an odd map",
# with the concrete next step the k-torus of the resonance lattice and the caveat
# that T^k != S^k blocks plain BU, so one needs the Z/2-index form.
#
# THAT DEBT IS PAYABLE, AND THE ANSWER CLOSES THE ROUTE.  The Fadell-Husseini /
# Yang index of a free Z/2-space X is ind(X) = max{ m : w_1^m != 0 } in
# H*(X/Z2; F_2).  The covering consequence is:
#
#     X free Z/2, ind(X)=m, covered by closed sets none of which contains an
#     antipodal pair  ==>  need at least m+2 sets.
#
# So the index is EXACTLY the strength of every Borsuk-Ulam-family argument here.
# Compute it.
from fractions import Fraction as F
import itertools, random

print("="*74)
print("(1) ind(T^k, half-translation) -- symbolic, in the exterior algebra over F_2")
print("="*74)
print("  s = translation by (1/2,...,1/2) is FREE on T^k, and T^k/s is again a")
print("  k-torus (lattice Z^k + Z*(e/2), index 2).  So")
print("      H*(T^k/s ; F_2) = Lambda_{F2}(x_1,...,x_k),   x_i^2 = 0.")
print("  Write w_1 = sum a_i x_i.  Then")
print("      w_1^2 = sum_{i,j} a_i a_j x_i x_j = sum_{i<j} a_i a_j (x_i x_j + x_j x_i) = 0")
print("  because x_i x_j = -x_j x_i = x_j x_i over F_2, so each bracket is 2*(..)=0.")
# brute-force confirm the exterior-square vanishing over F2 for small k
def sq_is_zero(k):
    for a in itertools.product([0,1],repeat=k):
        acc={}
        for i in range(k):
            for j in range(k):
                if i==j: continue          # x_i^2 = 0
                key=(min(i,j),max(i,j))
                acc[key]=(acc.get(key,0)+a[i]*a[j])%2
        if any(acc.values()): return False
    return True
print(f"\n  brute force over all w_1 in H^1, k=1..8: "
      f"{[sq_is_zero(k) for k in range(1,9)]}")
print("  => w_1^2 = 0 for EVERY w_1, every k.   ind(T^k, translation) = 1.")
print("  (contrast ind(S^k, antipodal) = k, where w_1^k != 0 in H*(RP^k))")
print("  THE INDEX DOES NOT GROW WITH k.  Raising the dimension buys NOTHING.")

print()
print("="*74)
print("(2) The parity dichotomy: which combs are antipodal-free?")
print("="*74)
print("  s: t -> t+1/2.  D_v = {t : ||vt|| < lam}.  Recall ||v(t+1/2)|| = ||vt|| (v even),")
print("  = 1/2 - ||vt|| (v odd).  So D_v contains an s-antipodal pair {t,t+1/2} iff ...")
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
lam=F(1,14)
for v in range(1,15):
    has=False
    for num in range(0,2*7*v*2):
        t=F(num,2*7*v*2)
        if fn(t*v)<lam and fn((t+F(1,2))*v)<lam: has=True;break
    par="even" if v%2==0 else "odd "
    print(f"    v={v:2d} ({par})  contains an antipodal pair: {str(has):5s}"
          f"   {'(D_v is s-INVARIANT)' if v%2==0 else '(D_v is s-FREE)'}")
print("\n  PROVED: for lam < 1/4, both ||vt||<lam and 1/2-||vt||<lam is impossible")
print("  (it would need ||vt|| < lam AND ||vt|| > 1/2-lam, empty when lam<1/4).")
print("  So: v EVEN => D_v is a union of antipodal pairs;  v ODD => D_v antipodal-FREE.")

print()
print("="*74)
print("(3) THE DECISIVE COMPARISON: index bound vs measure bound")
print("="*74)
print("  Lusternik-Schnirelmann-Borsuk, index form: a free Z/2-space X with")
print("  ind(X)=m covered by closed ANTIPODAL-FREE sets needs >= m+2 of them.")
print("  Our configuration space is the circle: ind(S^1, s) = 1.  So the ENTIRE")
print("  Borsuk-Ulam family can only ever say:  >= 3 combs are needed.")
print()
print("  The trivial UNION BOUND says: each comb has measure 2*lam, so covering")
print("  the circle needs >= 1/(2 lam) = n/2 combs.")
print()
print("     n     BU/LSB bound      measure bound      speeds available (n-1)")
for n in [4,5,6,7,8,10,14,20]:
    print(f"    {n:3d}         3               {n/2:>5.1f}                 {n-1}"
          f"   {'<- BU wins' if 3 > n/2 else '<- measure STRICTLY DOMINATES BU'}")
print()
print("  CROSSOVER at n/2 = 3, i.e. n = 6.  For every n >= 7 the measure bound is")
print("  strictly stronger than anything Borsuk-Ulam can give.  At n = 14: BU says")
print("  >= 3, measure says >= 7, and we must rule out 13.  BU is DOMINATED by a")
print("  method THM-1185 already showed is insufficient.")
print()
print("  Worse: LSB constrains only the ANTIPODAL-FREE sets, i.e. only the ODD")
print("  combs (part 2).  {1,...,13} has 7 odd speeds, and 7 >= 3 comfortably.")
print("  So the bound is not merely weak -- it is vacuous on the extremal family.")

print()
print("="*74)
print("(4) Sharpness: 3 antipodal-free closed arcs DO cover S^1 (so 3 is optimal)")
print("="*74)
arcs=[(F(0),F(17,50)),(F(33,100),F(67,100)),(F(66,100),F(1))]
cov=True
for k in range(2000):
    t=F(k,2000)
    if not any(a<=t<=b for a,b in arcs): cov=False;break
free=all(b-a<F(1,2) for a,b in arcs)
print(f"    three arcs of length ~0.34: cover S^1 = {cov}; each length < 1/2 "
      f"(=> antipodal-free) = {free}")
print("    => ind+2 = 3 is attained.  The bound cannot be improved.")

print()
print("="*74)
print("(5) WHY NO ENLARGEMENT ESCAPES: monotonicity of the index")
print("="*74)
print("  If there is a Z/2-equivariant map X -> Y then ind(X) <= ind(Y).")
print("  The LRC data lives on S^1 (index 1).  Any auxiliary free Z/2-space X that")
print("  maps equivariantly to the circle inherits ind(X) <= 1.  You cannot import")
print("  a high-index space into the problem; and mapping OUT (S^1 -> X) leaves the")
print("  conclusion on X, where it says nothing about the combs.")
print("  => The cap of 3 is structural, not an artifact of a bad choice of space.")
