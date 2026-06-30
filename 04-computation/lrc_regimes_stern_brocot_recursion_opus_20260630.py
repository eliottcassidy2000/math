"""
Thread 2: continued-fraction / Stern-Brocot refinement of the regime count, recursively.
 - Mode A recursion: O(n)-O(n-1) = phi(n-1) (the new Farey fractions of denom n-1) -- vertex insertion.
 - regimes = Farey intervals = Stern-Brocot NODES; the three-distance gaps per regime come from the
   endpoints' denominators (CF convergents). Mediant = the recursive refinement.
"""
import math
from fractions import Fraction
def totient(d): return sum(1 for k in range(1,d+1) if math.gcd(k,d)==1)
def Phi(N): return sum(totient(d) for d in range(1,N+1))
print("(A) Mode-A recursion: O(AP_n)=Phi(n-1); O(n)-O(n-1)=phi(n-1) (the phi(n-1) NEW regimes per vertex):")
print(f"   {'n':>3} {'O=Phi(n-1)':>11} {'O(n)-O(n-1)':>12} {'phi(n-1)':>9} {'match':>6}")
for n in range(4,13):
    O=Phi(n-1); d=O-Phi(n-2); ph=totient(n-1)
    print(f"   {n:>3} {O:>11} {d:>12} {ph:>9} {str(d==ph):>6}")
print("   => each vertex-insertion (n-1->n) splits phi(n-1) regimes via Farey MEDIANTS (Stern-Brocot growth).")
print()
# Farey F_{N} and the mediant tree; gaps at a regime from endpoint denominators
def farey(N):
    fs=sorted(set(Fraction(a,b) for b in range(1,N+1) for a in range(0,b+1)))
    return fs
print("(B) regimes = Farey intervals; three-distance gaps come from the endpoint denominators (q,q',q+q'):")
for n in [6,8]:
    N=n-1; fs=farey(N)
    print(f"   n={n} (N={N}): {len(fs)-1} Farey intervals (=O); sample intervals (p/q,p'/q') -> mediant, gaps:")
    for i in range(min(4,len(fs)-1)):
        lo,hi=fs[i],fs[i+1]; med=Fraction(lo.numerator+hi.numerator, lo.denominator+hi.denominator)
        t=med  # representative in the interval
        xs=sorted((k*t)%1 for k in range(n))
        gaps=sorted(set((xs[(j+1)%n]-xs[j])%1 for j in range(n)))
        print(f"      ({lo},{hi}) mediant={med}: denoms({lo.denominator},{hi.denominator},{med.denominator}); gaps={[str(g) for g in gaps]}")
print("   => the 3 gaps are governed by the endpoint denominators q,q' and the mediant q+q' (Stern-Brocot);")
print("      the regime's gap-pattern IS the local continued-fraction state. CF doesn't change the COUNT")
print("      (still Phi(n-1)) but gives each regime a Stern-Brocot ADDRESS + the recursive mediant refinement.")
print()
print("(C) the covering-min optimum t=1/n has CF [0;n] -- the SIMPLEST regime (a leaf at the n-th level);")
print("    M=1/n is the gap at the equally-spaced node. The escape lives at a Farey/Stern-Brocot vertex.")
