"""
Inhomogeneous LRC (observer at c): M_c(S)=max_t min_s ||st-c||. Complement = time-reversal t->-t sends
c->-c. Self-antipodal observers c=0,1/2 (the 2-torsion) <-> the SC (R-fixed) tournaments. Poke: c at the
6-torsion (k/6) <-> 14a torsion Z/6=units; the AP is time-reversal-invariant (SC).
"""
import math
from fractions import Fraction
def Mc(S,c,Qmax):
    # max over t=a/q of min_s ||s*a/q - c||, c=Fraction
    best=Fraction(-1)
    for q in range(1,Qmax+1):
        for a in range(0,q):
            mind=None
            for s in S:
                x=Fraction(s*a,q)-c
                d=x-x.__floor__()  # frac
                d=min(d,1-d)
                if mind is None or d<mind: mind=d
            if mind>best: best=mind
    return best
n=14; AP=list(range(1,n))  # the covmin core (M_0 = 1/14)
print("Inhomogeneous LRC M_c(AP {1..13}) at the 6-torsion c=k/6 (14a torsion Z/6 = units mod 14):")
for k in range(6):
    c=Fraction(k,6); M=Mc(AP,c,60)
    print(f"  c={c}: M_c = {M} = {float(M):.5f}   (c<->-c: -c={(-c)% 1})")
print(f"  c=0 (standard LRC) and c=1/2 (antipode) = the 2-torsion = self-antipodal (c=-c). M_0=1/14 the LRC bound.")
print()
# verify c <-> -c symmetry (complement / time reversal)
print("c <-> -c symmetry (= tournament complement = time reversal t->-t):")
for c in [Fraction(1,6),Fraction(1,3),Fraction(2,7)]:
    print(f"  M_{c} = {Mc(AP,c,40)} ; M_{(-c)%1} = {Mc(AP,(-c)%1,40)} ; equal: {Mc(AP,c,40)==Mc(AP,(-c)%1,40)}")
print()
# the AP is time-reversal-invariant (SC): {1..n-1} reversed mod n = {n-1..1} = same set
print(f"AP {{1..{n-1}}} reversed mod {n}: {sorted((-x)%n for x in AP)} = AP? {sorted((-x)%n for x in AP)==sorted(x%n for x in AP)}")
print("  => the covmin core (AP) is SELF-COMPLEMENTARY (time-reversal-invariant) = the SC/cusp tournament analog;")
print("     the transitive H=1 (ordered baseline) <-> the AP (equally-spaced ordered). Observer c=0 = the SC fixed pt.")
print()
print("CREATIVE CONNECTIONS (poked):")
print("  * complement R (tournament) = time-reversal t->-t (LRC); both act on the observer as c->-c.")
print("  * self-antipodal observers c in {0,1/2} (the 2-TORSION) <-> SC tournaments (R-fixed); c-pairs <-> R-pairs.")
print("  * observer at the 6-torsion c=k/6 <-> 14a torsion Z/6 = units mod 14 = the razor's-edge/binding points.")
print("  * the AP/covmin is time-reversal-invariant (SC) = the 'ordered' fixed config <-> transitive tournament H=1.")
print("  * inhomogeneous M_c over c = a 'defect spectrum'; c=0 gives the LRC bound 1/n, c=1/2 the antipodal bound.")
