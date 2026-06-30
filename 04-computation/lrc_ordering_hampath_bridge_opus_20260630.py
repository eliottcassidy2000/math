"""Honest bridge. Recompute O across sets: does AP minimize or maximize O? Then build the bridge on facts."""
import math
from fractions import Fraction
from itertools import permutations, combinations
def O_complexity(speeds,n):
    pts=[0]+list(speeds)
    cr=set([Fraction(0)])
    for i in range(n):
        for j in range(i+1,n):
            d=abs(pts[i]-pts[j])
            if d:
                for k in range(d): cr.add(Fraction(k,d))
    cr=sorted(cr); seen=set()
    for a in range(len(cr)):
        t=(cr[a]+(cr[a+1] if a+1<len(cr) else 1))/2
        order=sorted(range(n), key=lambda x:(pts[x]*t)-(pts[x]*t).__floor__())
        z=order.index(0); seen.add(tuple(order[z:]+order[:z]))
    return len(seen)
def Phi(N): return sum(sum(1 for k in range(1,d+1) if math.gcd(k,d)==1) for d in range(1,N+1))
n=5
print(f"n={n}: O(S) for AP {{1..{n-1}}} vs other {n-1}-sets of distinct positive integers:")
AP=list(range(1,n)); print(f"   AP {AP}: O={O_complexity(AP,n)} = Phi(n-1)=Phi({n-1})={Phi(n-1)}")
for S in [[1,2,3,5],[1,2,4,5],[1,3,5,7],[2,3,4,5],[1,2,3,6]]:
    print(f"   {S}: O={O_complexity(S,n)}")
# is AP the MIN over all sets with max speed = n-1? and over small sets?
mn=99; mnS=None
for S in combinations(range(1,9),n-1):
    o=O_complexity(list(S),n)
    if o<mn: mn=o; mnS=S
print(f"   MIN O over {n-1}-subsets of [1..8] = {mn} at {mnS} (AP? {list(mnS)==AP})")
print()
print("CORRECTED FACT: the AP MINIMIZES O (smallest max speed n-1 => smallest Farey F_{n-1} => O=Phi(n-1)).")
print("  My earlier 'AP maximizes O / doubly extremal' was WRONG -- AP is doubly MINIMAL (min M, min O).")
print()
print("THE BRIDGE (honest, structural -- not a value equality):")
print("  O(S) = #orders the runners realize over TIME (the snapshot family) -- EVEN (time-reverse t<->1-t closed),")
print("         a closed walk on the permutohedron; for the AP O=Phi(n-1) (Farey), the MIN.")
print("  H(T) = #orders the arcs realize by CONSISTENCY (Ham paths) -- ODD (Redei), not reversal-closed.")
print("  BOTH count realized linear orders of n elements with a marked observer; the LRC by TIME, the")
print("  tournament by ARCS. Duality: TIME<->ARCS, EVEN<->ODD, Farey<->OCF, family<->consistency.")
print("  EXTREMAL MATCH: AP = min O = min M (simplest LRC) <-> TRANSITIVE = min H=1 (simplest tournament).")
print("  No clean O=H (parity); the bridge is the ANALOGY + the time/arcs duality + the matched minima.")
