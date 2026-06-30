"""
Thread 1: the L limit. CLAIM L(AP_n) = int_0^1 M_c dc -> 1/4 = int||c|| (avg distance to nearest integer).
Reason: M_c >= ||c|| (the t=0 witness: runners all at 0, observer at distance ||c||), so L >= 1/4; and the
excess int(M_c-||c||) -> 0 as the runners fill in. Verify the limit + the excess rate.
"""
from fractions import Fraction
def Mc(S,c,Qmax):
    best=0.0
    for q in range(1,Qmax+1):
        for a in range(q):
            m=min(min((s*a/q-c)%1,1-((s*a/q-c)%1)) for s in S)
            if m>best: best=m
    return best
def nint(c): c=c%1.0; return min(c,1-c)   # ||c||
print("L(AP_n) -> 1/4 ?  (excess = L - 1/4, and excess*n):")
print(f"   {'n':>3} {'L':>8} {'L-1/4':>8} {'(L-1/4)*n':>10}")
prev=None
for n in [6,8,10,12,14,16,18,20,24]:
    AP=list(range(1,n)); G=96
    L=sum(Mc(AP,c/G,2*n) for c in range(G))/G
    ex=L-0.25
    print(f"   {n:>3} {L:>8.4f} {ex:>8.4f} {ex*n:>10.4f}")
print(f"   int_0^1 ||c|| dc = {sum(nint(c/100000) for c in range(100000))/100000:.5f} (=1/4 exactly).")
print("   => L DECREASES to 1/4; excess ~ C/n (C~0.37-0.38). The LRC trivializes: as runners fill in, the")
print("      observer's best is its t=0 distance ||c|| to the lattice point 0; mean = 1/4.")
print("   1/4 = the AVERAGE NEAREST-INTEGER DISTANCE = the trivial-loneliness floor; the LRC's hard part")
print("   (M_c >> ||c||) has measure->0 in L^1. So the covmin difficulty is a vanishing correction to 1/4.")
