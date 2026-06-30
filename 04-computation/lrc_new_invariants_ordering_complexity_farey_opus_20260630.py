"""
Wield definitions. New LRC invariants from the observer:
 O(S) = ORDERING COMPLEXITY = # distinct circular orderings of {0,s_i t} as t varies (LRC analog of Ham-paths).
 L(S) = LONELINESS INTEGRAL = int_0^1 M_c(S) dc (translation-averaged loneliness).
 Conjecture: O(AP_n) = Phi(n-1) = sum_{d<=n-1} phi(d) = Farey-sequence length (totient summatory).
"""
import math
from fractions import Fraction
def totient(d): return sum(1 for k in range(1,d+1) if math.gcd(k,d)==1)
def Phi(N): return sum(totient(d) for d in range(1,N+1))
def ordering_complexity(speeds, n):
    # orderings change at crossings t=j/d, d=|s_i-s_j|. Enumerate exact crossing rationals + midpoints.
    pts=[0]+list(speeds)
    crossings=set([Fraction(0)])
    for i in range(n):
        for j in range(i+1,n):
            d=abs(pts[i]-pts[j])
            if d==0: continue
            for k in range(0,d): crossings.add(Fraction(k,d))
    cr=sorted(crossings)
    seen=set()
    for a in range(len(cr)):
        t=(cr[a]+(cr[a+1] if a+1<len(cr) else 1))/2  # midpoint of each interval
        order=sorted(range(n), key=lambda x:(pts[x]*t)-(pts[x]*t).__floor__())
        z=order.index(0); rot=tuple(order[z:]+order[:z])
        seen.add(rot)
    return len(seen)
print("O(S) = ordering complexity vs Phi(n-1) = Farey length (totient summatory):")
print(f"{'n':>3} {'O(AP_n)':>9} {'Phi(n-1)':>9} {'match':>6} {'O(scaled AP)':>13} {'(dilation-inv?)'}")
for n in range(4,11):
    AP=list(range(1,n)); o=ordering_complexity(AP,n); ph=Phi(n-1)
    p=2 if n%2==0 else (3 if n%3==0 else n); sc=[p*k for k in range(1,n)]; osc=ordering_complexity(sc,n)
    print(f"{n:>3} {o:>9} {ph:>9} {str(o==ph):>6} {osc:>13}   {osc==o}")
print("  => O(AP_n)=Phi(n-1)=FAREY LENGTH; dilation-invariant (O(pS)=O(S)); EVEN (snapshots come in time-")
print("     reverse pairs) -- vs tournament H ODD (Redei). The LRC 'Ham-paths' = Farey; the tournament's = H.")
print()
# Loneliness integral L(S)=int M_c dc (numeric)
def Mc(S,c,Qmax):
    best=0.0
    for q in range(1,Qmax+1):
        for a in range(0,q):
            mind=min(min((s*a/q-c)%1,1-((s*a/q-c)%1)) for s in S)
            if mind>best: best=mind
    return best
print("L(S) = loneliness integral int_0^1 M_c dc (translation-averaged loneliness):")
for n in [6,8,10,14]:
    AP=list(range(1,n)); G=60
    L=sum(Mc(AP,c/G,2*n) for c in range(G))/G
    print(f"  n={n}: L(AP) ~ {L:.4f}  (vs M_0=1/n={1/n:.4f} [hardest] and antipode 1/2; L is the mean)")
print("  => L(S) averages the loneliness over all observer positions; floor 1/n at c=0, the rest higher.")
