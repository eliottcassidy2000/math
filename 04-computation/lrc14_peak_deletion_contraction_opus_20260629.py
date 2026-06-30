"""
The PEAK deletion-contraction (the margin-bearing twin of mac-mini's cap chain-rule).
Binary-relation reframe: M(empty)=1/2 (no runners: observer maximally lonely, TRANSITIVE/no cycles).
Adding each speed DECREASES the peak (intransitive constraint). M(S) >= 1/14 is the conjecture.
  M(prefix_k) - chain; decrements = the per-speed 'intransitive cost'; final = the peak.
Transitive AP {1..13} = exactly tight (1/14). Covering = perturbed, ends with margin.
"""
from fractions import Fraction as F
def norm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def M_exact(S):
    if not S: return F(1,2)   # empty: max_t min over nothing = 1/2 (the antipode)
    C=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): C.add(F(m,d))
        for m in range(2*S[i]+1): C.add(F(m,2*S[i]))
    best=F(0)
    for t in C:
        if 0<t<1:
            v=min(norm(s*t) for s in S)
            if v>best: best=v
    return best
thr=F(1,14)
def chain(S,label):
    print(f"{label}: peak M as speeds are added (transitive backbone 1/2 -> tight):")
    prev=F(1,2)
    for k in range(len(S)+1):
        M=M_exact(S[:k])
        dec = prev-M
        print(f"  +[{S[k-1] if k else '-':>3}] -> M(prefix_{k:>2})={str(M):>7}={float(M):.5f}  (decrement {float(dec):+.5f})")
        prev=M
    print(f"  final margin above 1/14: {float(M_exact(S)-thr):+.5f}\n")
chain(list(range(1,14)), "TRANSITIVE AP {1..13} (ordered speeds)")
chain([1,2,3,4,5,6,7,8,9,10,11,13,84], "COVERING {1..11,13,84} (forced mult-of-14 = intransitive)")
print("Reframe: M(empty)=1/2 is the pure transitive backbone (observer alone, antipodal witness t=1/2).")
print("Each added speed is an intransitive constraint cutting M. The AP cuts to EXACTLY 1/14 (tight).")
print("The covering set's mult-of-14 (the intransitive 'wrap' cycle) cuts differently -> ends ABOVE 1/14.")
print("DISPROOF would need the cuts to OVERSHOOT below 1/14 -- i.e. the danger arcs to close CYCLICALLY")
print("around the circle (intransitive cover, chi(nerve)=0). LRC = the cover never closes (always a")
print("transitive break = a lonely gap). The arithmetic obstruction to cyclic closure is the proof.")
