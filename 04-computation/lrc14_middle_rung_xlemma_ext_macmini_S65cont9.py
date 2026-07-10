from itertools import combinations
from math import gcd
from functools import reduce
# extend the X-lemma check to D = 22..25 (beyond the rung; probes the B>=33 zone)
# + verify the OMEGA-MAP facts: (F1) low-missing s => s in H; (F2) high-missing D+u => u in H;
# (F3) isolated a (a != D/2) => its canonical partner 2a (a<D/2) / 2a-D (a>D/2) is a hole;
# (F4) count mixed omega-collisions (D even, a and a+D/2 both isolated).
violX = viol1 = viol2 = viol3 = 0
mixed = 0
tot = 0
for D in range(22, 26):
    h = D - 12
    for mid in combinations(range(1, D), 11):
        A = (0,) + mid + (D,)
        if reduce(gcd, A[1:]) != 1: continue
        tot += 1
        Aset = set(A)
        S = set(); AA = set()
        for i in range(13):
            for j in range(i, 13):
                AA.add(A[i] + A[j])
                if i < j: S.add(A[i] + A[j])
        X = sum(1 for i in range(1, 12) if 2*A[i] not in S)
        missing = (2*D + 1) - len(AA)
        if missing > h - X: violX += 1; print("X-VIOL", A)
        for s in range(1, 2*D):
            if s in AA: continue
            if s < D and s in Aset: viol1 += 1
            if s > D and (s - D) in Aset: viol2 += 1
        iso = [a for a in A[1:12] if 2*a not in S]
        for a in iso:
            if 2*a == D: continue
            part = 2*a if 2*a < D else 2*a - D
            if part in Aset or part == 0: viol3 += 1; print("F3-VIOL", A, a)
        if D % 2 == 0:
            for a in iso:
                if a + D//2 in iso: mixed += 1
print(f"D in [22,25]: {tot} sets; X-lemma violations {violX}; F1 {viol1}; F2 {viol2}; F3 {viol3}; mixed omega-collisions {mixed}")
