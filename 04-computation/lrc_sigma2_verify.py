from fractions import Fraction
from math import gcd
from itertools import combinations
import sys

def M_exact(S):
    """Exact max-min collar via crossing candidates t=m/d, d in {v_i +- v_j} (and v_i).
    Integer arithmetic in the hot loop; returns a Fraction."""
    n = len(S)
    ds = set()
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a = S[i] + S[j]
            b = S[i] - S[j]
            if a > 0:
                ds.add(a)
            if b > 0:
                ds.add(b)
    for v in S:
        ds.add(v)  # single-speed candidate

    best_num = 0
    best_den = 1
    for d in ds:
        for m in range(1, d):
            mn = d  # min-norm numerator in units of 1/d
            for v in S:
                r = (v * m) % d
                rr = r if r <= d - r else d - r
                if rr < mn:
                    mn = rr
                    if mn == 0:
                        break
            # mn/d > best_num/best_den  <=>  mn*best_den > best_num*d
            if mn * best_den > best_num * d:
                best_num, best_den = mn, d
    return Fraction(best_num, best_den)

def run(k, box):
    floor = Fraction(1, k + 1)
    mediant = Fraction(2, 2 * k + 1)
    best_below = None
    witness = None
    below_count = 0
    below_examples = []
    total = 0
    for combo in combinations(range(2, box + 1), k - 1):
        S = (1,) + combo
        g = 0
        for v in S:
            g = gcd(g, v)
        if g != 1:
            continue
        total += 1
        M = M_exact(list(S))
        if M == floor:
            continue
        if best_below is None or M < best_below:
            best_below = M
            witness = S
        if floor < M < mediant:
            below_count += 1
            if len(below_examples) < 5:
                below_examples.append((str(M), S))
    return best_below, witness, floor, mediant, below_count, below_examples, total

for k, box in [(6, 20), (7, 22), (8, 24), (9, 26)]:
    bb, wit, floor, med, bc, bex, total = run(k, box)
    print(f"k={k} box={box}: floor=1/{k+1}={floor}, mediant=2/{2*k+1}={med}, gcd1-sets={total}")
    print(f"   sigma_2 (min non-tight M) = {bb}  witness={wit}")
    print(f"   below-mediant sets count = {bc}  examples={bex}")
    print(f"   => sigma_2 {'<' if (bb is not None and bb < med) else '>='} mediant")
    sys.stdout.flush()
