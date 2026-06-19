from fractions import Fraction
from math import gcd
from itertools import combinations
import sys

def M_with_prune(S, floor_num, floor_den, med_num, med_den):
    """Return (M_exact_fraction_or_None_if_>=med, found_ge_med_bool).
    We only need M precisely when M < mediant. If at any candidate t the
    min-norm >= mediant, then M >= mediant: we can stop and report 'not below'.
    floor/med given as numerator/denominator."""
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
        ds.add(v)

    best_num = 0
    best_den = 1
    for d in ds:
        for m in range(1, d):
            mn = d
            for v in S:
                r = (v * m) % d
                rr = r if r <= d - r else d - r
                if rr < mn:
                    mn = rr
                    if mn == 0:
                        break
            # min-norm = mn/d. If mn/d >= med  => M >= med, abort (not below).
            # mn/d >= med_num/med_den  <=>  mn*med_den >= med_num*d
            if mn * med_den >= med_num * d:
                return None, True
            # track best
            if mn * best_den > best_num * d:
                best_num, best_den = mn, d
    # M = best_num/best_den, which is < mediant here
    return Fraction(best_num, best_den), False

def run(k, box):
    floor = Fraction(1, k + 1)
    mediant = Fraction(2, 2 * k + 1)
    fn, fd = floor.numerator, floor.denominator
    mn_, md = mediant.numerator, mediant.denominator
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
        M, ge = M_with_prune(list(S), fn, fd, mn_, md)
        if ge:
            continue  # M >= mediant, not relevant to sigma_2 < mediant question
        # here M < mediant. M could be == floor (tight) or strictly between.
        if M == floor:
            continue
        # non-tight AND below mediant
        below_count += 1
        if best_below is None or M < best_below:
            best_below = M
            witness = S
        if len(below_examples) < 10:
            below_examples.append((str(M), S))
    return best_below, witness, floor, mediant, below_count, below_examples, total

k, box = 9, 26
bb, wit, floor, med, bc, bex, total = run(k, box)
print(f"k={k} box={box}: floor=1/{k+1}={floor}, mediant=2/{2*k+1}={med}, gcd1-sets={total}")
print(f"   below-mediant non-tight sets count = {bc}")
print(f"   min below-mediant M = {bb}  witness={wit}")
print(f"   examples={bex}")
if bc == 0:
    print(f"   => NO below-mediant set; sigma_2 = mediant = {med}")
else:
    print(f"   => below-mediant set EXISTS; sigma_2 = {bb} < {med}")
sys.stdout.flush()
