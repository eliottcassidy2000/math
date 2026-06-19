from fractions import Fraction
from math import gcd
from functools import reduce
import random

def norm(x):
    # ||x|| = distance to nearest integer, x a Fraction
    f = x - int(x)  # fractional part toward 0; handle negatives
    f = x % 1
    return min(f, 1 - f)

def M_method1(S):
    """Exact M(S) via crossing candidates t = m/d, d in {v_i+-v_j}."""
    S = list(S)
    ds = set()
    for i in range(len(S)):
        for j in range(len(S)):
            if i == j:
                continue
            ds.add(abs(S[i] + S[j]))
            ds.add(abs(S[i] - S[j]))
    ds.discard(0)
    best = Fraction(-1)
    best_t = None
    for d in ds:
        for m in range(1, d):
            t = Fraction(m, d)
            val = min(norm(Fraction(v) * t) for v in S)
            if val > best:
                best = val
                best_t = t
    return best, best_t

def divides_some_pair(q, S):
    S = list(S)
    for i in range(len(S)):
        for j in range(len(S)):
            if i == j:
                continue
            if (S[i] + S[j]) % q == 0:
                return True
            if abs(S[i] - S[j]) % q == 0 and abs(S[i]-S[j]) != 0:
                return True
    return False

def lcm_gcd1(S):
    return reduce(gcd, S)

# ---- Part (a): denominator lemma on 30 random gcd-1 sets ----
random.seed(12345)
print("=== Part (a): denominator lemma ===")
violations = []
count = 0
attempts = 0
while count < 30 and attempts < 100000:
    attempts += 1
    k = random.randint(4, 8)
    speeds = random.sample(range(1, 41), k)
    if lcm_gcd1(speeds) != 1:
        continue
    S = sorted(speeds)
    Mval, t = M_method1(S)
    q = Mval.denominator
    ok = divides_some_pair(q, S)
    qbound = q <= 2 * max(S)
    if not (ok and qbound):
        violations.append((S, Mval, q, ok, qbound))
    count += 1
    print(f"S={S} M={Mval} q={q} divides_pair={ok} q<=2max={qbound}")

print(f"\nTotal sets tested: {count}")
print(f"Violations: {len(violations)}")
for v in violations:
    print("VIOLATION:", v)

# ---- Part (b): D_k = {1,2,...,k-1, 2k}, M = 2/(2k+1) ----
print("\n=== Part (b): doubled apex D_k ===")
bmismatch = []
for k in range(2, 21):
    S = list(range(1, k)) + [2*k]
    # ensure gcd 1
    g = lcm_gcd1(S)
    Mval, t = M_method1(S)
    expected = Fraction(2, 2*k+1)
    match = (Mval == expected)
    if not match:
        bmismatch.append((k, S, Mval, expected))
    print(f"k={k} D_k={S} gcd={g} M={Mval} expected=2/{2*k+1}={expected} match={match} t*={t}")

print(f"\nMismatches in (b): {len(bmismatch)}")
for m in bmismatch:
    print("MISMATCH:", m)
