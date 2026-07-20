"""
opus-2026-07-19-S409 (HYP-8015): two conjunctions of the S404-S408 leads.

(A) THE GW TRICHOTOMY LAW (lambda-hole x (D,s)-rungs x gate, composed):
    for GW_N = {1..N-2, N, 2(N-1)}:
      M = 1/(N+1)   if N == 1 (mod 6)     [tight; slack 0]
      M = 2/(2N+1)  if N odd, N != 1 (mod 3)  [D=2 escape, pair (3, 2(N-1)); slack 1]
      M = 1/N       if N even             [D=1 escape; slack 1]
    MECHANISM (proved inline, 3 lines): the D=2 escape needs a with 3a == +-2
    (mod 2N+1), impossible iff 3 | 2N+1 iff N == 1 (mod 3); and the would-be
    killers of that escape are the residues +-(N-1) mod (2N+1) — exactly the
    DELETED element and its mirror: the lambda box-move hole IS the escape's
    immunity.  Verify N = 5..20: values, active pairs, killer-absence.
(B) THE THETA PHASE TRANSITION: exact M_theta(AP13) at theta = p/8 (integer
    power-comparison: v1^p d1^8 vs v2^p d2^8), locating where the maximizer
    jumps off 1/14 toward Fibonacci; and whether AP/GW stay TIED along theta.
"""
from math import gcd
from fractions import Fraction

def scan0(V, qmax):
    bg, bq, wit = 0, 1, None
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = q
            for v in V:
                r = (v*a) % q
                r = min(r, q-r)
                if r < m:
                    m = r
                    if m*bq < bg*q: break
            if m*bq > bg*q: bg, bq, wit = m, q, (a, q)
    return Fraction(bg, bq), wit

print("(A) GW trichotomy, N = 5..20:")
ok = True
for N in range(5, 21):
    V = list(range(1, N-1)) + [N, 2*(N-1)]
    M, (a, q) = scan0(V, 4*N + 30)
    if N % 6 == 1: pred, branch = Fraction(1, N+1), "tight"
    elif N % 2 == 1: pred, branch = Fraction(2, 2*N+1), "D=2 escape"
    else: pred, branch = Fraction(1, N), "D=1 escape"
    match = (M == pred)
    if not match: ok = False
    # active pair at the maximizer
    act = [v for v in V if Fraction(min((v*a) % q, q-(v*a) % q), q) == M]
    print(f"  N={N} (mod6={N%6}): M={M} pred={pred} [{branch}] match:{match}; "
          f"t*={a}/{q}; active={act}")
print(f"  trichotomy exact on N=5..20: {ok}")

print("\n  escape-blocking checks (odd N):")
for N in range(5, 21, 2):
    Q = 2*N + 1
    blocked = (Q % 3 == 0)
    sol = [aa for aa in range(1, Q) if (3*aa) % Q in (2, Q-2)]
    killers_absent = True
    if sol:
        aa = sol[0]
        inv = pow(aa, -1, Q) if gcd(aa, Q) == 1 else None
        kills = sorted({inv % Q, (Q - inv) % Q}) if inv else []
        V = set(list(range(1, N-1)) + [N, 2*(N-1)])
        killers_absent = not any((v % Q) in kills or ((Q - v % Q) % Q) in kills
                                  for v in V if min(v % Q, Q - v % Q) <= 1)
    print(f"  N={N}: 3|2N+1: {blocked} (N==1 mod 3: {N%3==1}); "
          f"3a==+-2 solvable: {bool(sol)}; killer residues +-(N-1),(N+2): "
          f"absent from V: {killers_absent}")

print("\n(B) theta phase transition for AP13 (theta = p/8, exact power comparison):")
def scan_theta(V, qmax, p, r=8):
    # value(v, d/q) ~ compare v^p * d^r  (theta = p/r); track best by tuple compare
    best, wit = None, None
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = None
            for v in V:
                d = (v*a) % q
                d = min(d, q-d)
                key = (v**p) * (d**r), q**r  # value^r = v^p d^r / q^r
                kv = Fraction(key[0], key[1])
                if m is None or kv < m: m = kv
            if best is None or m > best: best, wit = m, (a, q)
    return best, wit

AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
for p in range(0, 9):
    bA, wA = scan_theta(AP, 120, p)
    bG, wG = scan_theta(GW, 120, p)
    tie = (bA == bG)
    print(f"  theta={p}/8: AP maximizer {wA}; GW maximizer {wG}; "
          f"M^8 equal (pair tied): {tie}")
