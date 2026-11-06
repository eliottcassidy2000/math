# Independent audit of THM-4887, part 8: reproduce the checks-section totals of THM-4887 lines 176-177, which the
# committed d_lemma.py does not print.  "lemma totals" = sum over ordered pairs h != h' of the strip-lemma bound
# (min-b Dirichlet choice); "actual-shift totals" = sum over pairs of Z(h,h') at the actual shift D = g(h') - g(h).
import numpy as np, math

def g(v):
    G, A = 1, 0
    for ch in bin(v)[2:]:
        if ch == "1": A = 2 * A + 1
        else: G += A + (A & 1); A *= 2
    return G

def beta_minb(h, hp, B):
    for b in range(1, B + 1):
        a = (2 * b * h + hp) // (2 * hp)
        if B * abs(b * h - a * hp) < hp:
            return b // math.gcd(a, b)

for (n, m) in [(10, 2), (11, 2), (12, 3), (13, 3)]:
    C = 1 << m; P = 1 << (n - m - 1); B = P // C
    c = np.arange(C, dtype=np.int64)
    out = {}
    for W in (C, 2 * C):
        act = 0; lem = 0.0
        for h in range(P, 2 * P):
            gh = g(h)
            for hp in range(P, 2 * P):
                if hp == h: continue
                D = g(hp) - gh
                act += int((np.abs(h * c[:, None] - hp * c[None, :] - D) <= W).sum())
                b = beta_minb(h, hp, B)
                lem += (2 * C * C / P + 3 * C / b + 5) if W == C else (4 * C * C / P + 3 * C / b + 7)
        out[W] = (act, lem)
    print(f"(n,m)=({n},{m}): W=C: actual-shift total {out[C][0]}, lemma total {out[C][1]:.0f};  "
          f"W=2C: actual-shift total {out[2*C][0]}, lemma total {out[2*C][1]:.0f}", flush=True)
