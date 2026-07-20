#!/usr/bin/env python3
"""
graded_ladder_repair_boxeph_S183.py  (HYP-8510, THM-1680 section 5)

REPAIRED STAGE B — the sqrt(m)-graded ladder, demonstrated end-to-end on
EXACT geometry data (MISTAKE-204 discipline: every input is an exact moment
of an actual P under the circular-pair functional E[Z^aW^b] = delta 2^a a!;
nothing is synthesized from the model being tested).

Model: P_b = Z + b + W. Exact moments E[P_b^m] = sum_k m!/(k!(m-2k)!) 2^k
b^{m-2k}; EGF e^{2t^2 + bt}; saddle t* ~ sqrt(m)/2 gives the honest scale:
    E[P_b^m]/E[P_0^m] ~ e^{(b/2) sqrt(m)} (1 + O(m^{-1/2})),
i.e. e^{a sqrt(m)} with a = b/2 — the half-step germ dressing, exact.
NORMALIZATION NOTE: E[Z^aW^a] = 2^a a! forces w = sqrt(2s) in Lambda; germs
of P_b are v = b +- 2 sqrt(2s) and the P^2 germ below is v = 8s. (The
stacking-witness script used w = sqrt(s): its Lambda is the Lambda of a
coefficient-rescaled P; all structural conclusions are unaffected.)

All sums are computed by downward recurrence from the balanced term
(k = floor(m/2)), normalized so the b = 0 even sum is exactly 1:
    T_{k_max} = b^{m-2k_max};  T_{k-1} = T_k * b^2 * k / (2 (m-2k+1)(m-2k+2))
    E[P_b^m] = (m! 2^{k_max} / k_max!) * sum_k T_k.

TESTS
(B1) grading law: log(sum_b(m))/sqrt(m) -> b/2 (complex: Re = grade,
     Im = the sqrt(m)-drifting phase the old 1/m ladder missed).
(B2) parity/vacuity (referee obligation (3) instance): b = 0 kills ALL odd
     moments exactly; b != 0 restores them.
(B3) TIE RECOVERY: S_m = E-hat[P_b1^m] + e^{i phi} E-hat[P_b2^m], two germs
     tied at leading order. Grade -> b1; deflate -> b2. Phase unwrapped
     along the even-m grid.
(B4) CONJUGATE PAIR: T_m = E-hat[P_b^m] + E-hat[P_conj b^m] real; envelope
     slope -> Re b/2, sign-crossing spacing in sqrt(m) -> |Im b|/2 (the van
     der Corput stage made visible).
(B5) mirror-merge -> P^2 rotation stack: E[P_0^{2k}]/(8^k k!) = C(2k,k)/4^k
     ~ 1/sqrt(pi k): the v = 8s germ + the exact -1/2 prefactor rung.

boxeph-2026-07-20-S183. Pure python, no deps.
"""
import math
import cmath


def Ehat(b, m):
    """sum_k T_k (normalized moment). Even m, b=0 -> exactly 1."""
    kmax = m // 2
    T = b ** (m - 2 * kmax)
    tot = T
    k = kmax
    while k >= 1:
        T = T * b * b * k / (2.0 * (m - 2 * k + 1) * (m - 2 * k + 2))
        tot += T
        k -= 1
    return tot


def logE0(m):
    k = m // 2
    return math.lgamma(m + 1) + k * math.log(2) - math.lgamma(k + 1)


print("=" * 78)
print("REPAIRED STAGE B on exact circular-pair moments of P_b = Z + b + W")
print("=" * 78)

print("\n(B1) grading law  log(Ehat_b(m))/sqrt(m) -> b/2  (even m):")
for b in (complex(0.8, 0.0), complex(0.5, 0.6), complex(-0.3, 1.0)):
    line = []
    # unwrapped log along the even grid
    ph_prev, wraps = None, 0
    for m in range(16, 601, 2):
        z = Ehat(b, m)
        ph = cmath.phase(z)
        if ph_prev is not None:
            while ph + 2 * math.pi * wraps - ph_prev > math.pi:
                wraps -= 1
            while ph + 2 * math.pi * wraps - ph_prev < -math.pi:
                wraps += 1
        ph_prev = ph + 2 * math.pi * wraps
        if m in (64, 144, 256, 400, 576):
            est = complex(math.log(abs(z)), ph_prev) / math.sqrt(m)
            line.append("m=%d: %.4f%+.4fj" % (m, est.real, est.imag))
    print("  b=%s (target %.4f%+.4fj): %s" % (b, b.real / 2, b.imag / 2, "; ".join(line)))

print("\n(B2) parity vacuity: |E-hat[P_b^m]| for odd m:")
for b in (0.0, 0.4):
    vals = [abs(Ehat(complex(b, 0), m)) for m in (11, 21, 31)]
    print("  b=%.1f : %s -> odd conditions %s" %
          (b, ["%.3e" % v for v in vals],
           "VACUOUS (parity family)" if max(vals) == 0 else "NONVACUOUS"))

print("\n(B3) tie recovery: S_m = Ehat_b1 + e^{i phi} Ehat_b2 (tied leading order):")
b1 = complex(0.9, 0.5)
b2 = complex(0.3, -0.7)
phi = 1.1
print("    b1=%s (Re a1=%.3f), b2=%s (Re a2=%.3f), phi=%.1f" %
      (b1, b1.real / 2, b2, b2.real / 2, phi))
S = {m: Ehat(b1, m) + cmath.exp(1j * phi) * Ehat(b2, m) for m in range(16, 801, 2)}


def graded_fit(seq_fn, mlist, mgrid):
    ph_prev, wraps, logs = None, 0, {}
    for m in mgrid:
        z = seq_fn(m)
        ph = cmath.phase(z)
        if ph_prev is not None:
            while ph + 2 * math.pi * wraps - ph_prev > math.pi:
                wraps -= 1
            while ph + 2 * math.pi * wraps - ph_prev < -math.pi:
                wraps += 1
        ph_prev = ph + 2 * math.pi * wraps
        logs[m] = complex(math.log(abs(z)), ph_prev)
    out = []
    for m in mlist:
        out.append((m, logs[m] / math.sqrt(m)))
    m1, m2 = mlist[-2], mlist[-1]
    two_pt = (logs[m2] - logs[m1]) / (math.sqrt(m2) - math.sqrt(m1))
    return out, two_pt


grid = list(range(16, 801, 2))
tab, tp = graded_fit(lambda m: S[m], [144, 256, 400, 576, 800], grid)
for m, est in tab:
    print("      m=%3d  a-hat = %.5f%+.5fj  (target b1/2 = %.5f%+.5fj)" %
          (m, est.real, est.imag, b1.real / 2, b1.imag / 2))
print("    two-point a-hat = %.5f%+.5fj -> b1-hat = %.4f%+.4fj (true %s)" %
      (tp.real, tp.imag, (2 * tp).real, (2 * tp).imag, b1))
print("    deflate (exact b1 rung), re-grade the remainder:")
S2 = {m: S[m] - Ehat(b1, m) for m in grid}
tab2, tp2 = graded_fit(lambda m: S2[m] / cmath.exp(1j * phi), [256, 400, 576, 800], grid)
for m, est in tab2:
    print("      m=%3d  a-hat = %.5f%+.5fj  (target b2/2 = %.5f%+.5fj)" %
          (m, est.real, est.imag, b2.real / 2, b2.imag / 2))
print("    two-point a-hat = %.5f%+.5fj -> b2-hat = %.4f%+.4fj (true %s)" %
      (tp2.real, tp2.imag, (2 * tp2).real, (2 * tp2).imag, b2))

print("\n(B4) conjugate pair T_m = Ehat_b + Ehat_conj(b) (REAL data):")
bc = complex(0.6, 0.8)
print("    b=%s: targets Re b/2 = %.3f, |Im b|/2 = %.3f" %
      (bc, bc.real / 2, abs(bc.imag) / 2))
T = {}
for m in range(16, 2001, 2):
    T[m] = (Ehat(bc, m) + Ehat(bc.conjugate(), m)).real
lastsign, crossings = 0, []
for m in range(16, 2001, 2):
    sgn = 1 if T[m] > 0 else -1
    if lastsign and sgn != lastsign:
        crossings.append(math.sqrt(m))
    lastsign = sgn
gaps = [crossings[i + 1] - crossings[i] for i in range(len(crossings) - 1)]
avg_gap = sum(gaps[-8:]) / len(gaps[-8:])
print("    sign-crossing spacing in sqrt(m) = %.4f -> extracted |Im b|/2 = %.4f" %
      (avg_gap, math.pi / avg_gap))
win = {}
for m in range(16, 2001, 2):
    x = math.sqrt(m)
    v = abs(T[m])
    key = int(x)
    if key not in win or v > win[key][1]:
        win[key] = (x, v)
pts = sorted(win.values())
x1, v1 = pts[len(pts) // 2]
x2, v2 = pts[-2]
print("    envelope log-slope = %.4f (target Re b/2 = %.4f)" %
      ((math.log(v2) - math.log(v1)) / (x2 - x1), bc.real / 2))

print("\n(B5) mirror-merge -> P^2 rotation stack (v = 8s germ + the -1/2 rung):")
print("    R_k = E[P_0^{2k}]/(8^k k!) = C(2k,k)/4^k; check R_k*sqrt(pi k) -> 1:")
for k in (10, 40, 160, 640):
    lR = logE0(2 * k) - k * math.log(8) - math.lgamma(k + 1)
    val = math.exp(lR) * math.sqrt(math.pi * k)
    print("      k=%4d   R_k*sqrt(pi k) = %.6f" % (k, val))

print("\nDONE.")
