#!/usr/bin/env python3
"""
HYP-6994 RESONANCE STRESS TEST (boxeph-2026-07-16-S25)

Question: is the UNIFORM sup-norm lemma max_{n!=0}|S(n)|^2 = O(M) (HYP-6994,
klein THM-881/assault) actually TRUE asymptotically on the scanned family
E_t = {1..6, t}, or is the measured C <= 14 a small-t artifact?

Prediction (owner-resonance mechanism): at n = t the big owner's entry/exit
sums have zeta = e(n/t) = 1 (no intra-owner cancellation), leaving
  |S_t(t)| = | e(1/7)|Y_t| - |X_t| | + O(sum of small owners)
           >= sin(2pi/7)|Y_t| - O(1) ~ 0.78 * d * t   (grows linearly),
where X_t, Y_t are the ACTIVE entry/exit index sets of owner t (t-grid samples
of the {1..6}-frame's swing condition, density d > 0). Hence
max|S|^2 / M should grow ~ c * t while the WEIGHTED form Q_s stays O(M)
(the resonance sits under khat(t) ~ 1/t^2 weight).

Outputs, per t:
  M, max_{n != 0} |S(n)|^2 / M (full scan when P manageable),
  |S(q t)|^2 / M for q = 1..3 (the predicted resonance),
  Q_s(w)/M for a few coprime w via THM-880's exact bilinear form.

R_s convention exactly as klein's scripts: R_s = {x : the 7e-section occupancy
of E has len(occ) == 6 and s not in occ}. Signs: +1 entering R_s, -1 leaving.

Pure Python.
"""

from fractions import Fraction as Fr
from math import gcd, lcm, pi, sin, cos
import cmath, sys

def endpoints(E, s):
    """R_s endpoints as (position Fr, sign, owner). Position in [0,1)."""
    bps = sorted(set(Fr(k, 7 * e) for e in E for k in range(7 * e)) | {Fr(0), Fr(1)})
    pts = []
    prev_in = None
    # evaluate occupancy on each open cell
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        occ = set(int((e * mid % 1) * 7) for e in E)
        cur = (len(occ) == 6) and (s not in occ)
        if prev_in is None:
            first_in = cur
        else:
            if cur and not prev_in:
                pts.append([bps[i], +1])
            elif prev_in and not cur:
                pts.append([bps[i], -1])
        prev_in = cur
    # wrap at 0/1
    if prev_in != first_in:
        if first_in and not prev_in:
            pts.append([Fr(0), +1])
        else:
            pts.append([Fr(0), -1])
    # owner attribution: the runner whose section boundary sits at p
    out = []
    for p, sg in pts:
        own = [e for e in E if (p * 7 * e).denominator == 1]
        out.append((p, sg, min(own)))
    return out

def spectrum_scan(E, s, nmax=None, extra=()):
    """max_{n!=0}|S(n)|^2 over Z_P (full scan if nmax is None), plus values
    at requested extra frequencies. Returns (M, best_n, best_val, {n: |S|^2})."""
    P = 7 * lcm(*E)
    pts = endpoints(E, s)
    M = len(pts)
    pos = [int(p * P) for p, sg, o in pts]
    sgn = [sg for p, sg, o in pts]
    vals = {}
    best = (0, -1.0)
    rng = range(1, P) if nmax is None else range(1, min(nmax, P))
    for n in rng:
        z = sum(sg * cmath.exp(2j * pi * n * q / P) for sg, q in zip(sgn, pos))
        v = abs(z) ** 2
        if v > best[1]:
            best = (n, v)
    for n in extra:
        n = n % P
        if n == 0:
            continue
        z = sum(sg * cmath.exp(2j * pi * n * q / P) for sg, q in zip(sgn, pos))
        vals[n] = abs(z) ** 2
        if vals[n] > best[1]:
            best = (n, vals[n])
    return M, best[0], best[1], vals, P

def Qs_bilinear(E, s, w):
    """THM-880 exact form: Q_s(w) = -2 pi^2 sum_{k,k'} e_k e_k' {w D}(1-{w D})."""
    pts = endpoints(E, s)
    tot = 0.0
    for pk, sk, _ in pts:
        for pj, sj, _ in pts:
            d = (w * (pk - pj)) % 1
            f = float(d)
            tot += sk * sj * f * (1 - f)
    return -2 * pi * pi * tot

if __name__ == "__main__":
    print("family E_t = {1..6, t}; worst section reported")
    print(f"{'t':>5} {'M':>5} {'P':>7} | {'max|S|^2/M':>11} {'argmax':>7} | "
          f"{'|S(t)|^2/M':>10} {'|S(2t)|^2/M':>11} | {'Qs(1)/M':>8} {'Qs(w*)/M':>9}")
    for t in [6, 12, 25, 37, 50, 60, 120, 300, 600, 1200]:
        E = [1, 2, 3, 4, 5, 6, t]
        wrec = None
        rows = []
        for s in range(7):
            P = 7 * lcm(*E)
            full = P <= 25000
            M, bn, bv, vals, P = spectrum_scan(
                E, s, nmax=None if full else 3000,
                extra=[t, 2 * t, 3 * t, P - t])
            if M == 0:
                continue
            # a coprime w for Q_s
            w = 11
            while gcd(w, P) != 1:
                w += 2
            q1 = Qs_bilinear(E, s, 1)
            qw = Qs_bilinear(E, s, w)
            rows.append((s, M, P, bv / M, bn, vals.get(t % P, 0) / M,
                         vals.get((2 * t) % P, 0) / M, q1 / M, qw / M, full))
        r = max(rows, key=lambda x: x[3])
        s, M, P, ratio, bn, rt, r2t, q1, qw, full = r
        note = "" if full else " (restricted scan + resonances)"
        print(f"{t:>5} {M:>5} {P:>7} | {ratio:>11.2f} {bn:>7} | "
              f"{rt:>10.2f} {r2t:>11.2f} | {q1:>8.2f} {qw:>9.2f}{note}")

# ---------------------------------------------------------------- PART 2
# Lemma referees for THM-886 (boxeph-S25): the spectrum profile, the mode
# identity, the k-hat bound, and the assembled resonance law of Q_s(w).

def part2():
    from math import lcm, sin, cos
    print()
    print("=" * 74)
    print("PART 2 -- lemma referees (family E_t = {1..6, t})")

    # (L-khat) closed form / bound: |khat_P(n)| <= 1/(4 P^2 sin^2(pi n / P))
    # and khat_P(n) is real negative for n != 0; verify on a few P.
    for P in (420, 840):
        ok = True
        worst = 0.0
        for n in range(1, P):
            kh = sum((j / P) * (1 - j / P) * cmath.exp(-2j * pi * n * j / P)
                     for j in range(P)) / P
            bound = 1 / (4 * P * P * sin(pi * n / P) ** 2)
            if abs(kh) > bound * (1 + 1e-9):
                ok = False
            worst = max(worst, abs(kh) * 16 * min(n, P - n) ** 2)
        print(f"  (L-khat) P={P}: |khat(n)| <= 1/(4P^2 sin^2(pi n/P)): {ok}; "
              f"max |khat|*16*ntilde^2 = {worst:.3f} (<= 1 iff the 1/(16n^2) bound holds)")

    t = 600
    E = [1, 2, 3, 4, 5, 6, t]
    P = 7 * lcm(*E)
    s = 0
    pts = endpoints(E, s)
    M = len(pts)
    pos = [int(p * P) for p, sg, o in pts]
    sgn = [sg for p, sg, o in pts]
    own = [o for p, sg, o in pts]

    def S(n):
        return sum(sg * cmath.exp(2j * pi * (n % P) * q / P)
                   for sg, q in zip(sgn, pos))

    # frame constants: C_f = # frame-owned endpoints. The t-owned active
    # boundaries j/(7t) carry ALL residues c = j mod 7 (the exactly-six
    # occupancy condition reacts to t vacating/entering ANY section): per-class
    # signed counts N_c and per-class index sets in Z_t with run counts R_c.
    Cf = sum(1 for o in own if o < t)
    N = [0] * 7
    idx = {c: {+1: [], -1: []} for c in range(7)}
    for q, sg, o in zip(pos, sgn, own):
        if o == t:
            j = q * 7 * t // P          # position = j/(7t)
            c = j % 7
            N[c] += sg
            idx[c][sg].append(j // 7)   # period index m in Z_t

    def runs(lst, mod):
        st = set(lst)
        return sum(1 for m in st if (m - 1) % mod not in st)

    R = sum(runs(idx[c][sg], t) for c in range(7) for sg in (+1, -1))
    Mt = sum(1 for o in own if o == t)
    print(f"  frame endpoints C_f = {Cf}; t-owned M_t = {Mt}; "
          f"signed residue counts N = {N}; total runs R = {R}")

    # (L-mode) EXACT 7-residue mode identity:
    #   S(ta) = S_frame(ta) + sum_c e(ac/7) * N_c   for every integer a
    okm = True
    for a in (1, 2, 3, 4, 5, 6, 7, 11):
        Sf = sum(sg * cmath.exp(2j * pi * (a * t % P) * q / P)
                 for sg, q, o in zip(sgn, pos, own) if o < t)
        pred = Sf + sum(N[c] * cmath.exp(2j * pi * a * c / 7) for c in range(7))
        err = abs(S(a * t) - pred)
        okm = okm and err < 1e-8
    print(f"  (L-mode) S(ta) = S_frame(ta) + sum_c e(ac/7) N_c EXACT (a=1..7,11): {okm}")
    nu = [abs(sum(N[c] * cmath.exp(2j * pi * a * c / 7) for c in range(7)))
          for a in range(1, 7)]
    print(f"  mode amplitudes t|nu^(a)|-part: {[f'{v:.1f}' for v in nu]} "
          f"(vs measured |S(ta)|: {[f'{abs(S(a*t)):.1f}' for a in range(1,7)]})")

    # (L-profile) run bound off-mode: |S(m)| <= C_f + |sum over X runs| + |Y runs|
    # <= C_f + (R_X + R_Y) / (2 ||m/t||); verify pointwise over all m.
    viol = 0
    worst_ratio = 0.0
    for m in range(1, P):
        d = abs((m / t) - round(m / t))
        if d == 0:
            continue
        bound = Cf + R / (2 * sin(pi * min(d, 1 - d)))
        v = abs(S(m))
        if v > bound + 1e-6:
            viol += 1
        worst_ratio = max(worst_ratio, v / bound)
    print(f"  (L-profile) |S(m)| <= C_f + R/(2 sin(pi ||m/t||)) off-mode: "
          f"violations = {viol}, worst ratio = {worst_ratio:.3f}")

    # (L-low) low band: |S(n)| <= 2 pi n lambda(R_s) + ... use |S(n)| <= min(M, 2 pi n lam)
    lam = sum(float(b - a) for a, b in [])  # placeholder not needed; compute:
    # reconstruct R_s measure from endpoints: sum over arcs
    arcs = []
    cur = None
    for p, sg, o in sorted(pts):
        if sg == +1:
            cur = p
        elif cur is not None:
            arcs.append((cur, p)); cur = None
    lam = float(sum(b - a for a, b in arcs))
    okl = True
    for n in range(1, 200):
        if abs(S(n)) > min(M, 2 * pi * n * lam) + 1e-6:
            okl = False
    print(f"  (L-low) |S(n)| <= min(M, 2 pi n lambda) for n < 200: {okl} (lambda = {lam:.4f})")

    # (VALIDATE) the assembled law: Q_s(w) vs C1*M + resonance sum over hits
    # hits: (l, a, rho): l*w == a*t + rho (mod P), l <= 30, |rho| <= t/2, a != 0 mod 7*?
    print("  assembled resonance law check (upper-bound shape):")
    for w in (11, 1013, 601, 1207, 105):
        q = Qs_bilinear(E, s, w)
        res = 0.0
        for l in range(1, 31):
            m = (l * w) % P
            a = round(m / t)
            rho = m - a * t
            if a % 7 == 0 and a != 0:
                amp = Cf + abs(sum(N))
            elif a == 0:
                continue
            elif rho == 0:
                amp = abs(sum(N[c] * cmath.exp(2j * pi * a * c / 7) for c in range(7))) + Cf
            else:
                amp = Cf + R / (2 * sin(pi * min(abs(rho) / t, 1 - abs(rho) / t)))
            res += (1 / (16 * l * l)) * amp * amp * 2 * pi * pi * 2
        print(f"    w={w:>5}: Q_s = {q:9.1f};  C1*M + mode-hit sum = {12*M + res:9.1f}  "
              f"(cover: {q <= 12*M + res + 1e-6})")

part2()
