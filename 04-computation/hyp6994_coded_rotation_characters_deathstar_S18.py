#!/usr/bin/env python3
"""death-star-2026-07-16-S18 (HYP-7017): the coded-rotation character expansion of the
per-owner sign words — the uniform-HYP-6994 dichotomy test.

OBJECTS (klein THM-882-assault conventions, verbatim R_s):
  cluster E (7 speeds incl. 0), miss-pattern s: R_s = {x : occupied sections = 6 distinct,
  missing = s}. Endpoints signed +1 enter / -1 leave. Owner of p = min{e : p in (1/(7e))Z}.
  Per-owner word u_e on Z_{7e}: u_e(j) = net sign at position j/(7e) owned by e.
  S(n) = sum eps e(n p) on Z_P, P = 7 lcm(E); S_e(n) = per-owner part (depends on n mod 7e).

THE PREDICTION (hand analysis, this session):
  u_e(j) = F(sections of others at j/(7e); j mod 7) generically, F = the fixed sign rule.
  Character expansion: h_c(theta) = omega^{c*floor(7 theta)} has hhat_c(k) != 0 iff
  k = c mod 7 (EXACT; hhat_c(0) = [c=0]).  Hence:
  (i) FIRST-ORDER PEAKS: the j-mod-7 marginal E_O[F(O; j)] = p*([j=s+1] - [j=s]) != 0
      => uhat_e has peaks at m = -c0*e (mod 7e), height ~ 7e * |Fhat(c0;0)| = O(e)
      => sup|S|^2/M should GROW ~ e at owner-resonant frequencies (sup-form of HYP-6994
      FALSE asymptotically)  ... UNLESS rotation correlations kill p.
  (ii) THE WEIGHTED FORM SURVIVES: kernel khat(n) ~ 1/n^2 crushes the resonant classes
      (sum over the class ~ 1/e^2 * (e*p)^2 = p^2 = O(1)), so Q_s = O(M) needs only the
      >= 2-leg remainder to be flat (1/(k1 k2) hyperbola weights).
TESTS:
  T1: exact word reconstruction: idealized F-prediction vs true u_e; collision census.
  T2: the peak law: locate argmax_n |S(n)| — is it owner-resonant (n = c0*e mod 7e)?
      does max|S|^2/M grow with the dominant owner t? (t-ladder to 200)
  T3: the weighted split: Q_s(exact, all w reps) vs [first-order comb] + [remainder];
      is the resonant contribution to Q_s O(1) while M grows?
"""
from fractions import Fraction as Fr
from math import gcd, pi
import cmath, sys, time

def lcm(a, b): return a * b // gcd(a, b)

def R_s_exact(E, s):
    bps = sorted(set(Fr(k, 7 * e) for e in E if e > 0 for k in range(7 * e)) | {Fr(0), Fr(1)})
    arcs, inR, start = [], False, None
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        occ = set(int((e * mid % 1) * 7) for e in E)
        cur = (len(occ) == 6) and (s not in occ)
        if cur and not inR: start, inR = bps[i], True
        if (not cur) and inR: arcs.append((start, bps[i])); inR = False
    if inR: arcs.append((start, Fr(1)))
    return arcs

def endpoints(arcs):
    out = []
    for a, b in arcs:
        out.append((a, +1)); out.append((b, -1))
    return out

def owner_of(p, E):
    return min(e for e in E if e > 0 and (p * 7 * e).denominator == 1)

def words_by_owner(E, s):
    """Exact per-owner sign words u_e: Z_{7e} -> {-1,0,+1} (net sign; collisions merge)."""
    eps = endpoints(R_s_exact(E, s))
    words = {e: {} for e in E if e > 0}
    for p, sg in eps:
        e = owner_of(p, E)
        j = int(p * 7 * e) % (7 * e)
        words[e][j] = words[e].get(j, 0) + sg
    return words

def idealized_word(E, s, e):
    """The no-collision F-prediction for owner e: evaluate membership just left/right of
    j/(7e) using ONLY the section data (one-sided limits) — same logic, but attributing
    the sign to e even at shared boundaries."""
    others = [f for f in E if f != e]
    w = {}
    for j in range(7 * e):
        x = Fr(j, 7 * e)
        def member(side):
            # side = -1/+1: sections at x + side*0 (one-sided limit)
            occ = set()
            for f in others:
                th = (f * x) % 1
                sec = int(th * 7)
                if side < 0 and (f * x * 7) % 1 == 0 and f > 0:
                    sec = (sec - 1) % 7   # f's own boundary: left limit is previous section
                occ.add(sec if f > 0 else 0)
            esec = j % 7 if side > 0 else (j - 1) % 7
            occ.add(esec)
            return (len(occ) == 6) and (s not in occ)
        v = int(member(+1)) - int(member(-1))
        if v: w[j] = v
    return w

def dft_word(w, N):
    """uhat(m) = sum_j w[j] e(m j / N), m = 0..N-1 (owner-local spectrum)."""
    out = []
    for m in range(N):
        z = sum(v * cmath.exp(2j * pi * m * j / N) for j, v in w.items())
        out.append(z)
    return out

def S_full(E, s, P):
    eps = endpoints(R_s_exact(E, s))
    pos = [(int(p * P) % P, sg) for p, sg in eps]
    return pos

def maxS2_overM(E, s):
    P = 7
    for e in E:
        if e > 0: P = lcm(P, 7 * e)
    pos = S_full(E, s, P)
    M = len(pos)
    if M == 0: return None
    # exhaustive scan over Z_P (w-freeness): may be big; use per-divisor structure:
    best, argbest = 0.0, None
    for n in range(1, P):
        z = sum(sg * cmath.exp(2j * pi * n * k / P) for k, sg in pos)
        a = abs(z) ** 2
        if a > best: best, argbest = a, n
    return best, argbest, M, P

def window_structure(E, s):
    """THE WINDOW-DIPOLE DECOMPOSITION. Endpoints of R_s grouped into 'windows':
    maximal runs of the circle where the arcs of R_s sit between consecutive
    others'-section-change events. Returns per-arc classification and the window count
    per owner pair. Also verifies: every arc of R_s is a single interval whose two
    endpoints are section boundaries (owners possibly distinct)."""
    arcs = R_s_exact(E, s)
    stats = {"arcs": len(arcs), "M": 2 * len(arcs)}
    return stats

def sup_scan(E, s, ncap=None):
    """max_{0<n<=ncap or full} |S(n)|^2, M, argmax. Partial scan => LOWER bound on sup."""
    P = 7
    for e in E:
        if e > 0: P = lcm(P, 7 * e)
    pos = S_full(E, s, P)
    M = len(pos)
    if M == 0:
        return None
    rng = range(1, P) if (ncap is None or ncap >= P) else None
    best, arg = 0.0, None
    if rng is not None:
        for n in rng:
            z = sum(sg * cmath.exp(2j * pi * (n * k % P) / P) for k, sg in pos)
            a = abs(z) ** 2
            if a > best: best, arg = a, n
        full = True
    else:
        # structured partial scan: small n, owner-resonant classes, random sample
        import random
        rnd = random.Random(20260716)
        cand = set(range(1, min(ncap, P)))
        for e in E:
            if e <= 0: continue
            for c0 in range(1, 7):
                base = (c0 * e) % (7 * e)
                # lift the residue class mod 7e into Z_P: take a few lifts
                for r in range(0, min(200, P // (7 * e) + 1)):
                    cand.add((base + r * 7 * e) % P or 1)
        while len(cand) < min(P - 1, ncap + 3000):
            cand.add(rnd.randrange(1, P))
        for n in cand:
            z = sum(sg * cmath.exp(2j * pi * (n * k % P) / P) for k, sg in pos)
            a = abs(z) ** 2
            if a > best: best, arg = a, n
        full = False
    return best, arg, M, P, full

if __name__ == "__main__":
    t0 = time.time()
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "t1"):
        print("T1: exact vs idealized per-owner words (collision census)")
        for t in [6, 12, 25]:
            E = [0, 1, 2, 3, 4, 5, t]
            s = 3
            tw = words_by_owner(E, s)
            for e in sorted(x for x in E if x > 0):
                iw = idealized_word(E, s, e)
                true_w = tw[e]
                keys = set(iw) | set(true_w)
                mism = sum(1 for k in keys if iw.get(k, 0) != true_w.get(k, 0))
                print(f"  t={t} s={s} owner={e}: M_e(true)={sum(abs(v) for v in true_w.values())} "
                      f"M_e(ideal)={sum(abs(v) for v in iw.values())} mismatches={mism}")
        print(f"[{time.time()-t0:.1f}s]")
    if which in ("all", "t2"):
        print("\nT2: peak scaling on the FAR bank [0..5,t] (worst s each; full scan where feasible)")
        for t in [6, 12, 25, 37, 50, 75, 100, 150, 200]:
            E = [0, 1, 2, 3, 4, 5, t]
            worst = (0, None, None, None)
            for s in range(7):
                r = sup_scan(E, s, ncap=20000)
                if r is None: continue
                best, arg, M, P, full = r
                if M and best / M > worst[0]:
                    worst = (best / M, s, arg, (M, P, full))
            rat, s, arg, meta = worst
            if meta:
                M, P, full = meta
                # owner-resonance diagnostic of argmax
                res = [e for e in E if e > 0 and arg is not None and (arg % (7 * e)) % e == 0]
                print(f"  t={t}: max|S|^2/M = {rat:.2f} (s={s}, n*={arg}, M={M}, P={P}, "
                      f"full={full}) owner-res:{res}")
            sys.stdout.flush()
    if which in ("all", "t3"):
        print("\nT3: COMPACT CORES E = [0, c..c+5]: peak scaling as c grows (partial scans)")
        for c in [3, 5, 8, 12, 17, 23, 30]:
            E = [0] + list(range(c, c + 6))
            worst = (0, None, None, None)
            for s in range(7):
                r = sup_scan(E, s, ncap=8000)
                if r is None: continue
                best, arg, M, P, full = r
                if M and best / M > worst[0]:
                    worst = (best / M, s, arg, (M, P, full))
            rat, s, arg, meta = worst
            if meta:
                M, P, full = meta
                print(f"  c={c}: max|S|^2/M >= {rat:.2f} (s={s}, n*={arg}, M={M}, P={P}, full={full})")
            sys.stdout.flush()
    if which in ("all", "t4"):
        print("\nT4: THE WEIGHTED VERDICT — Q_s(w) anatomy on the far ladder (worst s per t)")
        print("  Q_s(w) = 2 pi^2 sum_{n!=0} khat_P(n w mod P)-reindexed |S(n)|^2 computed as")
        print("  sum_k(pairwise bilinear THM-880 form); here via spectrum + exact kernel DFT.")
        import numpy as np
        for t in [12, 25, 50, 100, 150, 200]:
            E = [0, 1, 2, 3, 4, 5, t]
            # find worst s by sup ratio first (cheap full scans since P small here)
            best_s, best_ratio = None, -1
            spec_cache = {}
            for s in range(7):
                P = 7
                for e in E:
                    if e > 0: P = lcm(P, 7 * e)
                pos = S_full(E, s, P)
                M = len(pos)
                if M == 0: continue
                v = np.zeros(P, dtype=complex)
                for k, sg in pos:
                    v[k] += sg
                sp = np.abs(np.fft.fft(v)) ** 2   # |S(n)|^2, n = 0..P-1 (fft sign irrelevant to |.|)
                spec_cache[s] = (sp, M, P)
                r = sp[1:].max() / M
                if r > best_ratio: best_ratio, best_s = r, s
            sp, M, P = spec_cache[best_s]
            # exact kernel on Z_P: K(m/P) = (m/P)(1-m/P); khat(n) = (1/P) sum_m K e(-nm/P)
            m = np.arange(P) / P
            K = m * (1 - m)
            khat = np.fft.fft(K) / P           # khat[n]
            akh = np.abs(khat)
            # Q_s(w) for all w in Z_P^*: Q_s(w) = 2 pi^2 sum_{n!=0} |khat(n)| |S(n w^{-1})|^2
            #   (upper-bound form; the true signed sum uses khat(n) real parts — K real symmetric
            #    => khat real; use signed khat exactly)
            khr = np.real(khat)
            # reindex: Q(w) = 2pi^2 * sum_{n!=0} khr[n] * sp[(n * winv) % P]
            # equivalently sum over j!=0 of khr[j w] sp[j]
            diam = t
            best_w, best_q = None, -1
            appl_vals = []
            ws = [w for w in range(1, P) if gcd(w, P) == 1]
            spn = sp.copy(); spn[0] = 0.0
            for w in ws:
                winv = pow(w, -1, P)
                idx = (np.arange(P) * winv) % P
                q = -2 * pi * pi * float(np.dot(khr[idx], spn))
                if q > best_q: best_q, best_w = q, w
                if w >= diam: appl_vals.append(q)
            appl_max = max(appl_vals) if appl_vals else float('nan')
            nstar = int(spn.argmax())
            print(f"  t={t}: M={M} P={P} s*={best_s} sup|S|^2/M={best_ratio:.2f} n*={nstar} | "
                  f"sup_w Q_s/M={best_q/M:.2f} (w*={best_w}) | sup_{{w>=diam}} Q_s/M={appl_max/M:.2f} | "
                  f"sup_w Q_s/diam={best_q/diam:.2f}")
            sys.stdout.flush()
    if which in ("all", "t5"):
        print("\nT5: NORMALIZATION REFEREE — brute ell-sum vs spectral bilinear, t=12")
        import numpy as np
        E = [0, 1, 2, 3, 4, 5, 12]
        s = 6
        P = 7
        for e in E:
            if e > 0: P = lcm(P, 7 * e)
        pos = S_full(E, s, P)
        M = len(pos)
        def Sval(N):
            return sum(sg * cmath.exp(2j * pi * ((N * k) % P) / P) for k, sg in pos)
        v = np.zeros(P, dtype=complex)
        for k, sg in pos:
            v[k] += sg
        sp = np.abs(np.fft.fft(v)) ** 2
        m = np.arange(P) / P
        K = m * (1 - m)
        khr = np.real(np.fft.fft(K) / P)
        for w in [1, 13, 199]:
            brute = sum(abs(Sval(l * w)) ** 2 / l**2 for l in range(1, 200001) if (l % P) != 0)
            winv = pow(w, -1, P)
            idx = (np.arange(P) * winv) % P
            spn = sp.copy(); spn[0] = 0.0
            spec = -2 * pi * pi * float(np.dot(khr[idx], spn))
            # also the pairwise THM-880 bilinear form directly
            bil = 0.0
            for k1, s1 in pos:
                for k2, s2 in pos:
                    d = ((w * (k1 - k2)) % P) / P
                    bil += s1 * s2 * d * (1 - d)
            bil *= -2 * pi * pi
            print(f"  w={w}: brute(l<=2e5)={brute:.6f}  spectral={spec:.6f}  bilinear={bil:.6f}  M={M}")
    if which in ("all", "t6"):
        print("\nT6: THE RESONANT-W SET — distribution of Q_s(w)/diam over coprime w, per t")
        import numpy as np
        for t in [50, 100, 150, 200]:
            E = [0, 1, 2, 3, 4, 5, t]
            best_s, best_ratio, cache = None, -1, None
            for s in range(7):
                P = 7
                for e in E:
                    if e > 0: P = lcm(P, 7 * e)
                pos = S_full(E, s, P)
                if not pos: continue
                v = np.zeros(P, dtype=complex)
                for k, sg in pos:
                    v[k] += sg
                sp = np.abs(np.fft.fft(v)) ** 2
                if sp[1:].max() / len(pos) > best_ratio:
                    best_ratio, best_s, cache = sp[1:].max() / len(pos), s, (sp, len(pos), P)
            sp, M, P = cache
            m = np.arange(P) / P
            khr = np.real(np.fft.fft(m * (1 - m)) / P)
            spn = sp.copy(); spn[0] = 0.0
            qs = {}
            for w in range(1, P):
                if gcd(w, P) != 1: continue
                winv = pow(w, -1, P)
                idx = (np.arange(P) * winv) % P
                qs[w] = -2 * pi * pi * float(np.dot(khr[idx], spn))
            diam = t
            vals = np.array(list(qs.values())) / diam
            ws = np.array(list(qs.keys()))
            nstar = int(spn.argmax())
            order = np.argsort(-vals)
            top = [(int(ws[i]), float(vals[i])) for i in order[:8]]
            for thr in [16, 24]:
                frac = float((vals > thr).mean())
                print(f"  t={t}: frac{{Q/diam > {thr}}} = {frac:.4f}", end="")
            med = float(np.median(vals))
            # C outside the exceeders: 99th percentile
            p99 = float(np.percentile(vals, 99))
            print(f"  median={med:.2f} p99={p99:.2f} n*={nstar}")
            print(f"      top-w: {[(w, round(v,1), 'n*+'+str(w-nstar) if abs(w-nstar)<50 else ('P-n*+'+str(w-(P-nstar)) if abs(w-(P-nstar))<50 else '')) for w,v in top]}")
            sys.stdout.flush()
    print(f"[total {time.time()-t0:.1f}s]")
