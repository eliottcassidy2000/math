#!/usr/bin/env python3
"""ROTATION-GRAM FORCING TEST mod 13 for the LRC(14) 91-line  (compute-agent scratch referee).

Frame (HYP-9050 + HYP-9055 + THM-2356 + THM-2337):
  grid Z/91;  band B = {x : dist91(x) <= 6}, |B| = 13   (the 1/14 danger band: 6 = floor(91/14));
  13-speed family v;  incidence (i,j): residue x = v_i*j mod 91, in-band iff x in B;
  duty sectors by g = gcd(v_i,91) in {1,7,13,91};  profile z = (z1,z7,z13,z91).

  JET (THM-2337-style first divided difference, balanced gauge):
     <x>   = balanced rep of x mod 91 in [-45,45];   <x>_7 = balanced rep of x mod 7 in [-3,3];
     J(x)  = (<x> - <x>_7)/7  in {-6..6}  (an exact integer divided difference);
     beta  = J mod 13  (the Bockstein-style jet).  On the septimal-zero sector 7|x the section
     drops out:  J(x) = <x>/7,  beta(x) = 2x mod 13  (7^{-1} = 2 mod 13)  -- gauge-free.

  PAIRING (THM-2356 / HYP-9055 rotation-Gram style, per-incidence weights):
     T[w,w'](v) = sum_{i<i'} sum_{j in Z/91}  sym( w(x_i) w'(x_{i'}) );
     per-pair full-orbit sums = weighted Gram entries Gamma_{w,w'}(rho) = sum_x w(x) w'(rho x).

  Rotation covariance (precise): the weight of incidence (i,j) must factor through the residue
  x = v_i*j mod 91 via ONE fixed function (no dependence on i, on j, or on the integer lift).
"""

import random
from math import gcd, comb
from itertools import combinations

M = 91
FAILS = []


def check(cond, msg):
    tag = "OK " if cond else "FAIL"
    if not cond:
        FAILS.append(msg)
    print(f"  [{tag}] {msg}")


def bal(x, m):
    x %= m
    return x - m if x > (m - 1) // 2 else x


def dist(x, m=M):
    x %= m
    return min(x, m - x)


BR = 6
BAND = frozenset(x % M for x in range(-BR, BR + 1))
COV = [frozenset(j for j in range(M) if (r * j) % M in BAND) for r in range(M)]


def J(x):        # balanced-gauge divided-difference jet, integer in [-6,6]
    return (bal(x, M) - bal(x, 7)) // 7


def J_std(x):    # standard-gauge jet (section {0..6}) for the gauge demo
    return ((x % M) - (x % 7)) // 7


# ---------------- weights: integer-valued functions of the residue alone ----------------
def w_band(x):        return 1 if x % M in BAND else 0                       # projector (unweighted)
def w_strict_band(x): return (J(x) if x % 7 == 0 else 0) * w_band(x)         # strict jet * band
def w_jet_band(x):    return J(x) * w_band(x)                                # sectioned jet * band
def w_jet2_band(x):   return (J(x) ** 2) * w_band(x)                         # jet-square (chirp phi=x^2/2 flavor)
def w_proj7(x):       return 1 if (x % 7 == 0 and x % M != 0) else 0         # beta^12 Fermat projector = 1_{7Z\0}


def w_beta(k):        # de-banded strict jet moment beta^k on 7Z (F13 value as int 0..12)
    def w(x):
        if x % 7:
            return 0
        return pow(J(x) % 13, k, 13)
    w.__name__ = f"beta^{k}"
    return w


# ---------------- family machinery ----------------
def sectors(v):
    z = {1: 0, 7: 0, 13: 0, 91: 0}
    for vi in v:
        z[gcd(vi, M)] += 1
    return z


def covering(v):
    u = set()
    for vi in v:
        u |= COV[vi % M]
        if len(u) == M:
            return True
    return False


def first_uncovered(v):
    cov = set()
    for vi in v:
        cov |= COV[vi % M]
    for j in range(M):
        if j not in cov:
            return j
    return None


def singleton(v, w):
    """S[w] = sum_i sum_j w(x_i); returns (total, duty_part) exact ints."""
    tot = duty = 0
    for vi in v:
        s = sum(w((vi * j) % M) for j in range(M))
        tot += s
        if gcd(vi, M) > 1:
            duty += s
    return tot, duty


def pair_sums(v, w, wp=None):
    """Unordered pairing sum_{i<i'} sum_j sym(w(x_i) wp(x_i')); returns (total, unit_unit, duty_involving)."""
    same = wp is None
    if same:
        wp = w
    g = [gcd(vi, M) for vi in v]
    tot = uu = 0
    n = len(v)
    for j in range(M):
        xs = [(vi * j) % M for vi in v]
        wa = [w(x) for x in xs]
        wb = wa if same else [wp(x) for x in xs]
        for a, b in combinations(range(n), 2):
            val = wa[a] * wb[b] if same else wa[a] * wb[b] + wb[a] * wa[b]
            tot += val
            if g[a] == 1 and g[b] == 1:
                uu += val
    return tot, uu, tot - uu


def gram(w, wp=None):
    if wp is None:
        wp = w
    return [sum(w(x) * wp((rho * x) % M) for x in range(M)) for rho in range(M)]


def sector_referee(v, w):
    """Exact per-pair orbit-sum classification (same weight, unordered pairs)."""
    G = gram(w)
    ok = True
    W7 = sum(w((7 * t) % M) for t in range(13))
    W13 = sum(w((13 * s) % M) for s in range(7))
    for a, b in combinations(range(len(v)), 2):
        va, vb = v[a] % M, v[b] % M
        ga, gb = gcd(va, M), gcd(vb, M)
        direct = sum(w((va * j) % M) * w((vb * j) % M) for j in range(M))
        if ga == 1:
            pred = G[(vb * pow(va, -1, M)) % M]
        elif gb == 1:
            pred = G[(va * pow(vb, -1, M)) % M]
        elif (ga, gb) == (7, 7):
            s = ((vb // 7) * pow(va // 7, -1, 13)) % 13
            pred = 7 * sum(w((7 * t) % M) * w((7 * ((s * t) % 13)) % M) for t in range(13))
        elif {ga, gb} == {7, 13}:
            pred = W7 * W13
        elif (ga, gb) == (13, 13):
            s = ((vb // 13) * pow(va // 13, -1, 7)) % 7
            pred = 13 * sum(w((13 * t) % M) * w((13 * ((s * t) % 7)) % M) for t in range(7))
        elif {ga, gb} == {7, 91}:
            pred = 7 * W7 * w(0)
        elif {ga, gb} == {13, 91}:
            pred = 13 * W13 * w(0)
        else:
            pred = 91 * w(0) * w(0)
        if direct != pred:
            ok = False
    return ok


# ---------------- family bank ----------------
CANON = {
    "AP {1..13}":           list(range(1, 14)),
    "GW {1..11,13,24}":     list(range(1, 12)) + [13, 24],
    "deepwell {1..12,182}": list(range(1, 13)) + [182],
    "{1..11,13,36}":        list(range(1, 12)) + [13, 36],
}
EXTRA = {
    "loose {2..14}":  list(range(2, 15)),          # S145 control; covering status decided below
    "all-duty {91k}": [91 * k for k in range(1, 14)],
}
UNITS91 = [r for r in range(1, M) if gcd(r, M) == 1]
SEV = [7 * t for t in range(1, 13)]
THIR = [13 * s for s in range(1, 7)]


def sample_family(rng, prof):
    z1, z7, z13, z91 = prof
    res = ([rng.choice(UNITS91) for _ in range(z1)] +
           [rng.choice(SEV) for _ in range(z7)] +
           [rng.choice(THIR) for _ in range(z13)] +
           [0] * z91)
    rng.shuffle(res)
    return [r + M * (idx + 1) for idx, r in enumerate(res)]


def random_bank(seed=20260727):
    rng = random.Random(seed)
    wanted = [(11, 1, 1, 0)] * 5 + [(10, 2, 1, 0)] * 3 + [(10, 1, 2, 0)] * 3 + \
             [(9, 2, 2, 0)] * 2 + [(9, 3, 1, 0)] * 1 + \
             [(10, 1, 1, 1)] * 2 + [(11, 1, 0, 1)] * 2 + [(11, 0, 1, 1)] * 2 + \
             [(12, 0, 0, 1)] * 2 + [(10, 0, 2, 1)] * 1 + [(8, 2, 2, 1)] * 1
    bank = {}
    for prof in wanted:
        for _ in range(4000):
            v = sample_family(rng, prof)
            if covering(v):
                bank[f"rand{len(bank)+1:02d} z={prof}"] = v
                break
    return bank


# =========================================================================================
print("=" * 96)
print("STEP 0: structural facts (family-independent)")
print("=" * 96)
check(sorted(BAND & set(range(0, M, 7))) == [0], "B intersect 7Z = {0}   (no nonzero septimal residue is in the 1/14 band)")
check(sorted(BAND & set(range(0, M, 13))) == [0], "B intersect 13Z = {0}  (same for the tredecimal sector)")
check(J(0) == 0, "J(0) = 0 (a divided-difference jet vanishes at the basepoint)")
check(all(J((-x) % M) == -J(x) for x in range(M)), "J is odd (balanced gauge)")
check(all(J(x) == bal(x, M) // 7 and (J(x) - 2 * x) % 13 == 0 for x in range(0, M, 7)),
      "on 7Z the jet is section-free: J = <x>/7 and beta = 2x mod 13")
check(any((J(x) - J_std(x)) % 13 != 0 for x in range(M) if x % 7),
      "off 7Z the jet is gauge-dependent (balanced vs standard section differ) -- THM-2337 S5 sidecar")
check(all((J(x) - J_std(x)) % 13 == 0 for x in range(0, M, 7)),
      "on 7Z the two sections agree mod 13 (the strict Bockstein is canonical)")
check(all(w_strict_band(x) == 0 for x in range(M)),
      "STRICT jet x band weight is the ZERO function (B cap 7Z = {0}, J(0) = 0): F1 band-septimal collapse")
shell = sorted(bal(x, M) for x in range(M) if w_jet_band(x) != 0)
check(shell == [-6, -5, -4, 4, 5, 6],
      f"sectioned jet x band weight = +-1 exactly on the outer shell +-{{4,5,6}} (all units): {shell}")
check(sum(w_jet_band(x) for x in range(M)) == 0, "sum_B J = 0  (odd jet vs mirror-even band): F2 parity kill")
check(sum(w_jet2_band(x) for x in range(M)) == 6, "sum_B J^2 = 6 (the even repair survives, on the unit shell only)")
ps_ok = all((sum(pow(e, k, 13) for e in range(13)) - (12 if k % 12 == 0 else 0)) % 13 == 0 for k in range(1, 13))
check(ps_ok, "power sums: sum_{e in F13} e^k = 0 mod 13 for k = 1..11, = -1 for k = 12  (F3 orbit-moment kill)")
carry = lambda vi, j: ((vi * j) - ((vi * j) % M)) // M % 13
check((1 * 5) % M == (96 * 1) % M and carry(1, 5) != carry(96, 1),
      "integer-carry weight takes different values on incidences with the SAME residue "
      "-> NOT rotation-covariant; repair = reduce first (J of the residue)")

print()
print("=" * 96)
print("STEP 1: family bank (canonical + controls + random covering 13-speed families)")
print("=" * 96)
BANK = dict(CANON)
BANK.update(EXTRA)
rb = random_bank()
BANK.update(rb)
check(len(rb) >= 20, f"random covering families found: {len(rb)} (target >= 20)")
for name, v in BANK.items():
    z = sectors(v)
    cov = covering(v)
    extra = "" if cov else f"   first uncovered j = {first_uncovered(v)}"
    print(f"  {name:<24} z=({z[1]},{z[7]},{z[13]},{z[91]})  covering={cov}{extra}")
allunits = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 15, 16]
check(not covering(allunits) and first_uncovered(allunits) == 7,
      "all-unit family fails covering exactly at j = 7: covering FORCES z13+z91 >= 1 (j = 13 forces z7+z91 >= 1)")
COVBANK = {n: v for n, v in BANK.items() if covering(v)}
check(all((lambda z: z[13] + z[91] >= 1 and z[7] + z[91] >= 1)(sectors(v)) for v in COVBANK.values()),
      f"every covering family in the bank ({len(COVBANK)}) carries both covering duties (HYP-9050's duty carriers)")

print()
print("=" * 96)
print("STEP 2: sector-calculus referee (exact Z identities: per-pair orbit sums = classified Gram formulas)")
print("=" * 96)
for wname, w in [("1_B", w_band), ("J*1_B", w_jet_band), ("J^2*1_B", w_jet2_band), ("1_{7Z\\0}", w_proj7)]:
    ok = all(sector_referee(v, w) for v in BANK.values())
    check(ok, f"weight {wname}: all per-pair orbit sums match the sector/Gram classification (exact, all families)")

print()
print("=" * 96)
print("STEP 3: Gram spectra Gamma_w(rho) over the 72 unit ratios (mod-13 flatness = localization hypothesis)")
print("=" * 96)
for wname, w in [("1_B (unweighted)", w_band), ("J*1_B (jet)", w_jet_band),
                 ("J^2*1_B (jet-square)", w_jet2_band), ("1_{7Z\\0} (beta^12 proj)", w_proj7),
                 ("beta^1 de-banded", w_beta(1)), ("beta^2 de-banded", w_beta(2)),
                 ("beta^6 de-banded", w_beta(6))]:
    Gm = gram(w)
    vals = sorted({Gm[r] % 13 for r in UNITS91})
    print(f"  {wname:<26} distinct Gamma(rho) mod 13 over units: {vals}   "
          f"{'FLAT' if len(vals) == 1 else 'NON-FLAT'}")
Gm6 = gram(w_beta(6))
check(all(Gm6[r] % 13 == (12 * pow(r, 6, 13)) % 13 for r in UNITS91),
      "beta^6 x beta^6 spectrum: Gamma(rho) = -Legendre(rho mod 13) exactly (extremal non-flat case)")
check(all(gram(w_proj7)[r] % 13 == 12 for r in UNITS91),
      "beta^12 projector spectrum is FLAT = -1 mod 13 (the Fermat degeneration: the only flat escape)")
check(all(gram(w_beta(1))[r] % 13 == 0 for r in UNITS91) and
      all(gram(w_beta(2))[r] % 13 == 0 for r in UNITS91),
      "beta^1, beta^2 de-banded spectra are IDENTICALLY 0 mod 13 (power-sum kill: flat but ZERO)")

print()
print("=" * 96)
print("STEP 4: the pairings, family by family (duty-sector values = the forcing test)")
print("=" * 96)
hdr = (f"  {'family':<24} {'profile':<14} {'NB%13':>5} {'7z7':>4} | {'N2%13':>5} | "
       f"{'Gs':>3} {'Gj':>5} {'Gj%13':>6} {'dtyJ':>4} | {'Gj2%13':>6} {'dtyJ2':>5} | {'P12%13':>6} {'frm':>4}")
print(hdr)
print("  " + "-" * (len(hdr) - 2))
rows = {}
for name, v in BANK.items():
    z = sectors(v)
    prof = (z[1], z[7], z[13], z[91])
    NB, NB_duty = singleton(v, w_band)
    assert NB == 13 * z[1] + 7 * z[7] + 13 * z[13] + 91 * z[91]
    N2 = pair_sums(v, w_band)[0]
    Gs = pair_sums(v, w_strict_band)[0]
    Gj, Gj_uu, Gj_duty = pair_sums(v, w_jet_band)
    Gj2, Gj2_uu, Gj2_duty = pair_sums(v, w_jet2_band)
    P12 = pair_sums(v, w_proj7)[0]
    P12_form = (-comb(z[1], 2) - z[1] * z[7] + 6 * comb(z[7], 2)) % 13
    SJ, SJ_duty = singleton(v, w_jet_band)
    SJ2, SJ2_duty = singleton(v, w_jet2_band)
    assert SJ == 0 and SJ_duty == 0 and SJ2 == 6 * z[1] and SJ2_duty == 0 and Gs == 0, name
    rows[name] = dict(prof=prof, NB=NB, N2=N2, Gj=Gj, Gj_duty=Gj_duty, Gj2=Gj2,
                      Gj2_duty=Gj2_duty, P12=P12, P12_form=P12_form)
    print(f"  {name:<24} {str(prof):<14} {NB % 13:>5} {7 * z[7] % 13:>4} | {N2 % 13:>5} | "
          f"{Gs:>3} {Gj:>5} {Gj % 13:>6} {Gj_duty:>4} | {Gj2 % 13:>6} {Gj2_duty:>5} | {P12 % 13:>6} {P12_form:>4}")

check(all(r['NB'] % 13 == 7 * r['prof'][1] % 13 for r in rows.values()),
      "singleton projector localizes: N_B = 7 z7 mod 13 (13z1 + 13z13 + 91z91 = 0 mod 13)")
check(all(r['Gj_duty'] == 0 and r['Gj2_duty'] == 0 for r in rows.values()),
      "ORIGIN-PINNING: duty-involving part of every jet-weighted band pairing = 0 EXACTLY (in Z), every family")
check(all(r['P12'] % 13 == r['P12_form'] for r in rows.values()),
      "beta^12 projector pairing localizes: P12 = -C(z1,2) - z1*z7 + 6*C(z7,2) mod 13, every family")
check(all(pair_sums(v, w_strict_band)[0] == 0 for v in BANK.values()),
      "STRICT-jet rotation-Gram pairing == 0 identically (dies before any localization question)")

print()
mom_ok = True
for (k, m) in [(1, 1), (1, 2), (2, 2), (3, 4), (5, 6), (6, 7), (11, 11)]:
    for v in BANK.values():
        t = pair_sums(v, w_beta(k), None if m == k else w_beta(m))[0]
        if t % 13 != 0:
            mom_ok = False
check(mom_ok, "de-banded beta^k x beta^m pairings with 12 not| k+m vanish mod 13 for ALL families (F3)")
leg_vals = {n.split()[0]: pair_sums(v, w_beta(6))[0] % 13 for n, v in list(BANK.items())[:6]}
print(f"  beta^6 x beta^6 (Legendre pairing, 12 | k+m, NON-flat) values mod 13: {leg_vals}")

print()
print("=" * 96)
print("STEP 5: localization verdicts (is T mod 13 a function of the duty profile alone?)")
print("=" * 96)
by_prof = {}
for name in COVBANK:
    by_prof.setdefault(rows[name]['prof'], []).append(name)
for label, key in [("N2 (unweighted pair)", lambda r: r['N2'] % 13),
                   ("Gj (jet pair)", lambda r: r['Gj'] % 13),
                   ("Gj2 (jet-square pair)", lambda r: r['Gj2'] % 13),
                   ("P12 (projector pair)", lambda r: r['P12'] % 13)]:
    splits = [(p, {n: key(rows[n]) for n in names}) for p, names in by_prof.items()
              if len({key(rows[n]) for n in names}) > 1]
    if splits:
        prof, ex = max(splits, key=lambda t: len(set(t[1].values())))
        print(f"  {label:<24} NOT localized: profile {prof} gives values {sorted(set(ex.values()))}")
        for n, val in list(ex.items())[:4]:
            print(f"        {n:<28} -> {val}")
    else:
        print(f"  {label:<24} localized on this bank (constant within every duty profile)")

print()
print("=" * 96)
print("STEP 6: band-radius sweep -- where would jet duty content first appear?")
print("=" * 96)
sweep_names = [n for n in ["AP {1..13}", "GW {1..11,13,24}", "deepwell {1..12,182}", "loose {2..14}"] if n in BANK]
for br in [5, 6, 7, 8, 13]:
    Bx = frozenset(x % M for x in range(-br, br + 1))
    wj = lambda x: J(x) * (1 if x % M in Bx else 0)
    wj2 = lambda x: (J(x) ** 2) * (1 if x % M in Bx else 0)
    dj = {n.split()[0]: pair_sums(BANK[n], wj)[2] for n in sweep_names}
    dj2 = {n.split()[0]: pair_sums(BANK[n], wj2)[2] for n in sweep_names}
    print(f"  radius {br:>2}: jet-pair duty part {dj}   jet^2-pair duty part {dj2}")
print("  -> at the LRC(14) radius 6 the duty sector is identically 0.  Jet-square duty content first")
print("     exists at radius 7 = 91/13 (the LRC(13) band); the odd (jet-linear) layer stays mirror-")
print("     suppressed even there.  The 14-runner band is exactly one unit too narrow for any jet duty.")

print()
print("=" * 96)
print("STEP 7: antisymmetric difference-jet pairing (the other 'natural' jet weighting)")
print("=" * 96)


def diff_pair(v, ordered):
    tot = 0
    idx = range(len(v))
    for j in range(M):
        xs = [(vi * j) % M for vi in v]
        pairs = [(a, b) for a in idx for b in idx if a != b] if ordered else list(combinations(idx, 2))
        for a, b in pairs:
            if xs[a] in BAND and xs[b] in BAND and (xs[a] - xs[b]) % 7 == 0:
                tot += J((xs[a] - xs[b]) % M)
    return tot


check(all(diff_pair(v, True) == 0 for v in list(BANK.values())[:5]),
      "ordered difference-jet pairing == 0 identically (antisymmetry of the jet under pair swap)")
check(all(diff_pair(v, False) == 0 for v in list(BANK.values())[:8]),
      "i<i' difference-jet pairing == 0 identically as well: the mirror j -> -j sends each orbit "
      "contribution to its negative (odd jet, even band), and the fixed time j = 0 carries J(0) = 0")
vAP = BANK["AP {1..13}"]
mir_ok = True
for j in range(1, M):
    xs = [(vi * j) % M for vi in vAP]
    ys = [(vi * (M - j)) % M for vi in vAP]
    for a, b in combinations(range(13), 2):
        ca = J((xs[a] - xs[b]) % M) if xs[a] in BAND and xs[b] in BAND and (xs[a] - xs[b]) % 7 == 0 else 0
        cb = J((ys[a] - ys[b]) % M) if ys[a] in BAND and ys[b] in BAND and (ys[a] - ys[b]) % 7 == 0 else 0
        if ca + cb != 0:
            mir_ok = False
check(mir_ok, "mirror kill verified pointwise on AP: contribution(j) + contribution(-j) = 0 for every pair")

print()
print("=" * 96)
print("STEP 8: HYP-9050 13-grid tie-in and exceptional loci of the surviving projector forcings")
print("=" * 96)
ok50 = all(sum(1 for vi in v for j in range(13) if dist(vi * j, 13) <= 1) % 13 ==
           (10 * sum(1 for vi in v if vi % 13 == 0)) % 13 for v in BANK.values())
check(ok50, "HYP-9050 on the 13-grid: N = 10*z13(v) mod 13 for every banked family")
exc_NB, exc_P12, exc_SJ2, tot_prof = [], [], [], 0
for z1 in range(14):
    for z7 in range(14 - z1):
        for z13 in range(14 - z1 - z7):
            z91 = 13 - z1 - z7 - z13
            if z13 + z91 < 1 or z7 + z91 < 1:
                continue
            tot_prof += 1
            if 7 * z7 % 13 == 0:
                exc_NB.append((z1, z7, z13, z91))
            if (-comb(z1, 2) - z1 * z7 + 6 * comb(z7, 2)) % 13 == 0:
                exc_P12.append((z1, z7, z13, z91))
            if 6 * z1 % 13 == 0:
                exc_SJ2.append((z1, z7, z13, z91))
print(f"  realizable covering duty profiles (z sums to 13, both duties present): {tot_prof}")
print(f"  N_B = 7z7 = 0 mod 13 on {len(exc_NB)} profiles (exactly z7 = 0 or 13); e.g. {exc_NB[:4]}")
print(f"  P12 duty value = 0 mod 13 on {len(exc_P12)} profiles; e.g. {exc_P12[:6]}")
print(f"  S[J^2] = 6z1 = 0 mod 13 on {len(exc_SJ2)} profiles (exactly z1 = 0): {exc_SJ2[:3]}")
both = [p for p in exc_NB if p in exc_P12]
print(f"  profiles where BOTH projector forcings die mod 13 simultaneously: {len(both)}; e.g. {both[:6]}")

print()
print("=" * 96)
print("VERDICT")
print("=" * 96)
gj_zero = sum(1 for n in COVBANK if rows[n]['Gj'] % 13 == 0)
gj2_zero = sum(1 for n in COVBANK if rows[n]['Gj2'] % 13 == 0)
print(f"""  (1) DUTY-SECTOR VALUE of every jet-weighted band pairing: 0 for EVERY family -- exactly, in Z,
      for the strict jet, the sectioned jet, and the jet-square repair.  STRUCTURAL, two mechanisms:
      B cap 7Z = B cap 13Z = {{0}} pins every duty-involving in-band incidence to the origin, and
      every divided-difference jet has J(0) = 0.
  (2) The gauge-free (strict) jet weight is the zero FUNCTION on the band: the strict rotation-Gram
      pairing is identically 0 before any localization question.
  (3) De-banded jet moments die mod 13 by the orbit power-sum kill (12 not| k+m); the only flat
      escapes are the Fermat projectors beta^0, beta^12 = jet-free sector counts.
  (4) Jet-faithful weights have NON-FLAT unit Gram spectra: the rotation lemma does not localize
      them (Step 5 witnesses).  Empirically Gj = 0 mod 13 on {gj_zero}/{len(COVBANK)} covering families,
      Gj2 = 0 on {gj2_zero}/{len(COVBANK)}: family-sensitive, NOT forced.
  (5) What the rotation lemma DOES force on the 91-line is the projector layer:
      N_B = 7*z7 and P12 = -C(z1,2) - z1*z7 + 6*C(z7,2) (mod 13), generically nonzero.""")
print(f"\nTOTAL: {'ALL CHECKS PASSED' if not FAILS else 'FAILURES: ' + str(FAILS)}")
