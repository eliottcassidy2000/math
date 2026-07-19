#!/usr/bin/env python3
"""
lrc14_evaders_23_25_rungs_boxeph_S132.py  (HYP-7900 target; owner-directed)

RERUN THE S126 MOD-19 EVADERS AGAINST THE 23/25 RUNGS.

Bank: EXACT reproduction of lrc14_mod19_spread_kernel_boxeph_S126.py part (B)
(lcg seed 999, 6000 trials): 12-sets C = {3,35} u (10 from [3,35]), band-filled
at m=1 by construction, covering{2..12}; S126 found 990 such, of which 271 also
pass the mod-19 rung (19-mult present or full antipodal spread) = "the evaders".

Rung logic (general-p spread lemma, S115/S126/S131): a 3/38 candidate must have,
for every p in {13,17,19,23,25} (all with 2/p > 3/38):
    p-BLOCKED (multiple of p present; for 25 this means a multiple of 25)
 or p-SPREAD  (unit residues hit all antipodal unit pairs of Z/p).
Otherwise M >= 2/p > 3/38: the family is KILLED as a 3/38 candidate by rung p.

Outputs:
 (1) the 271 evaders classified at 23 and 25 (and 13, 17): kill fractions;
 (2) full-stack survivors among the 990: blocker census, exact M (integer pinch,
     q <= 70 complete), maximizer modulus + cheapest needle q;
 (3) fresh larger sample (seed 424242, 60000 trials) for tighter statistics,
     incl. empirical P(23-spread | no 23-mult) on band+covering families vs the
     1/1191 uniform-slot figure (HYP-7895);
 (4) explicit fully-blocked family: the bounded-rung stack provably cannot
     close the band (needles remain necessary).

boxeph-2026-07-19-S132. Pure Python, exact integers.
"""

from math import gcd

# ---------- exact copies of the S126 bank generators (determinism contract) ----

def lcg(seed):
    x = seed
    while True:
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        yield x

def sample(gen, lo, hi, k, exclude=()):
    out = []
    while len(out) < k:
        v = lo + next(gen) % (hi - lo + 1)
        if v not in out and v not in exclude:
            out.append(v)
    return out

def covering(sp):
    return all(any(v % q == 0 for v in sp) for q in range(2, 13))

# ---------- rung machinery ----------------------------------------------------

def upairs(p):
    return [j for j in range(1, (p + 1) // 2) if gcd(j, p) == 1 and gcd(p - j, p) == 1]

def rung_status(sp, p):
    """'blk' (p | some speed; p=25 needs 25|v), 'sprd', or 'KILL' (=> M >= 2/p)."""
    if any(v % p == 0 for v in sp):
        return 'blk'
    hit = {min(v % p, p - v % p) for v in sp if gcd(v % p, p) == 1}
    return 'sprd' if all(j in hit for j in upairs(p)) else 'KILL'

RUNGS = [13, 17, 19, 23, 25]

def M_int(sp):
    """exact M for speeds <= 35: all q <= 70 (pair sums) complete by Pinch.
    returns (num, den, qstar) with M = num/den at reduced denominator, and the
    attaining q; also cheapest needle = min q with min-dist/q > 3/38."""
    bn, bd, qstar = 0, 1, None
    cheapest = None
    for q in range(2, 71):
        best_here = 0
        for b in range(1, q // 2 + 1):
            md = q
            for c in sp:
                r = (c * b) % q
                d = r if r <= q - r else q - r
                if d < md:
                    md = d
                    if md * bd < bn * q and md * 38 <= 3 * q:
                        break
            if md > best_here:
                best_here = md
            if md * bd > bn * q:
                bn, bd, qstar = md, q, q
        if cheapest is None and best_here * 38 > 3 * q:
            cheapest = q
    g = gcd(bn, bd)
    return bn // g, bd // g, qstar, cheapest

# ---------- (1)+(2): the S126 bank, reproduced --------------------------------

def build_bank(seed, trials):
    g2 = lcg(seed)
    bank = []
    for _ in range(trials):
        rest = sample(g2, 3, 35, 10, exclude=(3, 35))
        C = tuple(sorted({3, 35} | set(rest)))
        if len(C) != 12 or not covering(C):
            continue
        bank.append(C)
    return bank

def classify(bank, label):
    n = len(bank)
    ev19 = [C for C in bank if rung_status(C, 19) != 'KILL']
    print("%s: band+covering families: %d; pass mod-19 rung (evaders): %d" % (label, n, len(ev19)))
    # classification of the evaders at each further rung
    for p in (13, 17, 23, 25):
        st = {'blk': 0, 'sprd': 0, 'KILL': 0}
        for C in ev19:
            st[rung_status(C, p)] += 1
        print("  evaders at p=%2d: blocked=%4d  spread=%4d  KILLED=%4d  (kill %.1f%%)"
              % (p, st['blk'], st['sprd'], st['KILL'],
                 100 * st['KILL'] / max(1, len(ev19))))
    k2325 = [C for C in ev19 if rung_status(C, 23) == 'KILL' or rung_status(C, 25) == 'KILL']
    kany = [C for C in ev19 if any(rung_status(C, p) == 'KILL' for p in (13, 17, 23, 25))]
    surv = [C for C in ev19 if all(rung_status(C, p) != 'KILL' for p in RUNGS)]
    print("  evaders KILLED by 23-or-25 rung: %d/%d (%.1f%%)"
          % (len(k2325), len(ev19), 100 * len(k2325) / max(1, len(ev19))))
    print("  evaders killed by ANY of 13/17/23/25: %d/%d (%.1f%%)"
          % (len(kany), len(ev19), 100 * len(kany) / max(1, len(ev19))))
    print("  FULL-STACK SURVIVORS (pass all five rungs): %d/%d of bank (%.2f%%)"
          % (len(surv), n, 100 * len(surv) / max(1, n)))
    # 23-spread coincidence rate among no-23-mult band+covering families
    no23 = [C for C in bank if all(v % 23 for v in C)]
    sp23 = [C for C in no23 if rung_status(C, 23) == 'sprd']
    print("  P(23-spread | no 23-mult, band+covering) = %d/%d = %.4f  [uniform-slot: 1/1191 = 0.00084]"
          % (len(sp23), len(no23), len(sp23) / max(1, len(no23))))
    return ev19, surv

print("=" * 88)
print("S126 BANK REPRODUCTION (seed 999, 6000 trials) + 23/25 RUNGS   [boxeph-S132]")
print("=" * 88)
bank0 = build_bank(999, 6000)
ev19, surv = classify(bank0, "S126 bank")

print("\n--- survivor census (S126 bank): blockers carried + exact M + cheapest needle ---")
blocker_hist = {}
needle_hist = {}
worst = []
for C in surv:
    blk = tuple(p for p in RUNGS if any(v % p == 0 for v in C))
    blocker_hist[blk] = blocker_hist.get(blk, 0) + 1
    nu, de, qs, ch = M_int(list(C))
    needle_hist[ch] = needle_hist.get(ch, 0) + 1
    worst.append((nu / de, nu, de, qs, ch, C))
    assert (nu, de) != (3, 38), ("3/38 REALIZED?!", C)
for blk, c in sorted(blocker_hist.items(), key=lambda x: -x[1]):
    print("  blockers %-20s: %d" % (str(blk), c))
print("  cheapest-needle spectrum (q: count): %s"
      % dict(sorted(needle_hist.items())))
worst.sort()
print("  lowest-M survivors (closest to 3/38=0.0789):")
for Mv, nu, de, qs, ch, C in worst[:5]:
    print("    M=%d/%d=%.4f  maximizer q=%d  cheapest needle q=%s  C=%s" % (nu, de, Mv, qs, ch, C))

print("\n" + "=" * 88)
print("FRESH SAMPLE (seed 424242, 60000 trials) -- statistics at scale")
print("=" * 88)
bank1 = build_bank(424242, 60000)
ev19b, survb = classify(bank1, "fresh bank")
blocker_hist = {}
for C in survb:
    blk = tuple(p for p in RUNGS if any(v % p == 0 for v in C))
    blocker_hist[blk] = blocker_hist.get(blk, 0) + 1
print("  fresh-bank survivor blocker census (top 6):")
for blk, c in sorted(blocker_hist.items(), key=lambda x: -x[1])[:6]:
    print("    blockers %-20s: %d" % (str(blk), c))

print("\n" + "=" * 88)
print("STRUCTURAL: the bounded-rung stack CANNOT close the band (explicit family)")
print("=" * 88)
X = [3, 9, 10, 11, 13, 17, 19, 23, 24, 25, 33, 35]
ok_band = all(3 <= v <= 35 for v in X) and 3 in X and 35 in X
print("X = %s" % X)
print("  band[3,35] + straddle: %s ; covering{2..12}: %s" % (ok_band, covering(X)))
for p in RUNGS:
    print("  rung p=%2d: %s" % (p, rung_status(X, p)))
nu, de, qs, ch = M_int(X)
print("  exact M(X) = %d/%d = %.4f  (maximizer q=%d; cheapest needle q=%d)"
      % (nu, de, nu / de, qs, ch))
print("  => X passes band + covering + ALL five spread rungs, yet M(X) != 3/38:")
print("     killed only CROSS-MODULUS. The in-band certificate residual is NONEMPTY;")
print("     medium-needle arguments remain NECESSARY after the full bounded-rung stack.")
print("\nDONE.")
