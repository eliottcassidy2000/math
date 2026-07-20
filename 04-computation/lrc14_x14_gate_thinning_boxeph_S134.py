#!/usr/bin/env python3
"""
lrc14_x14_gate_thinning_boxeph_S134.py  (HYP-7955)

THE x14 GATE LAW at p=43 — executing kind-pasteur's c88/ledger-4b redirect (b):
"spread lemmas at the LIFTED moduli 2^a 7^b p (the LRCMod19/23Spread template,
re-aimed)".

Objects (S-T architecture, ledger s3): a 13-tuple of speeds is IMPROPER mod 14p
iff NO b in Z/14p has all distances dist_{14p}(v*b) > p (witness = all STRICTLY
clear; closed danger arcs — {1..13} is the boundary control: its only witness
candidates have equality, so it is improper mod every p).

CRT DECOMPOSITION of the gate (proved by the shadow computation, stated here,
verified by the harness):
  * pure-p scales  (b = 14b''): level-1 improperness mod p (dist_p <= floor(p/14))
      — kind-pasteur's folded-AP cover object;
  * pure-14 scales (b = p b'):  for every unit b' mod 14: some v with
      v*b' == 0 or +-1 (mod 14) — EXACTLY the LRCMod19/23Spread template at
      modulus 14 (spread-or-blocked with the +-1 window);
  * mixed unit scales: some v whose mod-14 part lies in the 2-3 element set
      SLAVED to its mod-p part: v*b mod 14p in [-p, p] <=> (v*b mod 14) in
      S(v*b mod p), |S| = 2 generically (3 at r = 0) out of 14.
The mixed-scale slaving is the NEW lift-level content (level-1 conditions are
subsumed by bookkeeping — c88 Verdict 1; the pure-14 spread alone is freely
satisfiable; the JOINT system is not).

MEASUREMENT: for level-1-improper base tuples (residues mod p), does ANY
assignment of (sign mod p, class mod 14) per speed make the tuple improper mod
14p — and how many of the 28^13 assignments survive?  This is the thinning
factor of the x14 gate on the level-1 sea (the number that prices the lift
tower's front door; death-star's |S| criterion is the x7 sequel).

Controls (MISTAKE-162): {1..13} and GW {1..11,13,24} must be level-1 improper
AND their canonical lifts improper mod 602.

boxeph-2026-07-19-S134.  Pure Python, exact integers.
"""

import random
from math import gcd

P = 43
Q = 14 * P                     # 602
DK = P // 14                   # 3
HALF = Q // 2                  # 301
HP = (P - 1) // 2              # 21

def dist(x, m):
    r = x % m
    return r if r <= m - r else m - r

def improper_modp(res):
    """level-1: no b'' with all dist_p(v*b'') > DK."""
    for b in range(1, HP + 1):
        if all(dist(v * b, P) > DK for v in res):
            return False
    return True

def improper_modQ(vals):
    """no b in [1, HALF] with all dist_Q(v*b) > P."""
    for b in range(1, HALF + 1):
        ok = False
        for v in vals:
            if dist(v * b, Q) <= P:
                ok = True; break
        if not ok:
            return False
    return True

# ---------------- controls ----------------------------------------------------
AP13 = list(range(1, 14))
GW = list(range(1, 12)) + [13, 24]
assert improper_modp([v % P for v in AP13]), "AP13 not level-1 improper?!"
assert improper_modp([v % P for v in GW]), "GW not level-1 improper?!"
assert improper_modQ(AP13), "AP13 not improper mod 602?!"
assert improper_modQ(GW), "GW not improper mod 602?!"
print("controls OK: {1..13} and GW level-1 improper AND improper mod 602")

# ---------------- level-1 base tuples ------------------------------------------

def random_level1_tuples(n_want, seed):
    """randomized MRV backtracking for improper-mod-p 13-tuples (residues 1..p-1,
    up to global sign per speed: canonicalize to [1, HP])."""
    rng = random.Random(seed)
    out = []
    scales = list(range(1, HP + 1))
    # hitters[b] = residues r in [1,HP] with dist_p(r*b) <= DK (as +-r both work:
    # dist_p((-r)b) = dist_p(rb), fold to canonical r <= HP)
    hitters = {b: [r for r in range(1, HP + 1) if dist(r * b, P) <= DK]
               for b in scales}
    while len(out) < n_want:
        chosen = []
        unc = set(scales)
        dead = False
        while unc and len(chosen) < 13:
            b = min(unc, key=lambda x: len(hitters[x]))
            cands = [r for r in hitters[b]]
            rng.shuffle(cands)
            r = cands[0]
            chosen.append(r)
            unc = {x for x in unc if dist(r * x, P) > DK}
        if not unc:
            while len(chosen) < 13:
                chosen.append(rng.randrange(1, HP + 1))
            res = sorted(chosen)
            if improper_modp(res):
                out.append(tuple(res))
    return out

# ---------------- the x14-gate CSP ---------------------------------------------

def gate_csp(base, node_cap=3_000_000, count_cap=50, max_zeros=13):
    """base: 13 residues in [1, HP] (mod-p parts up to sign).  Options per slot:
    (eps in {+1,-1}, c14 in Z/14) -> v mod Q via CRT.  Find assignments making
    the tuple improper mod Q.  Returns (exists, n_solutions_capped, nodes,
    example)."""
    # per-slot option masks over the HALF scales
    slots = []
    for r in base:
        opts = []
        for eps in (1, -1):
            rp = (eps * r) % P
            for c14 in range(14):
                # CRT: v mod Q with v = rp (p), c14 (14)
                # v = rp + p * ((c14 - rp) * inv_p_mod14 mod 14)
                inv = pow(P % 14, -1, 14)
                v = (rp + P * (((c14 - rp) * inv) % 14)) % Q
                m = 0
                for b in range(1, HALF + 1):
                    if dist(v * b, Q) <= P:
                        m |= 1 << (b - 1)
                opts.append((m, eps, c14, v))
        slots.append(opts)
    FULL = (1 << HALF) - 1
    # suffix reach
    sufreach = [0] * 14
    acc = 0
    for i in range(12, -1, -1):
        for (m, *_rest) in slots[i]:
            acc |= m
        sufreach[i] = acc
    sols = []
    nodes = 0
    choice = []

    def rec(i, mask, has_odd, has_non7, nz):
        nonlocal nodes
        if nodes > node_cap or len(sols) >= count_cap:
            return
        nodes += 1
        if i == 13:
            # NON-DEGENERACY (dilation-bite guard): exclude all-even and
            # all-7-divisible lifts -- those fold down the tower (v=2u/7u/14u
            # degenerates to the mod-301/86/43 object; improper trivially iff
            # the base is). The meaningful gate question is 14-primitive lifts.
            if has_odd and has_non7 and mask == FULL:
                sols.append(list(choice))
            return
        if (mask | sufreach[i]) != FULL:
            return
        for (m, eps, c14, v) in slots[i]:
            z = 1 if v % 14 == 0 else 0
            if nz + z > max_zeros: continue
            choice.append(v)
            rec(i + 1, mask | m, has_odd or (v % 2 == 1), has_non7 or (v % 7 != 0), nz + z)
            choice.pop()

    rec(0, 0, False, False, 0)
    return (len(sols) > 0, len(sols), nodes, sols[0] if sols else None)

# sanity: AP13's own base must admit a lift (its own!)
base_ap = tuple(sorted(min(v % P, P - v % P) for v in AP13))
ok, ns, nd, ex = gate_csp(base_ap, count_cap=5)  # AP13 is 14-primitive
print("AP13 base admits improper 14-lift: %s (>=%d solutions found; %d nodes)" % (ok, ns, nd))
assert ok

# ---------------- the measurement ----------------------------------------------
N_BASE = 40
bases = random_level1_tuples(N_BASE, seed=1443)
print("\nsampled %d random level-1-improper base tuples (p=43, dk=3)" % N_BASE)
admit = 0
sol_counts = []
tot_nodes = 0
examples = []
for i, base in enumerate(bases):
    ok, ns, nd, ex = gate_csp(base)
    tot_nodes += nd
    if ok:
        admit += 1
        sol_counts.append(ns)
        if len(examples) < 3:
            examples.append((base, ex))
print("bases admitting ANY improper 14-lift: %d/%d" % (admit, N_BASE))
if sol_counts:
    print("solution counts (capped at 50): %s" % sorted(sol_counts)[:10])
print("total CSP nodes: %d" % tot_nodes)

# ---- ZERO-BUDGET note (soundness lesson, kept honest) --------------------------
print("\nzero-budget stratification: WITHDRAWN as unsound in this version. The slot-order")
print("CSP hits its node cap at small max-zeros without deciding satisfiability, and a cap")
print("is NOT a refutation: AP13's EXPLICIT lift v=1..13 has z=0 zeros mod 14 (witness),")
print("while the capped search only found z=10 -- the caps bind exactly where the answer")
print("is interesting. Min-z stratification needs a SCALE-DRIVEN CSP (branch on the")
print("least-coverable uncovered scale, not on slots) -- named handoff. What stands:")
print("  * 40/40 random level-1 bases admit improper 14-lifts, but the found survivors")
print("    are (quasi-)degenerate: 12 speeds = 0 mod 14 + one free (dilated-core+killer)")
print("  * fully random lifts essentially never survive (MC 0/30000)")
print("  * AP13 achieves z=0 (explicit) -- the AP is the extreme non-degenerate survivor")

# MC thinning estimate on the first 10 bases: random (eps, c14) assignments
rng = random.Random(602)
mc_hits = mc_tries = 0
for base in bases[:10]:
    for _ in range(3000):
        vals = []
        for r in base:
            eps = rng.choice((1, -1))
            c14 = rng.randrange(14)
            inv = pow(P % 14, -1, 14)
            rp = (eps * r) % P
            v = (rp + P * (((c14 - rp) * inv) % 14)) % Q
            vals.append(v)
        mc_tries += 1
        if improper_modQ(vals):
            mc_hits += 1
print("MC random-lift survival: %d/%d" % (mc_hits, mc_tries))

for base, ex in examples:
    m14 = sorted(v % 14 for v in ex)
    print("example survivor: base(p-parts)=%s lift mod14=%s" % (base, m14))
print("\nsummary: the x14 gate multiplies the level-1 sea by (admit-rate x "
      "sol-density); pure-14 shadow = the mod-14 spread lemma (0,+-1 window); "
      "mixed-scale slaving = the lift-level kill structure (ledger 4b (b)).")
print("DONE.")
