---
id: HYP-3702
title: CREATIVE TAXONOMY of covering-set special cases for the covering-min problem -- partition covering sets into EASY families (provable M >> 1/n: non-covering/q-witness incl. the consecutive cusp; SHIFTED CONSECUTIVE BLOCKS {m,..,m+n-2} with m>=2 which are covering and have M>=1/8 at n=14 with margin growing in m -- a clean 'renormalization flow' from the m=1 cusp to easy; spread/large-min-speed M~0.3-0.46), the CUSP (even-heavy/full-residue, measure 0, existence/units), and the HARD CORE (the binding family: '{1}/small dense core + block + outlier' near-construction, M=2/(2n-1) for n<=6 / n/Phi_6(n) for n>=7, transition at n=7=apex, opus-confirmed). The hardness axis is LOWNESS (small speeds, esp. containing 1): dense-LOW is binding, dense-HIGH(shifted)/spread is easy. Novel families to test: {1}+AP (three-distance), Beatty/Sturmian, Sidon, lacunary, divisor-closed. The proof decomposes: discharge the easy families, isolate the '{1}+block+outlier' hard core, attack it via the hexagonal/Eisenstein (Kershner) + three-distance structure
status: COMPUTED + verified (covering-set M-classification n=5,6 exhaustive; shifted-block M(m) curve n=14). The easy families (non-covering, shifted-blocks m>=2, spread) have large margin; the hard core matches the n-dependent covering-min (HYP-3701, opus-confirmed). Taxonomy + novel-family proposals; the easy-family bounds are computationally clear, proofs sketched not done.
source: mac-mini-2026-06-30-S43
related:
  - HYP-3701  # the covering-min extremal family transitions (the hard core)
  - HYP-3706  # klein-S24: hexagonal/Eisenstein/Kershner (the hard-core attack route)
  - THM-523   # the covering reduction + q-witness (the non-covering easy family)
  - HYP-3700  # the gap-edge (isolated); this is the M-edge taxonomy
results:
  - 04-computation/covering_min_extremal_family_macmini_20260630.py  # (related)
---

# HYP-3702 -- a taxonomy of covering-set special cases

The owner asked for useful subconditions / special cases (dense covering sets, etc.). The covering
reduction (THM-523, LRC <=> M(S)>=1/n for covering S) is one big quantifier; the right move is to PARTITION
covering sets into families, discharge the easy ones, and isolate the hard core. Grounded by exhaustive
M-classification (n=5,6) and the shifted-block curve (n=14).

## EASY families (provable / large margin) -- discharge these
1. **NON-COVERING (q-witness, THM-523):** misses a multiple of some `q<=n` => `t=1/q` is lonely,
   `M>=1/q>=1/n`. This includes the consecutive CUSP `{1,..,n-1}` (`M=1/n` exactly but NOT covering -- no
   multiple of `n`). The tightest M-config is non-covering, hence off-path.
2. **SHIFTED CONSECUTIVE BLOCKS `{m,..,m+n-2}`, `m>=2`:** a clean "renormalization flow." `m=1` is the cusp
   (`M=1/n`, non-covering); every COVERING shift (`m>=2`) has `M` jumping up and growing with `m`. Verified
   n=14: `M = 1/8, 1/6, 5/22, 2/7, ..., ~5/13, ~7/17` for `m=2,3,5,8,...,28` -- MIN over covering blocks is
   `1/8 >> 1/14`, margin grows monotonically. So consecutive blocks are entirely SAFE (the only tight one is
   the non-covering cusp). Provable target: `{m,..,m+n-2}` covering (`m>=2`) `=> M>=1/8`.
3. **SPREAD / large-min-speed:** speeds spread out (large `min`, large gaps) `=> M ~ 0.3-0.46` (avg over
   covering sets, n=5,6). The vast majority of covering sets. Provable via a packing/measure argument (few
   small speeds => coarse danger => large safe region).

## The HARD CORE (binding, `M -> 1/n`) -- isolate and attack
4. **"{1}/small dense core + block + outlier" near-construction:** the binding covering sets contain a small
   speed (`1`, sometimes `2`) plus a near-consecutive block, plus an outlier covering the top `q`'s.
   Examples (covering-min): `{1,3,4}` (n=4), `{1,3,4,5}`/`{1,4,5,6}` (n=5), `{1,3,4,5,18}`/`{1,4,5,6,7}`
   (n=6) `= 2/(2n-1)`; at `n>=7` the `n(n-1)`-construction `{1,..,n-2,n(n-1)} = n/Phi_6(n)` wins (n=14:
   `14/183`, opus-confirmed unique hardest). The EXTREMAL FAMILY TRANSITIONS at `n=7=apex` (HYP-3701). This
   is the only family near `1/n`; the proof's whole difficulty lives here.

## The CUSP (measure 0) -- existence, not measure
5. **EVEN-HEAVY / full-residue:** deep 2-adic descent / the full-`Z_p` core; lonely MEASURE `-> 0`, carried
   by EXISTENCE (the `phi(n)` units, HYP-3615; the gap-edge isolated, HYP-3700). Distinct from the M-edge.

## The hardness axis and the dichotomy
The data shows **LOWNESS is the hardness axis**: dense-LOW (small speeds, especially containing `1`) is
binding (`M ~ 1/n`); dense-HIGH (shifted blocks, large speeds) and spread are easy (`M ~ 1/2`,`~0.35`). The
smallest speed and the presence of a small consecutive core near `t=1/n` create the tightness; pushing the
block UP (the shift `m`) or SPREADING it relaxes `M`. So the covering constraint is only dangerous when it
is satisfiable with SMALL, DENSE, LOW speeds -- which is exactly the near-consecutive core.

## Novel special cases worth testing (creative)
- **`{1} + AP` / `{1} + block`** (the binding shape): the gaps of an AP under `t` are governed by the
  THREE-DISTANCE (Steinhaus) theorem -- a possible CLOSED FORM for `M` of the hard core, hence a direct
  `>=1/n` proof.
- **Beatty / Sturmian covering sets** (`{floor(k*alpha)}`): three-distance gap structure; the consecutive
  is `alpha=1`.
- **Sidon / `B_h` covering sets** (distinct pairwise differences): danger zones don't resonantly overlap
  `=> M large` (likely easy).
- **Lacunary / geometric** (`{1,r,r^2,..}`): equidistribution `=> M` structure.
- **Divisor-closed / lcm-structured** (`S subset divisors of lcm(1..n)`): each `q` covered by `q` itself;
  highly resonant -- test whether these are hard.

## The proof decomposition (what this buys)
LRC(n) covering reduction splits into: (i) NON-COVERING (q-witness, done); (ii) SHIFTED BLOCKS `m>=2`
(provable `M>=1/8`); (iii) SPREAD (packing argument); (iv) the CUSP (existence/units); (v) the HARD CORE
(`{1}+block+outlier`, the `n`-dependent construction). Only (v) is open, and it is the hexagonal/Eisenstein
covering-optimality problem (klein-S24/opus: the construction's `speeds*t* = ` an AP mod `Phi_6` with
min-distance `n`, `t*=zeta_6`; Kershner thinnest covering). So the special-case taxonomy reduces the whole
conjecture to ONE family -- the near-consecutive `{1}+block+outlier` core -- and one tool (hexagonal
covering optimality + three-distance for the AP gaps).
