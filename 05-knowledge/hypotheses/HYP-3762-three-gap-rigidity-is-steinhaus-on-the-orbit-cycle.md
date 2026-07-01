---
id: HYP-3762
title: THE THREE-GAP RIGIDITY IS STEINHAUS ON THE ROTATION ORBIT-CYCLE (+ a support/wide-hole reduction, + the T_p band skeleton). The general rigidity g(n)<=3 (tight => <=3 distinct gaps mod n, the open core HYP-+2913) attacked three ways. (1) BAND-TRANSVERSAL: tight sets (AP, GW) AND the covering-min obey the SAME conditions T_p (radius-floor(p/n) covering of Z/p) for every prime p<=23 -- resonance-kill p<14, +-transversal 14<=p<28 -- the COMMON skeleton of both extremal families; but prime-T_p (p<=23) is NECESSARY not SUFFICIENT (4-gap sets obey it but are all NON-tight, M=4/31,2/15,1/7). (2) CYCLE route (EXACT for diff-closed): the difference-closed tight sets are the dilated APs c*{1..n-1} (HYP-3750), whose residues mod n ARE the rotation orbit {kc mod n}; the Steinhaus three-gap theorem applies VERBATIM => <=3 arc-lengths. The orbit is a single n-cycle (Cayley graph of Z/n by c); the three lengths are consecutive CF convergents + their sum (the Farey/Ostrowski cycle HYP-3746). (3) SUPPORT/WIDE-HOLE route (general): combinatorial identity #gaps = 1 + #distinct-run-lengths(missing) => 4 gaps needs missing>=6 (three distinct run-lengths 1+2+3) => support>=n-5 => <=3 gaps; and tight => high support via the wide-hole inequality HYP-3749 (missing>=4 => deep hole M>1/n, verified n=14). So the OPEN rigidity REDUCES to the cleaner support bound "tight => support >= n-5"
status: PARTIALLY-TRUE / REDUCTION. VERIFIED: (1) T_p common to tight+covering-min p<=23, and prime-T_p necessary-not-sufficient (n=14); (2) diff-closed tight = dilated-AP orbit, Steinhaus <=3 gaps (n=5,7,8,11,13,14; PROVED, it IS the three-gap theorem); (3) the run-length identity (PROVED elementary) + wide-hole support bound (missing>=4 => M>1/14 on 24000 samples). The general rigidity g(n)<=3 is REDUCED to "tight => support >= n-5" (a wide-hole/HYP-3749 target), and is EXACT (=Steinhaus) on the difference-closed family. Full general proof still open (rests on the wide-hole support bound at all n).
source: klein-2026-06-30-S51
depends_on:
  - HYP-+2913  # the three-gap (Steinhaus) rigidity g(n)<=3 -- the open core
  - HYP-3750   # tight-set classification: difference-closed = dilated AP
related:
  - HYP-3761   # n=14 census {AP,GW} (rigidity verified there as a byproduct)
  - HYP-3749   # the punctured-core wide-hole inequality (the support-bound mechanism)
  - HYP-3746   # the Farey-grid reach / Ostrowski cycle (the three lengths)
  - HYP-3736   # dense-core +-transversal lemma (the T_p at band primes)
  - HYP-3741   # witness hierarchy (the radius-r covering conditions)
  - HYP-3728   # Ihara prime-cycle / metazeta (the cycle framing)
results:
  - 04-computation/three_gap_rigidity_cycles_klein.py
  - 05-knowledge/results/three_gap_rigidity_cycles_klein.out
---

# HYP-3762 — the three-gap rigidity is Steinhaus on the rotation orbit-cycle

## The target (the open core)
The **three-gap rigidity** (`HYP-+2913`, the open core of the tight-set census): a tight set `S`
(`M(S)=1/n`, `|S|=n-1`) has a residue config mod `n` with **`<=3` distinct cyclic gaps**
(`g(n)<=3`). Verified computationally `n=4..7,14` (`n=14` as a byproduct of the complete census
`HYP-3761`); the general-`n` proof is open. This session attacks it three ways and **reduces** it.

## (1) The band-transversal `T_p` skeleton (the owner's hint) — VERIFIED, necessary-not-sufficient
Define `T_p` = the **radius-`floor(p/n)` covering of `Z/p`**: every rotation `a` has a runner within
`floor(p/n)` of `0` mod `p`. For tight sets (`M=1/n`) this is resonance-kill for `p<n` (radius 0) and
a `+-`transversal for `n<=p<2n` (radius 1).

**Both the tight sets (`AP`, `GW`) AND the covering-min `{1..12,182}` obey the SAME `T_p` for every
prime `p<=23`** (resonance-kill `p in {2,3,5,7,11,13}`, `+-`transversal `p in {17,19,23}`). So `T_p`
(`p<=23`) is the **common skeleton** of both extremal families of the two problems.

But **prime-`T_p` (`p<=23`) is NECESSARY, not SUFFICIENT** for `<=3` gaps: e.g.
`{1,5,6,10,11,12,15,19,20,21,23,25,26}` obeys all prime `T_p` (`p<=23`) yet has **4** gaps mod 14 —
and **every such counterexample is NON-tight** (`M=4/31, 2/15, 1/7`, killed by a deeper hole at a
composite/large `D`). So the rigidity genuinely needs the FULL tightness (the deeper-`D` holes), not
just the finite band conditions.

## (2) The CYCLE route — EXACT for the difference-closed family (it *is* Steinhaus)
The **difference-closed** tight sets are exactly the **dilated APs** `c*{1,...,n-1}`, `gcd(c,n)=1`
(`HYP-3750`). Their residues mod `n` are the **rotation orbit** `{c, 2c, ..., (n-1)c mod n}` of the
circle rotation by `c/n`. The **Steinhaus three-gap theorem** says any orbit `{0,α,2α,...,(N-1)α}`
partitions the circle into arcs of `<=3` lengths — so it applies **VERBATIM**: the difference-closed
tight sets have `<=3` gaps *by the three-gap theorem itself*. Verified `n=5,7,8,11,13,14` (all give
exactly 2 gaps).

**The cycle:** the orbit is a single `n`-cycle — the Cayley graph of `Z/n` with generator `c`. The
three gap-lengths are consecutive **continued-fraction convergent** denominators of `c/n` and their
sum (the Ostrowski/Farey cycle, `HYP-3746`; the Ihara prime-cycle framing `HYP-3728`). So the
rigidity, on the difference-closed family, is a **statement about the rotation orbit-cycle**, and it
is a *theorem* (Steinhaus), not a conjecture.

## (3) The SUPPORT / WIDE-HOLE route — reduces the general rigidity
**Combinatorial identity (elementary, proved):** for any residue config,
```
    #distinct gaps  =  1 + #distinct run-lengths of the MISSING residue set.
```
(A "run" is a maximal block of consecutive missing residues; a run of length `L` makes a gap `L+1`;
adjacent present residues make gaps `1`.) Hence **`>=4` gaps requires `>=3` distinct run-lengths**,
i.e. the missing set has `>= 1+2+3 = 6` residues (`support <= n-6`). Contrapositive:
```
    support >= n-5   =>   <=3 gaps.
```
And **tight => high support** by the wide-hole inequality (`HYP-3749`): a config missing `>=4`
residues has a wide hole, forcing a deep runner gap `M>1/n`. Verified `n=14`: 24000 sampled configs
with `support<=10` (missing `>=4`) ALL have `M>1/14`. (The census `HYP-3761` shows even more: tight
=> support in `{n-2,n-1}`.)

So the **open rigidity `g(n)<=3` REDUCES to the cleaner support bound `tight => support >= n-5`** — a
wide-hole statement (`HYP-3749`), far weaker than the empirical `support >= n-2`.

## Net
The three-gap rigidity is **Steinhaus on the rotation orbit-cycle**: EXACT for the difference-closed
tight sets (dilated APs = orbits, `<=3` gaps by the three-gap theorem), and for the general tight sets
it reduces to the support/wide-hole bound `support >= n-5` (via `#gaps = 1 + #distinct-run-lengths` +
`HYP-3749`). The band-transversal `T_p` conditions (`p<=23`) are the common skeleton of the tight sets
and the covering-min, but are necessary-not-sufficient. The open core is now pinned to one tractable
target: **the wide-hole support bound at all `n`**.
