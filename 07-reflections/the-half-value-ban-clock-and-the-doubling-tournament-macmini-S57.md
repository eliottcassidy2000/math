# The half-value ban clock: doubling orbits, voltage lifts, and a tournament worth defining

**mac-mini-2026-07-05-S57 (HYP-4132).** Developing the owner's parity-split seed
(E_p = {0,±2} always; O_p = {±1} or {±1,±3}) into exact structure.  Data:
`05-knowledge/results/lrc_halfvalue_ban_macmini_S57.out`.  Complementary to
kps-4137 (their reduction/dividend program; my counting dual + the object).

## 1. The exact criterion (verified 7,001 set-prime pairs, zero mismatches)

At witness modulus s = 2p (p odd prime, margin 2/25, μ = ⌈4p/25⌉):
dilations are odd units mod 2p ↔ units mod p.  Define the HALF-VALUE
ρ(v) = v (v odd), v/2 (v even), taken mod p up to sign.  Then runner v bans
exactly the dilation ±class ρ(v)^{-1} (μ=3: p ∈ {13,17}), plus ±3ρ(v)^{-1}
for odd v when μ=4 (p ∈ {19,23}).  Odd p-multiples ban nothing (antipode
distance p); runners ≡ 0 mod 2p kill every dilation.

**CRITERION: a margin-2/25 witness at s = 2p exists  ⟺  no runner ≡ 0 (mod 2p)
and the ban classes fail to cover the (p−1)/2 dilation ±classes.**

Everything is a set-cardinality statement about {±ρ(v)^{-1} mod p} — decidable,
12 modular inversions and a set-size compare.  This is the atom kps's
TemplateWitness wants at the 2p moduli: no dilation search needed at all.

## 2. Why the owner's "0 ∈ E_p always" is the right lever

An even runner blocks through its HALF on the p-grid — including at 0.  So
2p | v kills the modulus outright, and the even/odd asymmetry (3 blocked
residues vs 2) is the "parity dividend" (kps).  The counting dual: each
runner is one BAN; bans coincide exactly when half-values agree up to sign —
same-parity: v ≡ ±w (mod p); mixed: w ≡ ±2v (mod p), the DOUBLING FOLD.
Doubled-residue families (the census's dominant stratum) are exactly the
ban-colliders — which is WHY the census's first-witness mass sits at q = 26.

## 3. The AP's rigidity, seen from the clock

At p = 13 the AP {1..12} has ban fibers {1,2},{3,6},{5,10},{4,11},{8,9},{7,12}:
each fiber = one odd + one even runner, paired by v ↦ 2v.  The load is
PERFECTLY UNIFORM (2 on each of the 6 classes) — the AP is the perfect double
round-robin of the ban schedule, and ×2 acts on the ±classes as a single
6-cycle (ord±(2) = 6).  Tightness of the AP at this modulus = perfect balance
of the ban tournament.  The four gears p = 13, 17, 19, 23 have doubling-cycle
types 6, 4+4, 9, 11 — a gear system the template family turns.

## 4. THE NEW OBJECT: the ban-load tournament B_p(W)

Vertices: the (p−1)/2 dilation ±classes on the doubling cycle σ: c ↦ 2c.
Load: ℓ_W(c) = #{v ∈ W : ±ρ(v)^{-1} = c} (the pushforward of W under the ban
map).  Orient each σ-edge c → σc iff ℓ(c) ≥ ℓ(σc), ties along the cycle
(this is a Tournament Analysis switchboard in the canon's exact sense —
THM-372: pairwise observable = load comparison, tie path = the doubling
cycle).  Facts:
  - clearing at 2p ⟺ ℓ has a zero ⟺ the empty class is a load-minimum sink;
  - the AP ⟹ all ties ⟹ B_p = the pure clock orientation;
  - dilating W permutes loads by the unit action: B_p(uW) = rotated B_p(W) —
    the tournament is a TEMPLATE invariant (depends on W mod p only);
  - the parity split is a VOLTAGE LIFT in the canon's THM-378 sense: the
    mod-2p dilation space double-covers mod-p; runners lift by parity;
    the ban structure descends through the cover.  The even-graph duality
    (E_n first-class) echo is exact: even runners live on the halved sheet.

Invariant to pursue: the load spectrum (ℓ sorted) across the four gears as a
TEMPLATE SIGNATURE; H(B_p) as a hardness measure.  Rigidity hypothesis
(logged as HYP-4132-H1): profile-passing families with UNIFORM load at p=13
are AP-shaped mod 26 — uniformity is the tight-side fingerprint.

## 5. What the tests settled (and rerouted)

- T3: 36.3% of 12-subsets of [1,30] are ρ-full at 13 — so q=26 clears ~64%
  alone; the multi-gear family is genuinely needed.
- T4: profile-passing families blocked at ALL FOUR gears {26,34,38,46} EXIST
  (e.g. {1,15,18,25,33,37,38,40,42,44,46,48}) — kps's (ii) 5-modulus hope is
  answered NEGATIVELY.  But each such family witnesses at DOZENS of other
  moduli in [28,60] (odd primes 29..59, 4p-types, ...): Q50 unthreatened.
  The right template family is the full gear-box, with the 2p-gears carrying
  the criterion above and the other moduli needing their own split rules
  (at s = 4p: runners split by v mod 4; at odd prime s: no split — the
  gcd-type stratification is the general parity dividend).

## 6. Why this matters for the finish

The template dichotomy's verification cost is dominated by dilation searches
per modulus.  The criterion removes the search at the 2p-gears entirely
(set-cardinality), and the ban-load view suggests each gear-type has its own
search-free criterion (the gcd-type split).  If every gear in [26,50] gets
one, TemplateDichotomy becomes a conjunction of counting statements over
residue data — the cheapest possible decidable form, and the load-uniformity
rigidity (H1) is the tight-branch's fingerprint inside the same language.
