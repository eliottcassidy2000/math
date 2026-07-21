---
id: THM-2033
title: "THE NC2 WALL IS THE CONFLUENT LIMIT OF THE TRANSITIVITY VANDERMONDE — the tournament↔NC2 bridge, completed and verified. NC2 noncancellation reduces (THM-1815) to det[(a_i+k)!] over the radial channel degrees, and det[(a_i+k)!] = ∏a_i!·Vandermonde(a) = the SIGNED TOURNAMENT SUM Σ_T sgn(T)x^{score(T)} (THM-1805; transitive tournaments the surviving terms) — VERIFIED exactly. DISTINCT degrees ⟹ Vandermonde≠0 ⟹ transitive channel (death-star HYP-8772) ⟹ noncancellation (codex THM-2017 degree-gap). The RESONANCE WALL = REPEATED degrees: the ordinary Vandermonde VANISHES, but the CONFLUENT Vandermonde (a derivative/Wronskian row replaces the coincident one) SURVIVES nonzero (verified) — the exact 'one order down' that codex's hyper-Bessel boundary and my Laguerre-Pólya reduction (HYP-8775, ODE θ²Φ=ξΦ) supply, and the mechanism behind MISTAKE-212 (tied channels are CONFLUENCE, not intransitivity; the degree preorder ties, and the derivative order breaks the tie). So death-star's channel-tournament lens, THM-1815's transitivity Vandermonde, codex's hyper-Bessel boundary, and my L-P boundary are ONE object — the confluence structure of the tournament sign-sum."
status: >
  VERIFIED (boxeph-2026-07-21-S203). (1) det[(a_i+k)!]=∏a_i!·Vandermonde(a) exact over 7 degree sets;
  (2) Vandermonde(x)=Σ_T sgn(T)x^{score(T)} exact (n=3:30, n=4:1440) — the transitivity Vandermonde
  IS the tournament sign-sum (THM-1805/1925); (3) the confluent limit: as degrees coincide the det → 0
  through its Vandermonde factor, and the derivative-row (confluent) det survives nonzero. The bridge
  is a SYNTHESIS/identification, not a proof of NC2: the uniform non-vanishing (distinct-degree case
  for all patterns) is still opus THM-1710's multinomial-ratio step (shared with TNC), and the wall
  (confluent case) is codex's hyper-Bessel resonance band + my L-P reduction — both OPEN in general.
  Respects MISTAKE-212 (transitivity is a sufficient domination certificate / regime classifier, not
  necessary for noncancellation; the tied core needs the confluent order).
source: boxeph-2026-07-21-S203 (owner: find connections between tournaments and NC2, long session, push/pull often)
depends_on: []
related:
  - THM-1815  # NC2 pair-in-radical = transitivity Vandermonde = moment-matrix det (mac-mini)
  - THM-1805  # Vandermonde = signed tournament sum, transitive survive (klein)
  - THM-1925  # my signed partition = ∏(x_j-x_i) = product of sines (the same Vandermonde)
  - THM-2017  # codex degree-gap = transitive channel (Vandermonde nonzero) = noncancellation
  - HYP-8772  # death-star channel-tournament lens (+ MISTAKE-212)
  - HYP-8775  # my Laguerre-Pólya hyper-Bessel boundary = the confluent structure
  - "07-reflections/the-vandermonde-is-the-bridge-tournaments-and-nc2-boxeph-S203.md"
script: 04-computation/confluent_transitivity_vandermonde_boxeph_S203.py (+ .out)
---

# THM-2033 — the NC2 wall is the confluent transitivity Vandermonde

## The bridge in one diagram

```
   NC2 noncancellation resultant  =  det[(a_i+k)!]  =  ∏a_i! · Vandermonde(a)  =  Σ_T sgn(T) x^{score(T)}
        (THM-1815)                     (verified)          (verified)              (THM-1805/1925)
                                                                                    ↑ transitive T survive
   distinct radial degrees  ⟺  Vandermonde ≠ 0  ⟺  TRANSITIVE channel  ⟺  noncancellation (THM-2017)
   repeated radial degrees  ⟺  Vandermonde = 0  ⟺  the WALL / tied core  ⟶  CONFLUENT Vandermonde
                                                                            (derivative/Wronskian, nonzero)
                                                                            = codex hyper-Bessel = my L-P (HYP-8775)
```

## 1. The transitivity Vandermonde is the tournament sign-sum (verified)

`det[(a_i+k)!]_{i,k=0..r-1} = (∏_i a_i!)·∏_{i<j}(a_j−a_i)` (verified on 7 degree sets), and
`∏_{i<j}(x_j−x_i) = Σ_{tournaments T} sgn(T)·x^{score(T)}` (verified n=3,4), with the **transitive**
tournaments the surviving terms (THM-1805; my THM-1925 leg (b): this Vandermonde on the unit circle is
a product of sines). So the NC2 moment-matrix determinant that must be nonzero (THM-1815) *is* the
signed tournament sum over the radial channel degrees.

## 2. The two regimes = transitive vs confluent

- **Distinct degrees** → `Vandermonde ≠ 0` → the channel-degree "tournament" is **transitive** (a
  clean total order, one source = the dominant channel) → the source's sign survives →
  `E[P^m] ≠ 0` = **noncancellation**. This is codex THM-2017's degree-gap regime and death-star
  HYP-8772's transitive channel.
- **Repeated degrees** (the resonance **wall**) → `Vandermonde = 0`. The degree preorder **ties**
  (MISTAKE-212: not a genuine tournament, and noncancellation can still hold — e.g. THM-2021's
  `E[P²]=24` on fully tied channels). The resolution is **confluence**: replace the coincident row by
  its derivative (a Wronskian). The confluent determinant is nonzero (verified) — the "detected one
  order down." That derivative order is exactly codex's `1/m` hyper-Bessel correction and my
  Laguerre–Pólya boundary ODE `θ²Φ = ξΦ` (HYP-8775).

## 3. What it unifies

Four threads are one object — the **confluence structure of the tournament sign-sum**:
`det[(a_i+k)!]` (THM-1815) / the signed tournament Vandermonde (THM-1805/1925) / death-star's
channel-tournament lens (HYP-8772, corrected by MISTAKE-212) / codex's hyper-Bessel boundary
(THM-2017) + my Laguerre–Pólya reduction (HYP-8775). Transitive = the cold, distinct-degree,
noncancelling regime; the wall = the confluent, tied, regular core. The NC2 residual is therefore:
prove the transitivity Vandermonde is uniformly nonzero at distinct degrees (opus THM-1710,
multinomial-ratio, shared with TNC), and control its confluent (derivative) limit at the tied core
(the hyper-Bessel / L-P boundary) — both the *same discriminant*, seen at distinct vs coincident nodes.
