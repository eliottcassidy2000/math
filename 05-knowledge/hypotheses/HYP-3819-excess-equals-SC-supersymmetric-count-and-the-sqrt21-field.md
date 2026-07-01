---
id: HYP-3819
title: THE FLIP-RANK EXCESS = #{SC classes with |Aut|>n} (verified n=4..7 = 0,0,1,3, exact) -- each self-complementary super-symmetric class forces a covering dimension because SC = complement-FIXED (can't pair with a complement to share a dimension), while non-SC super-symmetric classes come in complement pairs {T,T^op} (the antipodal map on the arc-cube, THM-584) and share; the SC condition is ESSENTIAL (at n=7 only 3 of the 5 |Aut|>7 classes are SC, and excess=3). CORRECTION: my S81 formula kappa(7)=1+C(5,2)=11 is WRONG -- opus/klein found kappa(7)=12 (the lazy-caterer 1+C(n-2,2) BREAKS at n=7); excess(7)=12-9=3. + the iota-odd certificate lives in the BIQUADRATIC field Q(sqrt-3, sqrt-7), Galois group Z2xZ2 = the involution-atlas Klein four (S88), whose three quadratic subfields are Q(sqrt-3) [Eisenstein/hexagonal/Phi6/Dedekind margin], Q(sqrt-7) [apex-7/Klein-quartic-PSL(2,7)/X0(14)/cusp f14], and the REAL Q(sqrt21) with 21 = 3*7 = the forbidden H value {7,21} = C(7,2) -- so sqrt21 = the real coupling of the hexagonal(3) and apex(7) involutions; h(-3)=h(-7)=1 (both PID) matching X0(14) genus 1.
status: MIXED -- excess = #{SC & |Aut|>n} CONFIRMED n=4..7 (exact, 0,0,1,3); the SC-essential fact CONFIRMED (n=7: non-SC super-symmetric don't force excess); the complement-pairing rigorization is a PROOF SKETCH (mechanism identified, not fully rigorous). excess(8)=#{SC & |Aut|>8} = OPEN (the order-7 Paley family is non-SC => 0 contribution; the SC super-symmetric come from order-3/5 structures, full n=8 enumeration deferred). The Q(sqrt-3,sqrt-7)/sqrt21 Klein-four synthesis is exact number theory + a structural connection to the involution atlas, not a certificate proof.
source: mac-mini-2026-07-01-S90c
related:
  - HYP-3817   # S90 fixed-point instruments (this quantifies the |Aut| needles = the covering excess)
  - HYP-3810   # S85 T-join obstruction to SC covers (the SC classes are the fold-fixed rigid core)
  - HYP-3814   # S88 involution atlas (the Klein four; Q(sqrt-3,sqrt-7) is another realization)
  - HYP-3805   # opus/klein flip-rank k(n) + the S_n-folding excess (kappa(7)=12, not 11)
  - HYP-3798   # S81 kappa=1+C(n-2,2) -- CORRECTED here: breaks at n=7 (11 vs actual 12)
results:
  - 04-computation/excess_sc_supersymmetric_macmini_20260701.py
  - 05-knowledge/results/excess_sc_supersymmetric_macmini_20260701.out
---

# HYP-3819 -- the excess is the self-complementary super-symmetric count; and the sqrt21 field

## The main result: excess(n) = #{SC classes with |Aut| > n}
Per iso class (tiling model, exhaustive n<=7): compute `|Aut|` and self-complementary (SC) status. The
flip-rank excess (opus/klein `k(n) - `max classical bound `= 0,0,0,1,3` for `n=3..7`) equals exactly the
number of **self-complementary classes with `|Aut| > n`**:

| n | #classes | `|Aut|` distribution | `#{|Aut|>n}` | `#{SC & |Aut|>n}` | excess |
|---|---|---|---|---|---|
| 4 | 4 | `{1:2,3:2}` | 0 | 0 | 0 |
| 5 | 12 | `{1:7,3:4,5:1}` | 0 | 0 | 0 |
| 6 | 56 | `{1:41,3:12,5:2,9:1}` | 1 | **1** | **1** |
| 7 | 456 | `{1:399,3:47,5:4,7:1,9:4,21:1}` | 5 | **3** | **3** |

**The SC condition is ESSENTIAL.** At `n=7` the 5 super-symmetric (`|Aut|>7`) classes are: the Paley
`P_7` (`|Aut|=21`) and four `|Aut|=9` classes -- but only **3** of them are SC, and `excess = 3`, not 5.
So "all super-symmetric are SC" is FALSE (n=7); the non-SC super-symmetric classes do NOT force excess.

## The rigorization attempt (complement-pairing x rarity)
Why the SC ones force a covering dimension and the non-SC ones do not:
- **Rarity (the needle):** a class with `|Aut|=a` has `n!/a` labeled reps in the arc-cube `Q_m`
  (`m=C(n,2)`). A super-symmetric class (`|Aut|>n`) has `< n!/n = (n-1)!` reps -- too few to be reliably
  caught by an information-floor-sized subcube. Each such needle beyond the packing capacity costs `+1`
  covering dimension.
- **Complement pairing (why SC forces, non-SC shares):** the complement `T -> T^op` is the ANTIPODAL map
  on `Q_m` (flip all `m` bits, THM-584). A NON-SC class `C` has a distinct partner `C^op`, whose reps are
  the antipodes of `C`'s; a subcube and its antipode (subcube `+` all-ones) share free coordinates, so the
  pair `{C, C^op}` is a single antipodal-linked unit -- one dimension serves both needles. An SC class is
  complement-FIXED (`C = C^op`): it has **no partner**, so it stands alone and forces its own dimension.
- Hence `excess >= #{SC & |Aut|>n}` (each unpaired super-symmetric SC needle `= +1` dimension); empirically
  TIGHT (`=`) for `n<=7`. This is the covering face of the S85 T-join obstruction: the SC classes are the
  fold-fixed rigid spine, and the super-symmetric SC ones are its extreme needles. **PROOF SKETCH, not
  rigorous** (the exact tightness and the `|Aut|>n` threshold need the packing argument made precise).

## Corrections + excess(8)
- **CORRECTION:** my S81 `kappa(n)=1+C(n-2,2)` gives `kappa(7)=11`, but opus/klein found `kappa(7)=12` (the
  lazy-caterer formula BREAKS at `n=7`, the Paley-heptagon obstruction). So `excess(7) = 12 - 9 = 3`.
- **excess(8) = `#{SC & |Aut|>8}` (OPEN):** the order-7 Paley-extended family on 8 vertices is NON-SC (2
  classes, `|Aut|=21`, both `SC=False`), contributing `0`; the SC super-symmetric classes at `n=8` come
  from order-3/5 automorphism structures (full `n=8` enumeration deferred -- `2^21` tilings x `8!` canon).
  Prediction stands as `excess(8) = #{SC & |Aut|>8}`, value pending the harder enumeration.

## The iota-odd certificate field: Q(sqrt-3, sqrt-7) and sqrt21
The `iota`-odd LRC obstruction couples the two apex primes. The natural field is the **biquadratic**
`Q(sqrt-3, sqrt-7)`, `Gal = Z2 x Z2` = **the involution-atlas Klein four** (S88, HYP-3814). Its three
quadratic subfields ARE the three order-2 involutions:
- `Q(sqrt-3)` = **Eisenstein / hexagonal** (`Phi6 = n^2-n+1`, order-6 `n^3=-1`, the Dedekind-sum margin);
- `Q(sqrt-7)` = **apex-7** (Klein quartic `PSL(2,7)`, `X_0(14)`, the genus-1 cusp form `f_14`);
- `Q(sqrt21)` = the REAL coupling, `sqrt-3 * sqrt-7 = -sqrt21`, with **`21 = 3*7 = ` the forbidden `H`
  value `{7,21} = C(7,2)`** -- the TOURNAMENT side.
`h(-3) = h(-7) = 1` (both PID), a genus-1-clean biquadratic matching `X_0(14)` genus 1. So `sqrt21` is the
real geometric mean of the hexagonal (`3`) and apex (`7`) involutions, and its integer square `21` is the
forbidden Hamiltonian-path count -- the covering-min certificate's hexagonal/apex coupling made explicit.
NEXT STEP (contemplated): exhibit `sqrt21` in the covering-min certificate as the entry of the real subfield
where the `E_2`(`sqrt-3`) bulk and the `f_14`(`sqrt-7`) cusp meet.

## Tangential links (assessed)
- **Annals 2026 (Dinur-Evra-Livne-Lubotzky-Mozes), "Good Locally Testable Codes"** (`c^3` LTCs via
  LEFT-RIGHT Cayley complexes): tangential to the flip-rank COVERING CODE (opus HYP-3805) and to the
  two-sided-Cayley = `Z2 x Z2` structure; the Paley/QR tournaments are Cayley graphs. A codes-and-expanders
  lens on the covering excess.
- **Cornell CS 6840 L8 (Tardos), Price of Anarchy via smoothness + the Hotelling game**: the Hotelling
  facility-covering game (clients within radius `d`) is structurally the LRC covering (runners within `r`);
  no-regret dynamics -> CCE echoes kps's "OMWU frequencies = skew spectrum" and the Kaczmarz/POCS witness
  search (S80). A game-theory/learning-dynamics lens on the covering-min.
- **github Pengbinghui/pipeline-math** (AI prover-verifier + Lean-4 formalization of open problems):
  process-relevant -- a template for formalizing the project's LRC14 / tournament results (the repo has a
  Lean LRC14 skeleton).

## Honest scope
`excess = #{SC & |Aut|>n}` is EXACT `n<=7` and the SC-essential fact is confirmed; the rigorization is a
proof sketch (complement-pairing mechanism). `excess(8)` open. The `Q(sqrt-3,sqrt-7)/sqrt21` connection is
exact number theory (biquadratic Klein four, `21=3*7=`forbidden `H`) + a structural tie to the involution
atlas, NOT a certificate. The three links are tangential lenses, not deep dependencies.
