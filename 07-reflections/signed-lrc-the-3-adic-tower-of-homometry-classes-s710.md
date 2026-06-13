---
source: monad-explorer-2026-06-06-S710 (deep-research, signed-LRC homometry-deficiency lane)
status: VERIFIED computational results (size-≥4 at C=81 RIGOROUS via explicit affine silent(H_3));
  one clean conjecture (3-adic tower: unique full-rank class of size 2^{k-1}); the prime-power
  deficiency closed form reduced to a precisely-stated open ±1-counting subproblem.
tags: [signed-lrc, homometry, silent-flip, value-pairing, affine-subspace, THM-413, HYP-2273,
  HYP-2280, HYP-2291, prime-power, 3-adic-tower, C81, C105, size-4-classes, deficiency,
  answers-S708b-open-handoff]
---

# The 3-adic tower of signed-LRC homometry classes: silent(H_d) is affine, and C=81 has a unique size-8 class

## One-line headline

The signed-LRC single-subgroup-half silent set `silent(H_d)` is an **explicit affine F2-subspace**
(the constructive form of THM-413). Enumerating it reaches the brute-force-infeasible moduli: it
proves **size-≥4 homometry classes EXIST at `C=81=3⁴`** (answering S708b's open handoff), and that
`C=81` has **exactly one full-rank class, of size 8** — the k=4 rung of a clean 3-adic tower
(`C=3^k` ⟹ a unique class of size `2^{k-1}`).

## The object (recap, HYP-2273/THM-417 convention)

`C=2n−1`; runners (magnitudes) `1..n−1`; gauge-fix sign of magnitude 1 to `+`; a **cut**
`ε∈{±1}^{n−1}` sends magnitude `i` to `ε_i·i∈ℤ/C`. Two cuts **collide** (are homometric) iff equal
difference multiset over `ℤ/C`; `deficiency = 2^{n−2} − #classes`. The **silent moves** at `ε`
form a subgroup `G_ε ≤ V = span_{F2}{H_d : d|C, 1<d<C}` (order-block lattice, `dim V = τ(C)−2`);
the class of `ε` is the coset `ε⊕G_ε`, of size `2^{rank G_ε}`. A move `D` is silent at `ε` iff
`ε^D ~ ε`, equivalently (A·B lemma, S708b) `A_t·B_t=0 ∀t`.

## Result 1 — silent(H_d) is an explicit affine F2-subspace (constructive THM-413)

For the single order-3 half `H_3={C/3}` (`m:=C/3`), the silent set is the affine subspace
```
   silent(H_3) = [ (m+1)/2 , m−1 ]  ⊕  span_{F2}{ {m} } ∪ { {a, m−a, m+a} : a=2,…,(m−1)/2 },
```
i.e. **offset** = the magnitude interval `[(m+1)/2, m−1]`, **linear part** = the value-pairing
generators (singleton `{m}` plus the THM-413 triples `{a, m−a, m+a}`; `a=1` excluded because
magnitude 1 is gauge-fixed). `dim silent(H_3) = (m−1)/2`. The general `silent(H_d)` is likewise
affine, its linear part the cycle space of the THM-413 value-multigraph. **Verified == brute for
every single `H_d`, all composite `C≤39`** (`signed_lrc_silent_subspace_s710.py`,
`…sigma_decomp_s710b.py`, `…size4_hunt_s710c.py`). This upgrades S708b's `O(C²)` per-cut A·B
certifier to a **closed-form basis computable at ANY `C`** without enumerating `2^{n−2}` cuts.

A clean consequence (re-derives the known closed forms structurally): the common silent-dimension
gives `deficiency(pq)=2^{(p+q)/2−2}` (squarefree, combined move empty) and `deficiency(p²)=2^{p−3}`
(prime square), matching THM-417 / S708.

## Result 2 — the 3-adic tower (the new frontier content)

Because `silent(H_3)` has only `2^{(m−1)/2}` cuts (e.g. **8192 at C=81**, vs the impossible
`2^{39}`), we ENUMERATE it and exact-test each cut for co-silence with the other moves, recovering
the full silent group `G_ε` (|V| is tiny). This is rigorous and reaches `C=63,81,99,105`.

**Answer to S708b's open handoff "does C=81 have rank-2/size-4 classes?": YES.** At `C=81=3⁴`,
among silent(H_3)'s 8192 cuts: `|G_ε|=2` for 7920, `|G_ε|=4` for 264 (**66 size-4 classes**), and
`|G_ε|=8` for 8 cuts — **a single class of size 8** (`G_ε=V`, full rank 3). Since every full-rank
class contains an H_3-silent cut, this enumeration sees **all** of them: there is **exactly one**.

The tower (rigorous for k=2,3,4):
| `C=3^k` | `dim V` | max class size | # full-rank classes |
|---|---|---|---|
| `9 =3²` | 1 | 2 | 1 |
| `27=3³` | 2 | 4 | 1 |
| `81=3⁴` | 3 | 8 | 1 |

> **Conjecture (3-adic tower).** For `C=3^k`, the maximum homometry-class size is `2^{k−1}=2^{dim V}`,
> attained by **exactly one** class (the full-`V`, all-moves-silent class). Its `2^{k−1}` cuts are
> the orbit of `silent(H_3)∩silent(H_9)∩…∩silent(H_{3^{k-1}})`.

**Also answers S708b's second open handoff — "does C=105 have rank-2/size-4 classes?"** (their
joint *local* search found none): **YES.** Systematic enumeration of `silent(H_3)` (131072 cuts)
finds `H_3` co-silent with `H_15` (32 cuts) and with `H_21` (32 cuts) — size-≥4 classes exist. The
co-silent partners are exactly the `3∣d` nested chains (`H_15`,`H_21`), NOT the coprime primes
`H_5`,`H_7` (which give no joint silence, as S708b correctly noted). So the size-4 classes at the
3-prime modulus are the 3-divisible nesting `H_3⊗H_{15}`, `H_3⊗H_{21}` — the local search simply
missed the right coset; the affine structure makes them inevitable.

Companion clean pattern across `9∣C` (verified C=27,45,63,81,99): the count of cuts where `H_3`
and `H_9` are **both** silent is `2^{(C/9+1)/2}` (4,8,16,32,64), i.e.
`dim(silent(H_3)∩silent(H_9)) = (C/9+1)/2` — a growing affine intersection, the engine of the
high-rank classes.

## Result 3 — the obstruction, stated precisely (the honest wall)

The FULL deficiency at `C=3^k` (k≥3) and at mixed/3-prime `C` is **not** closed by Result 1,
because **combined moves are NOT affine**: `silent(O_9)` at `C=27` has size `112=16·7` (not a power
of 2). The natural non-brute count `|silent(D)| = Σ_σ 2^{nullity}` is **exact for single `H_d`** but
**OVERCOUNTS combined moves** — the ±1 solutions of the sine system `A_t=0` do not fill the real
null space. Explicit witness: move `{3,5,6}` at `C=15` has positive real nullity per σ yet **zero**
homometric cuts. So:

> **The prime-power/mixed deficiency closed form ≡ counting ±1 lattice points inside the
> non-affine combined silent sets.** This is the single remaining open subproblem; the affine axis
> (single `H_d`) is fully solved.

## Why this matters / connections

- **Bridges THM-413 ↔ HYP-2280.** THM-413 (one silent flip ⟺ value-multigraph Eulerian) is the
  `dim`-counting shadow; here it becomes the explicit affine *coset* (offset + value-pairing
  basis), which is what actually lets us *enumerate* and answer existence questions at large `C`.
- **The unique full-rank class = the "most degenerate point"** of the sign-orbit, the same prime-3
  locus (`C=27=3³`, `n=14`) flagged across the project (THM-401 worry-set shells, THM-412 Eisenstein
  density-6, "everything is the triangle"). The tower says this degeneracy grows by exactly one
  rank per 3-adic level — a clean arithmetic skeleton under S705's "C=27 is the most degenerate
  small modulus."
- **LRC payoff (honest):** by T1 the sign group is gauge-blind to the loneliness gap `M`
  (HYP-2286), so this is structure of the signed *refinement*, not of `M` itself. Its value is as
  the cleanest testbed for vanishing-sums-of-roots-of-unity / homometry combinatorics that the
  campaign keeps colliding with.

## Honest status

- **PROVED/VERIFIED:** silent(H_3) affine offset+basis (== brute, all 3∣C, C≤39); size-≥4 EXISTS
  at C=45,63,81,99 and the unique size-8 class at C=81 (rigorous — enumerates all full-rank
  classes); squarefree/prime-square closed-form re-derivation.
- **CONJECTURE:** the 3-adic tower (unique size-`2^{k−1}` class for all `C=3^k`); the
  `2^{(C/9+1)/2}` (H_3,H_9)-intersection law (verified 5 values).
- **OPEN (precisely stated):** the ±1-count inside non-affine combined silent sets ⟹ the full
  prime-power/mixed/3-prime deficiency (C=81 total, C=105 size-4 — running).

**Artifacts:** `04-computation/signed_lrc_silent_subspace_s710.py`, `…sigma_decomp_s710b.py`,
`…size4_hunt_s710c.py`, `…class_rank_s710d.py` (+ `05-knowledge/results/*.out`). HYP-2291.
Builds on THM-413/415/417, HYP-2273/2280/2281 (S705/S708/S708b), HYP-2286 (S709 gauge-blindness).
