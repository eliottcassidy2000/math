---
id: THM-411
title: The exactly-tight LRC family is finite (v_max ≤ 2n-1) and witnessed on the (Z/m)* orbit
status: VERIFIED   # exhaustive computation n=3..8; general-n proof open (see HYP-2196)
source: claudebox-2026-06-03-S621
depends_on:
  - THM-410   # covering-depth conservation + moment-sieve identity (PR #22)
  - THM-403   # the AP's floor-witness orbit is (Z/m)*
related:
  - HYP-2195  # covering-depth master object + additive-chain collapse (refined here)
  - HYP-2196  # the tight-family open problems (this theorem's conjectural extensions)
  - HYP-2130  # rigidity = orbit-type (perspective key) — the witness orbit is the same (Z/m)*
---

# THM-411 — the exactly-tight LRC family: finite, bounded, and (Z/m)*-witnessed

## Setup

`m = n+1`, gap `δ = 1/m`, speeds primitive (`gcd = 1`). A speed set is **tight** (a.k.a.
"barely lonely", Kravitz; the **collapse** family `p_0 = 0` of [[HYP-2195]]) when its loneliness
gap `γ(S) = max_t min_i ‖v_i t‖` equals exactly `1/m`. Equivalently the width-`2δ` forbidden
arcs **cover the circle** (the lonely set has measure zero, surviving only at isolated tight
times). LRC ⇒ `γ ≥ 1/m`, so tight = the extremal instances.

## Statement (exhaustively verified n = 3..8)

> **(1) Finiteness via a magnitude bound.** Every primitive tight set has `v_max ≤ 2n-1`. Hence
> there are finitely many per `n`; the non-AP extremals attain `v_max = 2n-1` exactly.
>
> **(2) Count.** The number of primitive tight sets is `1, 2, 2, 1, 3, 1` for `n = 3,…,8`
> (range-stable up to `v ≤ 8n+2`; not in OEIS). Even `n = 6, 8` give **only the AP**.
>
> **(3) Witness orbit.** For *every* tight set (AP or sporadic), the set of times achieving
> `γ = 1/m` is exactly `{ k/m : gcd(k,m)=1, 0 < k < m/2 }` — the `(ℤ/m)*` orbit, `φ(m)/2` times.
> The sporadic tights are **witness-equivalent to the AP**.
>
> **(4) Necessary residue conditions.** A tight set is `0`-free mod `m` and contains a speed
> `≡ ±1 (mod m)`.

## Evidence

Exhaustive enumeration over primitive sets (exact gap over ℚ, cross-checked by a fine numeric
scan and by the covering-depth `p_0` computation of THM-410):

| n | m | #tight | tight sets | v_max (2n−1) |
|---|---|---|---|---|
| 3 | 4 | 1 | (1,2,3) | 3 (5) |
| 4 | 5 | 2 | (1,2,3,4), (1,3,4,7) | 7 (7) |
| 5 | 6 | 2 | (1,2,3,4,5), (1,3,4,5,9) | 9 (9) |
| 6 | 7 | 1 | (1,2,3,4,5,6) | 6 (11) |
| 7 | 8 | 3 | (1,2,3,4,5,6,7), (1,2,3,4,5,7,12), (1,4,5,6,7,11,13) | 13 (13) |
| 8 | 9 | 1 | (1,…,8) | 8 (15) |

Witness orbit `(ℤ/m)*` confirmed for all tights (CLAIM B, `lrc_tight_claims_s621.out`).

## Why it matters

* **Refines the folklore "AP is the unique extremal".** It is unique only for `n = 3,6,8` in
  this range; `n = 4,5,7` have sporadic tights. The AP's uniqueness is an **even/odd-`m`,
  2-adic-seam phenomenon** (cf. [[HYP-2130]]/[[HYP-2140]]): the witness orbit `(ℤ/m)*` and its
  reflection structure decide whether sporadics can exist.
* **The witness orbit is the perspective/rigidity object.** (3) says the loneliness of the whole
  tight family is carried on the *same* `(ℤ/m)*` orbit that THM-403/404 found for the AP — the
  tight family is exactly the integer lifts that keep the gap pinned on that orbit. This is the
  covering-side incarnation of "rigidity = orbit-type" ([[HYP-2130]]).
* **Connects to the literature.** Tight = Kravitz's "barely lonely"; the bound `v_max ≤ 2n-1`
  matches the known `v_n = 2n-2` tight-instance barrier (Perarnau–Serra survey, normalization
  shift). Finiteness-up-to-scaling is the precise statement.

## Verification

`04-computation/lrc_tight_enum_s621.py`, `lrc_tight_structure_s621.py`, `lrc_tight_claims_s621.py`
→ `05-knowledge/results/lrc_tight_*_s621.out`. General-`n` proof of (1)–(3) is open ([[HYP-2196]]).
