---
id: THM-1475
title: "THE PFAFFIAN IS THE ODD FUNCTION, AND ITS SPECTRUM HAS GAPS. (1) Pf(S(T)) is ODD for every tournament of even order — one line: S ≡ J−I (mod 2), so Pf(S) ≡ Pf(J−I) = #{perfect matchings of K_n} = (n−1)!!, a product of odds. This is the general fact behind kind-pasteur's THM-1455 observation that 4-subsets have |Pf| ∈ {1,3}: det = Pf² ≤ 9 plus oddness forces the dichotomy. (2) Pf is literally an ODD (alternating) function and det = Pf² literally an EVEN (invariant) one — Pf(PᵀSP) = sign(σ)·Pf(S) — so the owner's odd-function/even-function dichotomy is carried by the sign character of S_n. (3) |Pf| is a SWITCHING-CLASS invariant, since Pf(DSD) = det(D)Pf(S); hence |Pf| is a function on the 2^{C(n−1,2)} TILINGS, not on all 2^{C(n,2)} tournaments, and descends to the A049313 classes. (4) COROLLARY: Pf odd ⟹ Pf ≠ 0 ⟹ Aut(T) ≤ A_n for n even — an independent one-line route to a classical consequence of Moon. (5) THE SPECTRUM, EXHAUSTIVE: {1,3} at n=4 and {1,…,9} at n=6, both GAP-FREE, but at n=8 the attained set is {1,…,27} ∪ {31,33,35} ∪ {49} — GENUINE GAPS at 29 and at 37,39,41,43,45,47. The maximum 49 = 7² is the skew-Hadamard value (SSᵀ = 7I ⟹ det = 7⁴) and it is ISOLATED, sitting above an eleven-wide run of unattainable values. So Pf joins hp as an odd-valued tournament invariant with a non-trivial spectrum — hp omits {7,21}, Pf omits {29,37,…,47} at n=8. (6) NEGATIVE: hp does NOT determine |Pf| mod 4, 8 or 16."
status: >
  (1) PROVED (one line) + VERIFIED n = 4,6,8.
  (2) PROVED (classical Pf(BᵀSB) = det(B)Pf(S)) + VERIFIED on 400 random relabellings.
  (3) PROVED (same identity with B = D diagonal ±1) + VERIFIED on 150 random switchings.
  (4) PROVED; the conclusion itself is classical via Moon's odd |Aut(T)|, the ROUTE is new.
      Verified: 810 automorphisms found at n = 4,6, none of odd sign.
  (5) EXHAUSTIVE AND DECISIVE.  n = 4,6 by full enumeration; n = 8 by full enumeration of
      all 2^21 = 2,097,152 switching classes, which is exhaustive for |Pf| BECAUSE of (3).
      The gaps are not search failures.
  (6) VERIFIED-EXHAUSTIVE at n = 6 (13 joint residue classes mod 8).  The apparent
      determination at n = 4 is a THREE-POINT small-case artefact and is flagged as such.
  SELF-CORRECTION ON RECORD: an earlier run of this session printed a VERDICT asserting the
  gaps were sampling artefacts.  That text was pre-written and printed regardless of the
  data, which in fact refuted it.  Retracted in place; see §5.
source: klein-2026-07-20-S339 (owner: each odd sized tournament corresponds to one of the natural numbers; chase the high leverage question; see the relation between odd valued functions and tournament adjacent ideas; they both relate also to even concepts like even graphs and even functions)
depends_on:
  - THM-1470  # klein: Pf(S(T)) odd, even tournaments, the skew-Seidel frame
related:
  - THM-1455  # kind-pasteur: the LOCAL Pfaffian expansion and the mod-16 law that (1) underwrites
  - THM-1440  # (seidel-spectra) the canonical Cut+Cycle splitting at ODD n — the "tournament = natural number" half
  - THM-476   # the skew Ehlich–Wojtas square law (the maximum, not the spectrum)
  - THM-343   # the hp-spectrum {odds} ∖ {7,21}
script: 04-computation/pfaffian_spectrum_klein_S339.py (+ .out)
---

# THM-1475 — the Pfaffian is the odd function, and its spectrum has gaps

## 0. What was already settled elsewhere

The owner's "each odd sized tournament corresponds to one of the natural numbers" is the
**canonical `Cut ⊕ Cycle` splitting at odd `n`** — the bicycle space `Cut ∩ Cycle` of `K_n` is
`0` exactly when `n` is odd, so `T ↦ (`tiling number`, `score-parity vector`)` is a canonical
bijection. That was established this same day by another agent's THM-1440 (seidel-spectra) and
is **not** re-derived here. This file works the *other* half of the directive: the odd-valued
functions, and the odd/even function dichotomy.

## 1. `Pf(S(T))` is odd

For `n` even, with `S` the skew adjacency (`S_ij = ±1`, `Sᵀ = −S`):

> Every off-diagonal entry is `±1 ≡ 1 (mod 2)`, so `S ≡ J − I (mod 2)`, hence
> `Pf(S) ≡ Pf(J−I) = #\{`perfect matchings of `K_n\} = (n−1)!!`, a product of odd numbers.
> **`Pf(S(T))` is odd — in particular never zero.** ∎

This is the general fact behind kind-pasteur's THM-1455, whose mod-16 law is stated in terms of
`k₄ = #\{`4-subsets with `|Pf| = 3\}`. For a 4-subset, `det = Pf² ≤ 9`, and oddness leaves only
`|Pf| ∈ \{1,3\}` — **the dichotomy their argument uses is forced, not observed.**

## 2. Odd function, even function

`Pf(BᵀSB) = det(B)·Pf(S)`. Two specialisations carry the owner's dichotomy exactly:

| `B` | consequence | character |
|---|---|---|
| permutation matrix `P_σ` | `Pf(P_σᵀ S P_σ) = sign(σ)·Pf(S)` | **odd / alternating** |
| — | `det S = Pf(S)²` is invariant | **even** |

So `Pf` is a *relative* invariant of weight `sign`, and `det` is a genuine invariant: the
odd-function/even-function pair, realised on tournaments by the sign character of `S_n` rather
than by `x ↦ −x`. *(Verified on 400 random relabellings at `n = 4,6`.)*

## 3. `|Pf|` is a switching-class invariant

Taking `B = D` diagonal `±1` gives `Pf(DSD) = det(D)·Pf(S) = ±Pf(S)`, so:

> **`|Pf|` is constant on Seidel switching classes.** It is therefore a function on the
> `2^{C(n−1,2)}` **tilings** — not on all `2^{C(n,2)}` tournaments — and descends to the
> `A049313` classes.

*(Verified on 150 random switchings.)* This is what makes §5 exhaustive at `n = 8`: `2^21`
switching classes instead of `2^28` tournaments.

## 4. Corollary: `Aut(T) ≤ A_n` for `n` even

An automorphism `σ` gives `Pf(S) = sign(σ)·Pf(S)`. If `sign(σ) = −1` then `Pf(S) = 0`,
contradicting §1. So **no tournament of even order has an automorphism of odd sign.**

The conclusion is classical — Moon's theorem that `|Aut(T)|` is odd gives it for every `n`,
since an element of odd order has sign `+1`. The **route** is new: it comes from the oddness of
a single integer. *(Verified: 810 automorphisms across `n = 4,6`, none of odd sign.)*

## 5. The spectrum — and the isolated maximum

| `n` | attained `\|Pf\|` | gaps below the max |
|---|---|---|
| 4 | `{1, 3}` | **none** |
| 6 | `{1, 3, 5, 7, 9}` | **none** |
| 8 | `{1,3,…,27} ∪ {31, 33, 35} ∪ {49}` | **`29`, and `37,39,41,43,45,47`** |

`n = 4, 6` are full enumerations; `n = 8` is a full enumeration of all `2^21 = 2,097,152`
switching classes, which by §3 is exhaustive for `|Pf|`.

> **The Pfaffian spectrum is gap-free at `n = 4, 6` and acquires gaps at `n = 8`.**

Two features worth naming:

- **The maximum is isolated.** `49 = 7²` is the skew-Hadamard value: a doubly regular tournament
  on 8 vertices has `SSᵀ = 7I`, so `det S = 7⁴` and `|Pf| = 49`. It sits above an *eleven-wide*
  run of unattainable values `37…47`. Nothing approaches the optimum gradually — this is a
  rigidity phenomenon at the extremal design, not a boundary effect. THM-476 governs the
  *maximum*; the *gap below it* is new here.
- **A second odd-valued invariant with a non-trivial spectrum.** `hp` attains `{odds} ∖ {7,21}`
  (THM-343/THM-115). `Pf` now joins it. Whether the `n = 8` gaps persist or close as `n` grows
  is open, and is the natural next computation (`n = 10` needs `2^36` switching classes, so it
  needs a better idea than enumeration).

## 6. The negative: `hp` does not control `Pf`

Two odd-valued invariants on the same object invite a congruence. There is none:

- `n = 6`, exhaustive: the joint `(hp, |Pf|)` residues occupy **13** classes mod 8, and `hp`'s
  residue does **not** determine `|Pf|`'s mod 4, 8 or 16.
- `n = 4` *appears* to determine it — but the entire joint support is `{(1,1),(3,3),(5,1)}`,
  three points. Flagged as a small-case artefact, not evidence.

## 7. Method note, on the record

An earlier run in this session printed a "VERDICT" declaring the `n = 8` gaps to be sampling
artefacts. That sentence was **written before the data existed** and printed unconditionally;
the run's own output — targeted hill-climbing failing on exactly `29` and `37…47` while
reaching `31,33,35` and `49` — refuted it. The lesson is narrower than "check your work": *do
not put a conclusion in a print statement*. A script should report what it measured and leave
the verdict to be written afterwards.

*Files: `04-computation/pfaffian_spectrum_klein_S339.py` (+ `.out`).*
