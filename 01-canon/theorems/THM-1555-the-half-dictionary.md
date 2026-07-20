---
id: THM-1555
title: "THE HALF DICTIONARY — [PRIORITY: parts (0) and (II) are klein's THM-1560, pushed 23 minutes earlier and stated more strongly — cite theirs; the net-new content here is (I), the spectral trichotomy, and (III), the charge/Newton-midpoint identification] — {−1,0,1} and {0,½,1} are one trichotomy in two coordinates, conjugate by x ↦ (1+x)/2, and the residue is always the ½. (0) THE CONJUGATION IS LITERAL, not analogical: for every tournament A = (J − I + S)/2 exactly, S the skew matrix with entries in {−1,0,1}, A the adjacency matrix with entries in {0,1} — verified 5096/5096 — and the −I/2 is the residue left by the diagonal, since (1+0)/2 = ½ rather than 0. (I) DIFFERENTIATING THE DICTIONARY GIVES A + Aᵀ = J − I, AND THAT ONE IDENTITY IS USED THREE TIMES. (I.a) on a general eigenvector: 2·Re(x*Ax) = |Σxᵢ|² − 1 ≥ −1, so EVERY eigenvalue of EVERY tournament has Re(λ) ≥ −½, equality iff the eigenvector sums to zero; minimum over 4096 tournaments attained at exactly −0.500000000000. (I.b) on the PERRON vector: 2ρ = (1ᵀx)² − 1 ≤ n − 1, so ρ ≤ (n−1)/2 with equality iff x is uniform iff T is REGULAR — 0 identity failures, and max ρ at n = 4,6,7 is 1.3953, 2.4340, 2.9294, strictly below (n−1)/2 exactly where no regular tournament exists. (I.c) with tr A = 0 the two ends CLOSE INTO AN EQUIVALENCE: T regular ⟺ ρ = (n−1)/2 ⟺ every non-Perron eigenvalue lies on the vertical line Re = −½; and the failure is an exact identity, not an inequality — Σ(Re λₖ + ½) = (n−1)/2 − ρ = (n − (1ᵀx)²)/2, so spectral height above the line = Perron deficiency = the Perron vector's deviation from uniform, one number with three readings. (II) THE CHARACTERISTIC-2 COLLAPSE EXPLAINS WHY THIS REPO CARRIES TWO MATRIX FORMALISMS: the dictionary is invertible iff 2 is, and mod 2 it fails ASYMMETRICALLY — S ≡ A + Aᵀ ≡ J + I mod 2 is INDEPENDENT OF T, so S mod 2 takes exactly ONE value at every n (verified n = 3..7), while A mod 2 separates everything. The GF(2) layer (cut/cycle, switching classes, even graphs, Rédei mod 2) MUST be written in A; the parity layer (Pfaffians, skew spectrum, THM-1455) MUST be written in S. Neither is a convention: each is unavailable on the other side. (III) THE SAME DICHOTOMY IN THIS WEEK'S GMC WORK: THM-1540's nullcone condition is that the charges a−b are all of one strict SIGN, while its step L2 asks whether the support of h straddles the Newton midpoint d/2. Since a − b = 2a − d these are one condition in two coordinates — agreeing on every case tested. Charge is the sign coordinate, exponent-with-midpoint is the affine coordinate, and the d/2 IS the ½"
status: >
  (0) VERIFIED-EXACT over 5096 tournaments, n = 3..6; algebraic, the check only
  confirms no sign convention slipped.
  (I.a) PROVED in one line from A + Aᵀ = J − I, equality case immediate.  CLASSICAL —
  presented as a recovery, not a discovery.
  (I.b) PROVED (Perron–Frobenius + Cauchy–Schwarz).  Also classical for tournament
  matrices; recovered here to show it is the SAME identity as (I.a) on a different
  vector.
  (I.c) PROVED both directions and self-contained: (⟸) trace 0 forces ρ = (n−1)/2,
  then Cauchy–Schwarz equality forces the Perron vector uniform, hence regular;
  (⟹) for regular T, 1 is both a left and a right eigenvector so 1^⊥ is A-invariant
  and A = (S − I)/2 there, and a real skew matrix has purely imaginary spectrum.  Note
  this needs only skew-symmetry, NOT THM-1455's parity refinement — the earlier draft
  overstated that dependency and it is corrected here.  The deficiency identity is
  immediate from tr A = 0.  Elementary throughout; the content is the consolidation.
  (II) PROVED in one line (−1 ≡ 1 mod 2) and VERIFIED n = 3..7: S mod 2 has exactly
  one value at each n.  NOT MINE — this is klein's THM-1560(B), reached independently
  and 23 minutes later, and theirs is the stronger form.  See the PRIORITY note below.
  (III) VERIFIED on the tested supports; a − b = 2a − d is algebraic.
  DETECTION FLOOR for the computations: exhaustive at n ≤ 5 (8, 64, 1024 tournaments);
  at n = 6, 7 the first 3000–4000 of 2^15 and 2^21 only.  Since (I) and (II) are proved,
  the computation is confirmation and not the evidence base.
  PRIORITY: klein-2026-07-20-S349 pushed THM-1560 covering (0) and (II) at 12:06:14;
  this file went in at 12:29:26.  Independently derived, but klein is first-pusher and
  their (B) is the stronger statement (every integer polynomial in S is constant mod 2,
  with THM-1475's Pfaffian-oddness as a corollary).  Those sections are retained only
  in compressed form, with the deferral stated in the body.  NET-NEW HERE: (I) and (III).
  This theorem settles no open problem.  It states a dictionary, proves what the
  dictionary implies, and shows two results from this session series are one result.
source: kind-pasteur-2026-07-20-S128c119 (owner: {−1,0,1} and {1,½,0} are functionally equivalent and both recur here; this is what I mean by even/odd versus positive/negative)
depends_on:
  - THM-1540    # the GMC(2) nullcone one-sidedness condition
related: [THM-1560, THM-1455, THM-1525, THM-895, THM-1440]
supersedes_priority_note: THM-1560 (klein-S349) owns the dictionary and the mod-2 collapse
script: 04-computation/the_half_dictionary_kps_S128c119.py, half_dictionary_consequences_kps_S128c119.py (+ .out)
---

# THM-1555 — the half dictionary

The owner's observation is that `{−1, 0, 1}` and `{1, ½, 0}` keep recurring here and are
"functionally equivalent". They are, by

> **`x ↦ (1 + x)/2` :  `−1 ↦ 0`,  `0 ↦ ½`,  `+1 ↦ 1`.**

The left set is the **sign** world — multiplicative, closed under negation, the odd/even
axis. The right is the **affine** world — convex, living in `[0,1]`, with `½` the
midpoint. The point of this note is that the repo's objects come in conjugate *pairs*
under exactly this map, and that **the number `½` is where the conjugation leaves its
residue**.

## 0. The conjugation is literal — but klein got here first

A tournament has two matrices: the skew `S` with entries in `{−1, 0, +1}` and the
adjacency `A` with entries in `{0, 1}`. They are related by the owner's map, exactly:

> **`A = (J − I + S)/2`,  `S = 2A − J + I`** — verified on 5096 tournaments, no failures.

And the `−I/2` is not bookkeeping: it is the residue left by the diagonal, because
`(1 + 0)/2 = ½` and the adjacency diagonal must be `0`.

> ### ⚠ PRIORITY: this section and §II below are **klein's THM-1560**, not mine.
>
> The owner dispatched this prompt to several machines. **klein-2026-07-20-S349 pushed
> THM-1560 ("the halving dictionary and the mod-two collapse") at 12:06:14; I pushed at
> 12:29:26 — klein is first-pusher by 23 minutes.** Their (A) is this section and their
> (B) is my §II, derived independently and *stated more strongly*: where I show only that
> `S mod 2` is constant, klein proves that **every invariant which is an integer
> polynomial in `S` is constant mod 2**, and derives THM-1475's Pfaffian-oddness as a
> one-line corollary. Their (D) then makes the point I missed entirely — that Rédei's
> `hp` is mod-2 blind *while needing the ½*, so its oddness is a genuinely special fact
> and not an instance of the collapse.
>
> **Cite THM-1560 for the dictionary and the collapse.** What remains net-new here is
> §I (the spectral trichotomy) and §III (the charge/Newton-midpoint identification).
> §II is retained in compressed form only because its verification counts differ from
> klein's and because §I depends on the identity it derives.

## I. Differentiate the dictionary and use the result three times

Adding the dictionary to its transpose kills `S` and leaves

> **`A + Aᵀ = J − I`.**

That single identity does all the spectral work below. What changes between the three
uses is only **which vector it is evaluated on**.

### I.a — on a general eigenvector: a half-plane

For any unit `x`, `2·Re(x*Ax) = x*(J − I)x = |Σᵢ xᵢ|² − 1 ≥ −1`. Hence

> **every eigenvalue of every tournament satisfies `Re(λ) ≥ −½`**,

with equality exactly when the eigenvector sums to zero. Verified over 4096 tournaments:
the minimum real part observed is **exactly `−0.500000000000`**, so the bound is attained.

### I.b — on the Perron vector: a ceiling on `ρ`

Let `x ≥ 0` be the unit Perron vector, `Ax = ρx`. Then `x*Ax = x*Aᵀx = ρ`, so the same
identity reads `2ρ = (1ᵀx)² − 1`, and Cauchy–Schwarz gives `(1ᵀx)² ≤ n‖x‖² = n`:

> **`ρ ≤ (n−1)/2`, with equality iff `x` is uniform iff `T` is regular.**

Verified with zero identity failures. The max `ρ` over all tournaments is `1, 1.3953, 2,
2.4340, 2.9294` for `n = 3..7` — attaining `(n−1)/2` exactly at odd `n`, where regular
tournaments exist, and falling strictly below at even `n`, where they do not.

Both I.a and I.b are classical; they are recovered here to exhibit them as **one
identity read on two vectors**.

### I.c — with `tr A = 0` the two ends close into an equivalence

`tr A = 0` gives `0 = ρ + Σ_{k≥1} Re λ_k`. Feeding in I.a recovers `ρ ≤ (n−1)/2` a
second, independent way — and running the chain backwards closes it:

> **`T` is regular  ⟺  `ρ = (n−1)/2`  ⟺  every non-Perron eigenvalue lies on the
> vertical line `Re = −½`.**

*Proof.* (⟸) If the `n−1` non-Perron eigenvalues all have `Re = −½`, then `tr A = 0`
forces `ρ = (n−1)/2`; then equality in I.b forces the Perron vector uniform, so all row
sums equal `(n−1)/2`. (⟹) If `T` is regular then `1` is both a left and a right
eigenvector, so `1^⊥` is `A`-invariant; there `J = 0` and the dictionary gives
`A = (S − I)/2`; a real skew matrix has purely imaginary spectrum, so those eigenvalues
are `−½ + it/2`. ∎

Only skew-symmetry is used, not THM-1455's parity refinement — an earlier draft of this
file claimed that dependency and it was too strong.

**Prior art, credited:** THM-1440(D) already established the line `Re = −½` **for
circulant tournaments**. What is added here is that it holds for *every* regular
tournament, that it is an *equivalence* (the line characterises regularity), and that
the non-regular case is governed by the exact deficiency identity below. THM-1440(D) is
the special case; this is not a rediscovery of it.

**And the failure of regularity is an exact identity, not an inequality:**

> `Σ_{k≥1} (Re λ_k + ½)  =  (n−1)/2 − ρ  =  ( n − (1ᵀx)² ) / 2`.

The **total height of the spectrum above the line** equals the **Perron deficiency**
equals the **Perron vector's deviation from uniform**. One number, three readings;
verified to `10⁻¹⁴` at every `n` tested.

Two further constraints come free and cost nothing: `tr A = 0` and `tr A² = 0` (because
`A_ij A_ji = 0` in a tournament), so `Σλ = Σλ² = 0` before any structure is imposed. The
dictionary supplies the third constraint, and it is of a different type — a **half-plane**
rather than a moment. Regularity is exactly the case of equality throughout.

## II. The characteristic-2 collapse — see THM-1560(B), which is stronger

Compressed to what §I needs, with priority to klein.

The dictionary is invertible exactly when `2` is. Over `GF(2)` it fails, and it fails
**asymmetrically**: since `−1 ≡ 1`, `S = A − Aᵀ ≡ A + Aᵀ = J − I (mod 2)`, which **does
not depend on the tournament at all**. My verification, `n = 3..7`: `S mod 2` takes
**exactly one value at every `n`**, while `A mod 2` separates all `8, 64, 1024, …`
tournaments tested — a distinctness count, where klein checked the identity itself.

So the `GF(2)` layer (cut/cycle, switching classes, even graphs, Rédei mod 2) must be
written in `A`, and the char-0 parity layer (Pfaffians, skew spectrum, THM-1455,
THM-1475) in `S`. Neither is a convention that could have gone the other way.

**THM-1560(B) generalises this from `S` to every integer polynomial in `S`, which is the
right statement; use theirs.** The only methodological remark I would add is that it
gives a cheap consistency check: any attempt to run a parity or Pfaffian argument through
a mod-2 reduction of `S` is guaranteed vacuous before it starts.

## III. The same dichotomy inside this week's GMC work

THM-1540 reduces GMC(2) to: the charges `a − b` are all of one strict **sign**. Its proof
step L2 works with `h(t)·t^{−d/2}` and asks whether the support straddles the **Newton
midpoint** `d/2`. Since

> `charge = a − b = 2a − d`,

those are one condition written twice: *"all charges of one sign"* and *"the midpoint
`d/2` lies outside the support"*. Verified to agree on every support tested:

| support of `h` | `d` | charges `2a − d` | one-sided? | midpoint outside? |
|---|---|---|---|---|
| `{1}` | 1 | `{1}` | yes | yes |
| `{0}` | 1 | `{−1}` | yes | yes |
| `{0,1}` | 1 | `{−1,1}` | no | no |
| `{2}` | 2 | `{2}` | yes | yes |
| `{0,1,2}` | 2 | `{−2,0,2}` | no | no |
| `{1,2,3}` | 3 | `{−1,1,3}` | no | no |

The charge is the sign coordinate; the exponent-with-midpoint is the affine coordinate;
and the `d/2` *is* the `½`.

## What this is and is not

**It settles no open problem, and most of the individual statements are classical.** The
half-plane `Re ≥ −½` and the ceiling `ρ ≤ (n−1)/2` are both standard facts about
tournament matrices and are recovered, not discovered. What is new here is the
consolidation and two of its consequences:

- the three spectral facts of §I are **one identity evaluated on three vectors**, and
  that is why the trichotomy in I.c closes with no extra input;
- the characteristic-2 collapse of §II is, as far as this repo's canon goes,
  unrecorded — and it *explains* an otherwise arbitrary-looking feature of the corpus;
- THM-1540's nullcone one-sidedness and the spectral line are the same dichotomy in
  the two coordinates.

The practical value is directional. When a statement here is awkward in one coordinate,
write down the conjugate: parity and divisibility arguments want the sign world
`{−1,0,1}`; convexity, Newton polytopes, measures and thresholds want the affine world
`{0,½,1}`. THM-1455's `k₄ ≡ 0 mod 2` argument is natural on the left; THM-1540's
midpoint-of-a-segment argument is natural on the right; each was hard to see from the
other side.

Note also what §II says about *method*: any attempt to run a parity or Pfaffian argument
through a mod-2 reduction of `S` is guaranteed to be vacuous. That is a cheap
consistency check worth applying to future skew-matrix work.

## Named next

- **The even-`n` deficiency.** At even `n` no tournament attains `ρ = (n−1)/2`, and the
  observed maxima `1.3953 (n=4)`, `2.4340 (n=6)` are the near-regular optima. By I.c the
  deficiency `(n−1)/2 − ρ` equals the spectral height above the line — so *"how close to
  regular can an even-`n` tournament be"* is literally *"how close to the line can the
  spectrum sit"*. Worth a closed form; the `n = 4` value is an algebraic number of small
  degree and should be identified.
- Three further `½`'s in the repo look like the same residue and were **not** checked:
  the `(u − iv)/2` in the de Bondt–van den Essen twist (THM-1430) — which is likely the
  **polarisation** constant, the other place characteristic 2 obstructs — the centres at
  permutations of `(¼, ½, ¾)` in the Kakeya bad set, and the incentre slacks `1/28`.
- The `GF(2)` table in §II predicts that every repo result stated in `S` should have an
  `A`-form that survives mod 2 and vice versa. Auditing THM-1440 and THM-1475 against
  that would either produce translations or locate a genuine obstruction.
