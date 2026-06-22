---
id: HYP-2866
title: LRC14 nu denominator-prefix majorization for Lemma A
status: OPEN / exact bounded evidence; individual q-bucket version refuted
source: codex-2026-06-22-S95
depends_on:
  - THM-530
  - THM-531
  - THM-566
  - HYP-2832
  - HYP-2865
related:
  - THM-567
  - THM-534
  - HYP-2603
  - HYP-2830
  - HYP-+2866
  - HYP-2867
  - OPEN-Q-108
results:
  - 05-knowledge/results/lrc14_nu_denom_center_budget_codex_s95.out
---

# HYP-2866 -- Denominator-Prefix Majorization For The LRC14 `nu` Lemma

## Claim Shape

Let

```text
D(E) = meas{x in [0,1): maxgap({frac(e*x): e in E}) <= 1/7}
nu(E) = 1 - D(E).
```

Lemma A in the Bonferroni route is equivalent to

```text
D(E) <= D(C_k),  C_k={0,1,...,k-1},  8 <= k <= 13.
```

The new candidate proof carrier is not raw runner compression and not a fixed
bounded-denominator witness.  Decompose the exact dense components of `D(E)` by
the least denominator `q` of a reduced rational `p/q` strictly inside the
component.  Write `B_Q(E)` for the total width of components with least
denominator `q <= Q`.

The conjectural useful strengthening is a Landau-style prefix inequality:

```text
B_Q(E) <= B_Q(C_k)          for 7 <= Q <= k,
tail_{q>k}(E) is bounded by the remaining Farey/three-gap width budget.
```

The individual bucket claim `mass_q(E) <= mass_q(C_k)` is false.

## Exact Evidence

The script `04-computation/lrc14_nu_denom_center_budget_codex_s95.py` uses an
exact rational engine.  On each phase-order cell, every cyclic gap is affine;
the dense interval is the intersection of the halfspaces `gap <= 1/7`.  This
avoids the earlier mistake of sampling only collision cells while still giving
the same wall set as the `(7m±1)/(7d)` refinement.

Bounded primitive anchored scan in `[0,14]`:

```text
k=8:  checked 3431, D(C_k)=44/735,     best nonconsec=62/2205
k=9:  checked 3003, D(C_k)=47/294,     best nonconsec=169/1470
k=10: checked 2002, D(C_k)=11/49,      best nonconsec=47/294
k=11: checked 1001, D(C_k)=824/2205,   best nonconsec=1363/4410
k=12: checked 364,  D(C_k)=10432/24255,best nonconsec=824/2205
k=13: checked 91,   D(C_k)=601/1078,   best nonconsec=2851/5390
```

Consecutive is the unique dense-mass leader in all six bounded banks.  The
tightest displayed nonconsecutive gap is at `k=13`:

```text
601/1078 - 2851/5390 = 77/2695.
```

The exact prefix audit found no prefix violations for any `7 <= Q <= k` in the
same banks.  This is stronger and more structured than simply observing that
the total `D` is smaller.

## Negative Findings That Shape The Proof

1. **Raw local compression is false.**  There are one-step-left compression
   dead ends at `k=8`, for example

```text
(0,1,2,3,11,12,13,14): D=40/1617, widths q15=1/49, q23=1/231
(0,1,3,5,7,9,11,13):  D=1/91,   width q15=1/91
```

   These are real local obstructions to a monotone compression proof, but they
   sit far below consecutive and live in high-denominator tails.

2. **Individual q-bucket domination is false.**  Examples from the exact scan:

```text
k=9:  q=9 bucket 41/490 beats C_9's 17/245.
k=10: q=10 bucket 94/1617 beats C_10's 2/49.
k=11: q=11 bucket 3953/28665 beats C_11's 4/105.
k=13: q=13 bucket 2/49 appears while C_13 has no q=13 bucket.
```

   Each such improvement is paid for by losing lower-denominator mass.  The
   surviving invariant is cumulative prefix mass, not pointwise bucket mass.

3. **Fixed bounded-denominator witnesses are not enough.**  THM-566 and
   HYP-2865 show that covering hard cores can kill every denominator up to any
   prescribed `B`.  The denominator object here must be neighborhoods and
   component widths, not a finite list of rational witnesses.

## Tournament Analogy

This is the LRC analogue of the tournament score-sequence lesson.  A single
vertex score can rise under local edge flips, and a nontransitive tournament can
beat the transitive tournament at a selected coordinate.  The robust invariant
is a prefix-majorization ledger, as in Landau's inequalities.

Here the "scores" are denominator bucket widths ordered by increasing `q`.
Nonconsecutive rows can beat a later bucket, but the exact data say they cannot
beat the cumulative low-denominator prefixes.  The proof should therefore
route through a majorization/exchange inequality on component budgets, not
through runner-wise compression.

## Proposed Proof Route

1. **Component existence.**  Use THM-530's denominator-net characterization:
   a dense component near reduced `p/q` requires the residues
   `{e_i p mod q}` to form a strict `1/7`-net.

2. **Local width.**  On a fixed residue pattern, the component width is governed
   by affine gap slack divided by boundary slopes `e_b-e_a`.  Consecutive blocks
   have the smallest slopes among complete low-denominator lifts, explaining
   their wide low-q components.

3. **Prefix exchange.**  Prove that any move increasing a later denominator
   bucket must delete at least as much earlier prefix budget.  This is the
   Landau-score part of the theorem and should be finite over residue patterns
   for `q <= k`.

4. **High-q tail packing.**  Bound `q>k` components by Farey/three-gap width
   packing.  This should connect to the existing EWLB and arc-complexity work:
   bounded component count plus `1/q` local-width decay prevents high-q tails
   from compensating for the missing low-q prefix.

## Impact On LRC14

If the prefix-majorization plus tail-packing theorem is proved, then

```text
D(E) <= D(C_k)
nu(E) >= nu(C_k)
rho*_glob(P,E) >= nu(C_k) + cap_k - 1
```

and HYP-2832's Bonferroni bridge supplies the positive witness floor required
by the current LRC14 skeleton.  This is an alternate route to the p0 cap proof,
but it also clarifies why the p0 and witness routes keep finding the same
extremal AP/consecutive rows.

## Integration With S94 Resonance-Channel Work

After this audit, S94 landed THM-567 and HYP-2867.  That work attacks the
`R'` decorrelation floor through low Fourier resonance channels; this note
attacks Lemma A through low denominator components.  They are different proof
nodes, but they share the same architecture:

```text
finite low channels/components + signed or prefix control + analytic tail.
```

THM-567's QR/NQR even-residue pairing is not a direct proof of this `D(E)`
majorization statement.  Its relevance is diagnostic: it shows that when the
carrier is a residue or channel ledger, the right invariant can be a balanced
aggregate rather than a pointwise bucket.  HYP-2867's correction to the
Paley-7 overclaim is also a warning for this route: do not promote a single
apex-prime bucket to the full proof.  The exact S95 data already says the same
thing on the denominator side: individual q-buckets are false, cumulative
prefixes are the live object.

## Tournament Analysis And Assumption Challenge

Vertices considered:

```text
runners, gaps, fixed 1/7 sections, section boundaries, wall-crossing events,
residues, cover arcs, Fourier modes, rational denominator buckets,
local-slope pairs, matroid-circuit-like relations, and proof obligations.
```

Chosen quotient: denominator-center component budgets plus proof carriers.  It
preserves the exact `D(E)` measure by partitioning it into component widths,
but it destroys runner identity and the full cyclic phase order.

The proof-carrier tournament in the script is transitive:

```text
denominator-center
> EWLB-window
> finite-atlas
> Fourier-spectrum
> runner-compression
> raw-gap-moment
```

Fingerprint: score histogram `{0,1,2,3,4,5}`, no directed 3-cycles, singleton
SCCs, and one Hamiltonian path.

Challenged assumption: the relevant tournament vertices for LRC should be
runners or arcs.  In this problem, runner-level local moves hide the invariant;
the useful object is the sorted denominator-budget ledger.
