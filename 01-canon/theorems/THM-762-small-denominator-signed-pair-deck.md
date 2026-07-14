---
id: THM-762
title: The small-denominator signed-pair deck criterion (15 <= q <= 28), and exact counterexamples to the proposed uniform q <= 25 good-period bound
status: PROVED (criterion and explicit counterexamples); the historical S105 q <= 38 statement is FINITE-EXACT FOR THAT CAPPED BANK ONLY
source: codex-2026-07-14-S3
depends_on:
  - THM-366
related:
  - THM-566
  - THM-758
  - THM-760
  - HYP-6820
  - MISTAKE-143
---

# THM-762 — the small-denominator signed-pair deck

This theorem replaces the false reading of the S312 observation that every
covering residual has a lonely rational time of denominator at most `25`.
It identifies the exact finite object seen by denominators `15,...,28` and
exhibits both coherent and gcd-incoherent counterexamples to the proposed
uniform bound.

## 1. Exact criterion

Let `S` be a finite set of positive integer speeds which is **covering**:
for every `d=2,...,14`, some `s in S` is divisible by `d`.  For an integer
`q` with `15 <= q <= 28`, write

```
U_q = (Z/qZ)^x / {+1,-1},
B_q(S) = { [s]_+- : s in S and gcd(s,q)=1 } subseteq U_q,
Z_q(S) = {s in S : q divides s}.
```

Here `[s]_+-` denotes the signed unit class `{s,-s}`.  Then the following
are equivalent:

1. there is an integer `a`, `1 <= a < q`, for which
   `||as/q|| >= 1/14` for every `s in S`;
2. `Z_q(S)` is empty and `B_q(S)` is a proper subset of `U_q`.

Thus a `q`-witness exists exactly when there is no zero owner and the signed
unit-pair blocker deck misses a card.

### Proof

First suppose `gcd(a,q)=g>1`, and reduce `a/q=a_0/q_0`.  Since
`q<=28`, we have `q_0=q/g<=14`.  Covering gives an `s in S` divisible
by `q_0`, so `a_0s/q_0` is an integer.  Hence every witness in this range
must have `a` a unit modulo `q`.

For a unit `a`, put
`r_s=min(as mod q, -as mod q)` in `{0,...,floor(q/2)}`.  Because
`15<=q<=28`, the strict failure inequality `r_s<q/14` is equivalent to
`r_s in {0,1}`.  The value `0` occurs exactly when `q|s`.  If there is no
zero owner, the value `1` occurs exactly when `s` is a unit and
`[s]_+-=[a^{-1}]_+-`.  Consequently `a` is a witness exactly when there
is no zero owner and the inverse signed class is absent from `B_q(S)`.
As inversion permutes `U_q`, such an `a` exists exactly when the blocker
deck is incomplete. ∎

The endpoint `q=28` is harmless: `r_s<2`, so the only failing integer
distances are still `0` and `1`.  For `q>=29`, distance `2` can also fail;
the signed-pair deck alone no longer contains the full predicate.

## 2. Uniform `q<=25` is false

Two exact thirteen-speed counterexamples are

```
V26 = {26,52,78,104,130,156,182,208,234,260,286,312,339},
S*  = {81,91,131,151,157,196,258,274,313,328,330,339,348}.
```

Both are primitive, covering, have thirteen speeds above `14`, diameter at
most `339`, and pass the exact rational analogue of every leave-one-out
capped-envelope test used by S312.  Their zero-owner/pair decks are complete
for every `q=15,...,25`; covering blocks every `q<=14`.  Therefore neither
has a witness with denominator at most `25`.

The examples rule out two possible repairs:

- `V26=26*{1,...,12} union {339}` is the transparent coherent obstruction.
  Its first rational witness is `2/27`, and `M(V26)=1/13`.
- `S*` has no prime dividing seven of its speeds and every leave-one-out gcd
  is one, so a twelve-speed common-factor dispatch does not explain it.  Its
  first rational witness is `3/26`, while its exact maximum is
  `M(S*)=101/470` at `t=167/470`.  Hence even a very loose,
  gcd-incoherent exact residual need not have a `q<=25` good period.

All claims in this paragraph are checked by exact integer or rational
arithmetic in the stored certificate artifact.

## 3. Why coherent blocks evade every `q<=25` claim

Let `C_c=c*{1,...,12}`.  If `a/q`, `q<=25`, is safe for `C_c`, reduce
`ac/q=b/Q`.  If `Q<=12`, the speed indexed by `j=Q` is at the origin.
If `15<=Q<=25`, one of the signed residues of `b^{-1}` has a representative
`j<=floor(Q/2)<=12`; then `||jb/Q||=1/Q<1/14`.  Therefore only
`Q in {13,14}` can survive.

If `C_c union {w}` is covering, a surviving denominator `q<=25` must itself
be `13` or `14`; the block then misses that divisor, so covering forces
`q|w`, which kills the proposed witness.  Thus every covering extension of
the dilated twelve-block has no rational witness of denominator at most 25.
This does not threaten LRC: THM-760 closes the one-coprime-exception sheet by
an adaptive lifted witness rather than a fixed raw denominator.

## 4. Exact scope of the historical banks

The 120 S312 rows establish only `120/120` for that random sample.  An exact
replay of the separate S105 generator gives `91/8260` rows with no
`q<=25` witness.  All 8,260 rows do have a witness by `q<=38`; the unique
row whose least denominator is `38` is

```
{1,2,...,11,338,420},  witnessed first at 3/38.
```

This is finite-exact for S105's stated interval cores, outlier pool, ordering,
and cap only.  It is not a uniform residual theorem; THM-566 already proves
that no global fixed raw denominator bound can hold on primitive covering
families.

## 5. Tournament-analysis guardrail

Runner tournaments and denominator tournaments are controlled forgettings of
the exact incidence object.  The proof object here is the blocker deck on
signed pairs, together with zero ownership; across several moduli it becomes
a witness-blocker hypergraph on `(q,a)` obligations.  Quotienting to one
vertex per modulus destroys which multiplier and which runner own a block,
and quotienting to runners destroys joint coverage of one denominator.

The certificate therefore reports tournament fingerprints only as telemetry:
vertices are moduli `15,...,25`; the pair-first and compression-first gauges
are both transitive, have score histogram `{0:1,...,10:1}`, no directed
3-cycles, singleton SCCs, and one Hamiltonian path, yet flip `47/55` edges.
That instability while the exact blocker verdict stays fixed is evidence that
the deck/hypergraph, not the tournament orientation, carries this theorem.

*Artifacts:* `04-computation/lrc14_q25_uniformity_refutation_codex_S3.py`
and its stored output; `04-computation/lrc14_s105_q25_bank_audit_codex_S3.py`
and its stored output.
