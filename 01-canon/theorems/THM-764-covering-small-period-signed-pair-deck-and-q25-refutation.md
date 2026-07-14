---
id: THM-764
title: The covering small-period signed-pair deck and the failure of the q<=25 bound
status: PROVED (elementary residue criterion; exact residual certificates are computer-checked over integers and rational intervals)
source: codex-2026-07-14-S3
depends_on:
  - THM-366   # covering/non-covering divisor gate
  - THM-668   # exact rational maximizer audit used by the companion script
  - THM-755   # definition of the uncapped residual tested below
related:
  - THM-762   # first-pushed statement of the same criterion; includes S105 replay
  - THM-758
  - HYP-6750
  - HYP-6780
  - HYP-6820
---

# THM-764 — The signed-pair deck and the failure of the `q<=25` bound

This is an independently developed expanded restatement of THM-762, using the
same exact counterexample artifact.  It is retained because its band-residual
and gcd-incoherence certificates are presented in greater detail.

## Statement

Let `S` be a finite set of positive integer speeds that is **covering** in the
LRC(14) sense: for every `d in {2,...,14}`, some `s in S` is divisible by
`d`.  For `15<=q<=28`, define

```text
U_q = (Z/qZ)^*/{+1,-1},
B_q(S) = { {+s,-s} : s in S and gcd(s,q)=1 } subset U_q,
Z_q(S) = {s in S : q divides s}.
```

Then `S` has a closed LRC(14) witness `a/q`, `0<a<q`, if and only if

```text
Z_q(S) is empty and B_q(S) is not all of U_q.            (1)
```

More precisely, every nonunit multiplier `a` is blocked by a zero phase, and,
for a unit multiplier, the only unsafe integer residue classes are
`0,+1,-1`.  Thus a missing signed unit pair `b in U_q\B_q(S)` gives the
witness `a=b^(-1) mod q`, up to sign.

The proposed uniform conclusion that every covering thirteen-speed residual
has such a witness with `q<=25` is false, even inside the exact S312
band-residual predicate.  Two counterexamples are

```text
V26 = (26,52,78,104,130,156,182,208,234,260,286,312,339),

S*  = (81,91,131,151,157,196,258,274,313,328,330,339,348).
```

Both are primitive, covering, have thirteen speeds above `14`, diameter at
most `339`, pass every leave-one-out S312 band-residual inequality, and
have no witness at any denominator `15<=q<=25`.  In fact:

```text
M(V26) = 1/13,       first rational witness = 2/27;
M(S*)  = 101/470,    first rational witness = 3/26.
```

The second example is deliberately gcd-incoherent: no prime divides seven of
its speeds, and every leave-one-out gcd is `1`.  Hence neither looseness nor
the absence of a large common-factor pack repairs the `q<=25` assertion.

## Proof of the deck criterion

Fix `a` with `0<a<q`.  If `g=gcd(a,q)>1`, reduce `a/q` to denominator
`q_0=q/g`.  Since `q<=28`, one has `2<=q_0<=14`.  Covering supplies a speed
divisible by `q_0`, and that speed has phase zero at `a/q`.  Therefore a
nonunit multiplier is never a witness.

Now suppose `gcd(a,q)=1`.  The phase of a speed `s` is safe exactly when

```text
14 min(sa mod q, -sa mod q) >= q.                       (2)
```

For `15<=q<=28`, one has `1<q/14<=2`.  Consequently, the integer distances
that violate (2) are exactly `0` and `1`; equivalently the unsafe residue
classes are `0,+1,-1`.  A zero residue occurs for a unit `a` exactly when
`q|s`.  If there is no zero owner, a residue `+1` or `-1` occurs exactly when
the signed pair of `a^(-1)` lies in `B_q(S)`.  Therefore some unit `a` is safe
for every speed exactly when there is no zero owner and some signed unit pair
is missing.  This proves (1).  ∎

The upper endpoint `28` is exact for this simple formulation.  At `q=29`,
integer residue distance `2` is also strictly below `q/14`, so the blocker
alphabet must be enlarged.

## Transparent counterexample

For any `15<=q<=25` and any multiplier `a`, put `g=gcd(26a,q)`.

- If `g>1`, then `i=q/g` is an integer in `{1,...,12}` and
  `q | 26ia`; the speed `26i` is a zero blocker.
- If `g=1`, choose the least signed representative `i` of
  `(26a)^(-1) mod q`.  Then `1<=i<=(q-1)/2<=12` and
  `26ia=+1` or `-1 mod q`; the speed `26i` is a sign-pair blocker.

Thus the first twelve speeds of `V26` alone block every multiplier at every
target denominator.  The family is primitive because `gcd(26,339)=1`, and
its displayed covering owners are checked directly.  Scaling the tight AP
core gives the upper bound `M(V26)<=1/13`.  At

```text
t = 27/338
```

the speed `26i` has the same phase as `i/13`, while the last runner has
distance `27/338>1/13`; hence `M(V26)=1/13`.  At `2/27` every speed has
clearance at least `2/27`, so `q=27` is the first rational escape.

## Exact residual and incoherence certificates

The companion program uses `Fraction` arithmetic and the repository's exact
safe-interval engine.  For each deletion `P=S\{w}` it computes the strict
good-set measure `|G'_P|` and component count `r_P`, then proves

```text
pi*w*|G'_P| < (22/7)*w*|G'_P| < r_P.                   (3)
```

This verifies the S312 band-residual predicate without floating point.
For `V26`, the largest exact value of `w|G'_P|` is

```text
51768233/1993320,
```

and the uniform comparison in (3) is below `290<=r_P`.  For `S*`, it is

```text
106329879949471805834083/2072849518344678893820,
```

and the corresponding comparison is below `468<=r_P`.

The program prints every zero-owner and signed-pair deck for `q=15,...,28`,
checks (1) against every multiplier, and evaluates `M` with the exact
pair-sum ruler engine.  For `S*`, the largest prime packets have sizes only
`6` (prime `2`), `5` (prime `3`), and `2` (prime `7`); all thirteen
leave-one-out gcds equal `1`.  Nevertheless `B_q(S*)=U_q` or `Z_q(S*)` is
nonempty for every `15<=q<=25`.  This completes both counterexample
certificates.  ∎

Companion artifacts:

```text
04-computation/lrc14_q25_uniformity_refutation_codex_S3.py
05-knowledge/results/lrc14_q25_uniformity_refutation_codex_S3.out
```

## What replaces the failed bound

The exact obstruction is a **blocker deck**, not a fixed denominator ceiling.
For `15<=q<=28`, each modulus carries a zero-owner bit and a subset of signed
unit pairs.  A proof by rational periods must either show that some adaptive
modulus has an empty pair, or derive structure from simultaneous completeness
of these decks.  The counterexamples show that raw looseness, bounded
diameter, the band-residual test, and the absence of a large gcd
packet do not force that empty pair by `q=25`.

## Tournament analysis and challenged vertex choice

Using runners as vertices loses the simultaneous covering predicate in (1).
For telemetry, the companion program instead uses the moduli `15,...,25` as
vertices.  The pairwise observable compares signed-pair retention and shortest
blocker-certificate cost; switching from the pair-first to the
compression-first gauge flips `47` of `55` edges for both examples.  Ties use
the fixed Hamiltonian path

```text
15 -> 16 -> ... -> 25.
```

Both gauges happen to be transitive: score histogram `{0:1,...,10:1}`, no
directed triangles, singleton SCCs, and one Hamiltonian path.  This dramatic
edge instability with no change in the proof verdict is the guardrail: the
tournament is a useful fingerprint, but the zero/sign-pair deck is the exact
sidecar that preserves the LRC predicate.  Alternate vertices considered
were runners, gaps, fixed sections, section boundaries, wall events,
residues, cover arcs, Fourier modes, matroid circuits, and proof obligations.
