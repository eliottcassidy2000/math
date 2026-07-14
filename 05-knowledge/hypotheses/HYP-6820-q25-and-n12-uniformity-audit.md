---
id: HYP-6820
title: Uniformity audit for the LRC(14) q<=25 good-period claim and the n=12 sporadic branch
status: PARTIALLY RESOLVED — uniform q<=25 is DISPROVED; n=12 is uniformly finite and tooth-stratified but branch emptiness remains OPEN
source: codex-2026-07-14-S3
renumber_note: reserved as HYP-6810 by codex-S3, which collided with opus-S298's earlier-pushed
  HYP-6810 claim (the assembly write-up); renumbered to HYP-6820 by opus-2026-07-14-S299 per the
  first-pusher protocol, as codex-S2's atlas itself requested; subsequently updated by codex-S3.
depends_on:
  - THM-758
  - THM-759
  - THM-762
  - THM-763
  - THM-764
  - THM-765
  - THM-766
  - HYP-6750
  - HYP-6775
related:
  - HYP-6780
  - MISTAKE-143
---

# HYP-6820 — uniformity audit (formerly HYP-6810; renumbered after collision)

This audit began with two concrete proof obligations requested together:

1. determine whether every stated residual family has a rational lonely witness
   of denominator `q <= 25`, and separate exact banks from uniform claims;
2. prove, or reduce without quantifier loss, that the super-lonely-core (sporadic)
   branch of primitive tight 12-speed families is empty.

They now have different outcomes.  The first statement is false.  The second
has acquired two uniform reductions, but the requested emptiness theorem is
still open.

## A. The `q<=25` claim is disproved

THM-762 (independently restated in THM-764) proves the exact criterion visible
in the range `15<=q<=28`.  For a
covering family `S`, a `q`-witness exists exactly when:

1. no speed is divisible by `q`; and
2. the signed unit-pair deck
   `B_q(S)={[s] modulo +/-1 : gcd(s,q)=1}` misses a class of
   `(Z/qZ)^x/{+1,-1}`.

The proof is elementary: a nonunit multiplier reduces to a covered denominator
at most 14; for a unit multiplier and `q<=28`, the only strictly unsafe integer
residue distances are `0` and `1`.

Two primitive covering exact residuals refute uniform `q<=25`:

```
26*{1,...,12} union {339}                         first witness 2/27,
{81,91,131,151,157,196,258,274,313,328,330,339,348} first witness 3/26.
```

The second has no prime common to seven speeds, every leave-one-out gcd is one,
and exact `M=101/470`; coherence or near-tightness therefore cannot rescue the
claim.  Both pass the rational analogue of S312's band-residual predicate.
MISTAKE-143 corrects the universal wording.

The historical banks now have exact scope:

- S312: `120/120` sampled rows with `q<=25`, not an exhaustive class;
- S105: exact replay of its capped interval-core generator gives `91/8260`
  rows with no `q<=25` witness, but every row witnesses by `q<=38`.  The unique
  least-denominator-38 row is `{1,...,11,338,420}`, first witnessed by `3/38`.

THM-566 remains the global guardrail: no fixed raw denominator ladder can cover
all primitive covering families.

## B. The `n=12` branch is finite and tooth-stratified, not proved empty

Let `A={a_1<...<a_12}` be primitive and tight, `M(A)=1/13`, and let
`P=A\{a_12}`.  The sporadic branch is `M(P)>1/12`.

THM-763 gives the first quantifier-honest global finite bound:

```
sum A <= 78^11 = 650190514836423555072.
```

Thus uniform emptiness is a finite decision problem, conditional only on the
settled lower-dimensional LRC input, but a naive enumeration is infeasible.

THM-766 supplies a scale-normal reduction.  For a tight `n`-set, with
`m=a_1`, `b=a_(n-1)`, and `w=a_n`,

```
w >= n m,
b/m >= n^2/(n+2).
```

If `b<nm`, the first core-safe interval must fit one `w`-danger tooth, so for
some `1<=k<=n-1`,

```
((n+1)k-1)m <= w,
n w <= ((n+1)k+1)b.
```

At `n=12`, every candidate therefore has `a_11/a_1>=72/7`; the sub-12 span
lies in only eleven explicit projective cones.  At an exact core maximum
`t_0`, the stronger band

```
||a_12 t_0|| <= 1/13 - a_12(M(P)-1/13)/a_11
```

retains the residue alignment lost by THM-759's scalar ratio cap.

The faithful completion predicate is

```
G_P={t:min_(p in P)||pt||>1/13} subseteq
D_w={t:||wt||<=1/13}.
```

Every connected component of `G_P` must fit one `w`-tooth simultaneously.
THM-765 strengthens this component language with a gcd-deck theorem: every
leave-one-out core of a primitive tight twelve-set is itself primitive.  It
therefore removes all imprimitive-core sporadic candidates at every height,
not just in a box.
Total measure alone cannot decide this: for the nonextremal core
`P={1,...,10,12}`, exact `|G_P|=461/8190<2/13=|D_w|`.

An exact regression on the historical 77 nonextremal cores in `{1,...,13}`
uses the true THM-759 cap rather than an arbitrary killer window.  It leaves
790 primitive completions: 40 narrow-core and 750 wide-core.  THM-766's cone
test eliminates all 40 narrow candidates; exact pair-sum/difference/half-turn
evaluation finds zero tight completions among all 790, with bank minimum
`M=1/12`.  This is finite-exact for that box only.

## C. Remaining proof obligation

The requested uniform emptiness theorem is still missing.  The residual is:

> Prove that every primitive nonextremal eleven-speed core below THM-763's
> finite height has at least one safe component whose rational endpoint band
> is incompatible with every top speed in its THM-759/766 cone.

Equivalently, classify the rare zero-defect endpoint-splice packets.  A useful
next lemma is splice-lattice coherence: when the twelve residues occupy every
nonzero class modulo 13 and every safe gap is closed by complementary-pair
splices, the six splice lattices should share one scale, forcing a dilate of
`{1,...,12}` and hence the primitive AP.

Concretely, let `kappa` be the number of components of the open danger-tooth
graph and `P_splice` the number of end/start coincidences with no third tooth
active.  Provided the open danger union is not the whole circle, exact
endpoint sweep gives

```
chi_13 = kappa-P_splice = number of open 1/13-safe components.
```

That proviso is automatic here: the settled twelve-speed LRC gives a point
outside the open danger union at threshold `1/13`.  Without it, the sweep's
component-start count is zero whereas the full circle has one topological
component.

For full nonzero residues, the nominal complementary-pair capacity is
`P*=2 sum_(r=1)^6 gcd(w_r,w_(13-r))`, so

```
chi_13=(kappa-P*)+(P*-P_splice),
```

separating overlap-rank shortage from third-runner blocker debt.  The exact
height-one lift cube has 4,085 rank-short rows, 9 blocker-debt rows, and one
zero-defect row, the nonprimitive doubled AP.  All 4,094 primitive rows have
`chi_13>=2`.  This is finite evidence, not the coherence lemma.

## D. Information-preservation / Tournament Analysis

The exact vertices are witness obligations `(q,a)`, safe components, or
endpoint splices—not runners by default.  Modulus and runner tournaments are
telemetry: changing gauges flips many edges while the blocker verdict stays
fixed, and pairwise component compatibility does not imply simultaneous tooth
containment.  The deciding objects are therefore:

- the zero-owner/signed-pair blocker deck for small periods;
- the component-tooth incidence hypergraph with endpoint widths, divisor pins,
  and core-maximizer residues for tight completions.

These objects preserve the LRC predicate.  Their tournament quotients destroy
joint blocker ownership, multiplier identity, scale, and simultaneous
alignment.
