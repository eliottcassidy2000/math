---
id: HYP-6820
title: Uniformity audit for the LRC(14) q<=25 good-period claim and the n=12 sporadic branch
status: PARTIALLY RESOLVED — uniform q<=25 is DISPROVED; n=12 is uniformly finite, binding-scale/sheet-stratified, and shallow-exact through height twelve, but branch emptiness remains OPEN
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
  - THM-768
  - THM-769
  - THM-770
  - HYP-6750
  - HYP-6775
related:
  - HYP-6780
  - HYP-6785
  - HYP-6800
  - HYP-6815
  - MISTAKE-143
---

# HYP-6820 — uniformity audit (formerly HYP-6810; renumbered after collision)

This audit began with two concrete proof obligations requested together:

1. determine whether every stated residual family has a rational lonely witness
   of denominator `q <= 25`, and separate exact banks from uniform claims;
2. prove, or reduce without quantifier loss, that the super-lonely-core (sporadic)
   branch of primitive tight 12-speed families is empty.

They now have different outcomes.  The first statement is false.  The second
has acquired several uniform reductions, but the requested emptiness theorem is
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

The exact max-peel tooth atlas gives a further guardrail.  It exhausts all
primitive twelve-subsets `A` of `{1,...,20}` with `w=max(A)`,
`M(A\{w})>1/12`, and `M(A)<=1/10`, plus eleven tight AP deletion controls.
The component criterion agrees with an independent exact-maximum calculation
on all `2,464` rows.  Every one of the `2,453` escaping rows nevertheless has
cyclic nearest-tooth winding one, and `1,972` have pure endpoint owners.  The
first liar to “pure endpoints plus winding one implies cover” is
`{1,...,10,12,13}`; even requiring a transitive component-phase tournament
with one Hamiltonian path fails for `{1,...,9,11,12,15}`.  Their minimum
slacks are respectively `-11/26` and `-43/91`.  Thus winding is a checksum,
not the missing separator: the width term in
`sigma_J(w)=1/13-||wc_J||-w h_J` must survive every quotient.

## C. Binding-scale recursion: what the residue picture was missing

THM-768 eliminates one tempting deep configuration: if the largest speed is
the unique multiple of 13, an explicit perturbation of a shallow prime-grid
point makes all twelve runners strictly safer than `1/13`.  Thus a tight set's
largest speed cannot be its unique 13-multiple.

THM-769 gives the scale-normal statement at **every** rational global maximum.
Write a reduced maximizer as `p/(13s)`.  Its multiplied residues lie in the
packet `[s,12s]`, both endpoints occur, and endpoint owners are divisible by
`s`.  With

```text
E={v:s|v}=sU,                  F=A\E,
```

the familiar complete nonzero residue system is exactly the shallow `s=1`
branch.  If a 13-multiple blocks all shallow points, every maximizer is deep
and `F` must cover all `s` lifts of every point in
`G_U={tau:phi_U(tau)>1/13}`.  Putting

```text
D_w=s/gcd(w,s),
```

gives the necessary capacity inequality

```text
sum_(w in F) (floor(2D_w/13)+1)/D_w >= 1.
```

In particular `|F|>=2`.  Equality at two exceptions forces `s=2`, ten even
speeds and two odd speeds; tightness is then **equivalent** to persistent
opposite nearest-integer parity over all of `G_U`.  At three exceptions,
either a half-sheet tightener occurs or the equality edge is `s=3`, with nine
multiples of three and three nonmultiples persistently owning all three sheet
colours.  For `r=|F|<=6`, some exception satisfies
`D_w<=13r/(13-2r)`.  These are uniform reductions, not empirical patterns.

THM-770 settles a very large but still bounded part of the shallow branch.
For the labelled packets

```text
W(k)={r+13k_r:1<=r<=12},       0<=k_r<=12,
```

an exact unique-owner CSP represents all `13^12=23,298,085,122,481` rows
without literal enumeration.  It has exactly thirteen zero-defect leaves,
`c*{1,...,12}` for `c=1,...,12,14`, and only the `c=1` leaf is primitive.
Consequently a primitive shallow tight set with `max A<=168` is the AP.
Dilation gives the exact defect law `chi_13(cW)=c chi_13(W)`, so gcd descent is
lossless.  This result does not supply an unbounded height descent.

## D. Remaining proof obligation

The uniform theorem now has two explicit residuals:

1. **Shallow descent.**  Prove that a primitive full-nonzero-residue packet
   with `chi_13=0` descends into THM-770's height-twelve box, or extend its
   owner-CSP by a scale-free coherence argument.
2. **Deep colour cover.**  Rule out THM-769's persistent folded parity cover
   at `s=2` and its higher-sheet analogues (or prove they descend to a smaller
   primitive tight packet).

Equivalently in the original top-peel language, prove that every primitive
nonextremal eleven-speed core below THM-763's finite height has at least one
safe component whose rational endpoint band is incompatible with every top
speed in its THM-759/766 cone.

For the shallow branch, let `kappa` be the number of components of the open
danger-tooth graph and `P_splice` the number of protected end/start
coincidences.  The settled twelve-speed LRC ensures that the open danger union
is not the whole circle, and exact endpoint sweep gives

```text
chi_13 = kappa-P_splice = number of open 1/13-safe components.
```

For full nonzero residues, the nominal complementary-pair capacity is
`P*=2 sum_(r=1)^6 gcd(w_r,w_(13-r))`, so

```text
chi_13=(kappa-P*)+(P*-P_splice).
```

This separates overlap-rank shortage from third-runner blocker debt.  The
height-one lift cube had 4,085 rank-short rows, 9 blocker-debt rows, and one
zero-defect row, the nonprimitive doubled AP; all 4,094 primitive rows had
`chi_13>=2`.  THM-770 supersedes that cube by an exact height-twelve theorem,
but the same warning remains: bounded exactness is not the global coherence
lemma.

## E. Information-preservation / Tournament Analysis

The exact vertices are witness obligations `(q,a)`, safe components, endpoint
splices, or the `s` sheet fibres—not runners by default.  Modulus and runner
tournaments are telemetry: changing gauges flips many edges while the blocker
verdict stays fixed, and pairwise component compatibility does not imply
containment.  The component-phase atlas makes the same loss explicit: a
transitive tournament with singleton SCCs and one Hamiltonian path can still
have negative slack at every component because phase order forgets interval
width.  In THM-770 all 66 residue-obligation comparisons tie, while the
endpoint-owner hypergraph still distinguishes thirteen solutions.  The
deciding objects are therefore:

- the zero-owner/signed-pair blocker deck for small periods;
- the component-tooth incidence hypergraph with endpoint widths, divisor pins,
  and core-maximizer residues for tight completions;
- the off-sheet-runner by sheet incidence cover, with effective orders
  `s/gcd(w,s)` and persistent colour ownership over the quotient loose set.

These objects preserve the LRC predicate.  Their tournament quotients destroy
joint blocker ownership, multiplier identity, scale, ramification, and
simultaneous alignment.
