# Minority-anchor owner renewal and endpoint transitions

**Status:** PROVED elementary component lemma + FINITE-EXACT hostile probe;
not canonized here, and not a proof of LRC(14).

This note records an orthogonal pass on the degree-twelve minority anchor in
[THM-4330](../01-canon/theorems/THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve.md).
The stochastic-process analogy is useful only as vocabulary: the actual object
is a deterministic finite interval-reachability process with exact rational
renewal tolls.  No independence, stationarity, or probability claim is used.

## Inheritance pass and concept board

The closest proved mechanisms are:

- THM-4330, MC2a--MC2d: the exact two-sheet owner-bit obstruction;
- [THM-1156](../01-canon/theorems/THM-1156-tooth-seam-chi7-bipartition.md):
  the directed tooth-seam numerator and the need for backup at a strict
  endpoint;
- [THM-1244](../01-canon/theorems/THM-1244-slowest-spoke-component-handoff-debt.md):
  irredundant chronological tooth chains on a protected component;
- [THM-592](../01-canon/theorems/THM-592-radius-derivative-structure.md):
  pair-sum merge radii and endpoint-owner evolution; and
- [THM-2075](../01-canon/theorems/THM-2075-safe-child-homeomorphism-and-wall-word-conjugacy.md):
  componentwise constant sheet addresses and the warning that forgetting
  addresses destroys circular placement.

The canonical hostiles are THM-4330's two `h=420` rows with `P=1287` and
`P=9009`.  The corrected near miss is to treat the owner bits as independent
Bernoulli blockers, or to infer danger from a failed marginal capacity test.
An individual bit is gauge-dependent, and both hostiles are safe.  The
least-used relevant sidecar is the pair consisting of the physical component
address `(j,epsilon)` and the actual tooth address `n`, not merely the owner
speed.

The live concept board was:

```text
missing owner bit | spanning tooth | collision-backed renewal
pair-sum merge slack | binary component address | endpoint slack
```

## The deterministic shortest-cover-chain lemma

Let

```text
S_h={2h} union W,
W a finite set of positive odd integers.
```

At radius `1/14`, the `2h` physical anchor-safe components are

```text
I_k=[L_k,R_k]
   =[(14k+1)/(28h),(14k+13)/(28h)],     0<=k<2h,       (1)
```

each of width `3/(7h)`.  Write the open danger tooth of the odd speed `w`
with address `n` as

```text
T(w,n)=((14n-1)/(14w),(14n+1)/(14w)).                  (2)
```

Fix `k`.  Exactly one of the following holds.

1. **Missing bit.**  Some `t in I_k` belongs to no tooth `(2)`.  Then `t`
   is a `1/14`-lonely witness for `S_h`.
2. **One-state span.**  One tooth `T(w,n)` strictly contains the closed
   interval `I_k`.  Necessarily

   ```text
   3w<h.                                               (3)
   ```

3. **Collision-backed renewal.**  There is a farthest-reach greedy covering chain

   ```text
   T(w_1,n_1),...,T(w_m,n_m),          m>=2,           (4)
   ```

   of minimum cardinality among all tooth chains covering `I_k`.  Its left and
   right walls increase strictly and its consecutive teeth cross properly.
   Put

   ```text
   D_i=n_(i+1)w_i-n_iw_(i+1),
   q_i=w_i+w_(i+1)-14D_i.                              (5)
   ```

   Then

   ```text
   D_i>0,
   q_i in 2*gcd(w_i,w_(i+1))*Z_(>0),                  (6)

   |T(w_i,n_i) intersect T(w_(i+1),n_(i+1))|
      =q_i/(14w_iw_(i+1))
      >=1/[7 lcm(w_i,w_(i+1))].                       (7)
   ```

Thus a covered component either has a constant owner tooth with strict
endpoint slack, or pays a same-bit collision of a quantized positive size.
If neither occurs, the missing owner bit is already a witness.
The collision bound is sharp: the probe finds `h=10`, `W={3,13}`, `k=6`
with chain `((13,4),(3,1))` and `q_1=2`.

There is an exact renewal budget, not just the inequalities above.  Define

```text
A=28h w_1(L_k-left(T(w_1,n_1)))
 =w_1(14k+1)-2h(14n_1-1)>0,

B=28h w_m(right(T(w_m,n_m))-R_k)
 =2h(14n_m+1)-w_m(14k+13)>0.                          (8)
```

Then `A,B` are positive odd integers and

```text
3/(7h)
 =sum_(i=1)^m 1/(7w_i)
  -A/(28h w_1)-B/(28h w_m)
  -sum_(i=1)^(m-1) q_i/(14w_iw_(i+1)).                (9)
```

In particular, every renewal path pays both endpoint tolls and every
collision toll.  Repeated speed labels are allowed in `(4)`; the states are
addressed teeth, not runners.

### Proof

If the open teeth do not cover `I_k`, compact one-dimensional geometry gives
case 1.  Otherwise choose, among all teeth containing `L_k`, one with the
farthest right wall.  If that wall does not pass `R_k`, choose among all teeth
strictly containing the current right wall one with the farthest right wall,
and iterate.  Coverage guarantees a choice; finiteness of the teeth meeting
`I_k` guarantees termination.

The greedy rule excludes containment.  Consecutive left and right walls
therefore increase strictly and cross properly.  It also implies that a
tooth two steps later cannot overlap the tooth two steps earlier: otherwise
it was already eligible and extended farther at the previous choice.  Hence
only consecutive chain teeth overlap.  The chain is shortest by the standard
frontier exchange: after `r` teeth, no `r`-tooth chain beginning at `L_k` can
reach farther than the greedy frontier.  Indeed, an interval extending past
that frontier would contain it and would have been eligible for the greedy
choice; intervals ending before it do not improve the frontier.

If the chain has one tooth, its width `1/(7w)` is strictly greater than the
width `3/(7h)` of `I_k`, proving `(3)`.  For a proper crossing, both tooth
centres increase, so

```text
D_i=n_(i+1)w_i-n_iw_(i+1)>0.
```

Direct subtraction, always in the order "exiting right wall minus entering
left wall", gives

```text
right(T(w_i,n_i))-left(T(w_(i+1),n_(i+1)))
 = [w_i+w_(i+1)-14D_i]/(14w_iw_(i+1))
 = q_i/(14w_iw_(i+1))>0.
```

Thus the sign in `(5)` is `q_i=w_i+w_(i+1)-14D_i`, not its negative.
Writing `g_i=gcd(w_i,w_(i+1))`, both `D_i` and `w_i+w_(i+1)` are divisible
by `g_i`; because `g_i`, `w_i/g_i`, and `w_(i+1)/g_i` are odd, `q_i/g_i`
is positive and even.  Hence `q_i in 2g_i Z_(>0)`, which proves `(6)` and
the lower bound in `(7)`.

At the left endpoint, `(8)` is odd minus even; at the right endpoint it is
even minus odd.  The strict containments of `L_k` and `R_k` give the displayed
positive signs.  Therefore `A,B in 2Z_(>=0)+1`, with in particular
`A>=1` and `B>=1`.  Finally, the union of the chain is one interval, only
adjacent overlaps have positive length, and removing its two positive
endpoint overhangs leaves `I_k`.  Inclusion--exclusion therefore gives `(9)`
with all three toll families subtracted.  QED.

## MC2 owner bits and the address word

Write `k=j+h epsilon`, where `0<=j<h` and `epsilon in {0,1}`.  Then `(1)` is
the physical lift

```text
I_k=(J_j+epsilon)/2
```

of the `j`th component `J_j` of THM-4330's quotient set `G_h`.  If
`T(w,n)` is active on this physical component, the quotient nearest integer
is

```text
N_w(x)=2n-epsilon*w,
N_w(x) mod 2=epsilon.                                  (10)
```

Thus every state in one chain `(4)` has the same MC2 owner bit.  Case 1 is
literally a missing owner bit in MC2d.  In case 3, at the endpoint of an
exiting tooth that owner is at equality and hence safe; the next tooth is
strictly active there.  The transition is therefore backed by another owner
of the same bit, exactly the strict-event repair mechanism of THM-771 and
THM-1156.

This is also the appropriate extension of THM-2075's address word.  The
binary sheet address `epsilon` stays constant on a component.  Unique-child
transport would give a one-state word; the minority-anchor obstruction can
instead renew through addressed owner teeth.  Quotienting `(w,n)` to `w`, or
quotienting `(j,epsilon)` to `j`, loses the transition placement.

The THM-592 connection is exact but limited.  An edge `(5)` has pair-sum
merge radius

```text
rho_i=D_i/(w_i+w_(i+1))<1/14,                          (11)
```

and

```text
q_i=14(w_i+w_(i+1))(1/14-rho_i)>0.                    (11a)
```

Thus `q_i` is precisely its integer distance from the `1/14` wall.  These
interior overlaps can be invisible to the global lonely-measure derivative;
they must not be mislabelled as exposed THM-592 kinks.  The renewal carrier
retains them because they preserve the pointwise cover predicate.

## Exact occupation identity

There is also a deterministic continuum renewal identity.  Let `E` be the
anchor-safe set, so `mu(E)=6/7`, and for `t in E` let `m(t)` be the number of
strict tail-danger teeth containing `t`.  Define

```text
F     =mu({t in E:m(t)=0}),
Omega =integral_E (m(t)-1)_+ dt,
A_0   =sum_(w in W) mu(D_(2h) intersect D_w).
```

Since each tail danger set has mass `1/7`, its load on `E` is
`|W|/7-A_0`.  Pointwise, `m=1_(m>=1)+(m-1)_+`; integration therefore gives

```text
F=(6-|W|)/7+A_0+Omega.                                (12)
```

For the degree-twelve tail used in the minority branch this becomes

```text
F=A_0+Omega-6/7.                                      (13)
```

This is the continuum analogue of THM-771's defect bookkeeping.  It is exact,
not a heuristic expectation: a hypothetical strict counterexample has `F=0`
and hence must satisfy `Omega=6/7-A_0`.  It does not itself force `F>0`, because
the anchor-tail overlap and tail-collision masses can balance exactly.

## Global collision debt

Let `P_h(W)` be the set of component addresses `k` for which some single
tail tooth strictly contains `I_k`.  If `S_h` were a strict counterexample,
case 1 would be impossible at every address.  Therefore at least

```text
2h-|P_h(W)|                                             (14)
```

component addresses would have to carry a collision-backed renewal.  This
is a deterministic count of proof obligations, not a claim that collision
events are independent or disjoint in physical time.

## Exact hostile probe

The companion artifacts are:

```text
04-computation/lrc14_minority_owner_renewal_probe_owner_walk_20260902.py
05-knowledge/results/lrc14_minority_owner_renewal_probe_owner_walk_20260902.out
```

They use `Fraction` arithmetic.  A cell-by-cell wall/midpoint implementation
independently agrees with greedy reachability on `46,124` small component
instances.  Normal, optimized, and hash-seeded runs reproduce the frozen
output.

The status split is essential: the shortest-cover-chain lemma and occupation
identity are **PROVED for every positive integer `h` and every finite odd
tail `W`**.  The numerical profiles below are **FINITE-EXACT controls for only
the two displayed `h=420` rows**; no census or extrapolation is asserted.

For each THM-4330 joint hostile (`h=420`, `P=1287` or `P=9009`) the exact
discrete `840`-component classification is identical:

```text
missing owner bit       726
one-state span          110   (all owned by w=11)
collision renewal         4
maximum chain length       3
minimum q_i                4
```

Their exact occupation masses are not identical and are recorded separately
in the frozen output; only the displayed component counts coincide.

Thus `|P_h(W)|=110`, and the nonspanning obligation count `(14)` is `730`.
Only four are paid by renewals; the other `726` expose a missing bit.  Both
rows share the first exact endpoint witness

```text
t=939/141274,
min_(v in S_h)||vt||=1/14,
unique binding speed=10091.                            (15)
```

This is an exact finite certificate for the already-known safety of those two
rows: failed component reachability returns the witness `(15)`.  It is a new
certificate *format*, not a new analytic sufficient condition or an LRC family
theorem.  In general, the lemma is a sharper obstruction/address compiler;
without a hypothesis forcing reachability failure at some `k`, it does not
prove that a missing bit exists.  The `P` coordinate is invisible in these two
controls because the failure occurs before either `1287` or `9009` is needed.

Replay:

```bash
python3 -B 04-computation/lrc14_minority_owner_renewal_probe_owner_walk_20260902.py \
  | diff -u 05-knowledge/results/lrc14_minority_owner_renewal_probe_owner_walk_20260902.out -
python3 -B -O 04-computation/lrc14_minority_owner_renewal_probe_owner_walk_20260902.py \
  | diff -u 05-knowledge/results/lrc14_minority_owner_renewal_probe_owner_walk_20260902.out -
PYTHONHASHSEED=937 python3 -B \
  04-computation/lrc14_minority_owner_renewal_probe_owner_walk_20260902.py \
  | diff -u 05-knowledge/results/lrc14_minority_owner_renewal_probe_owner_walk_20260902.out -
```

Raw-LF SHA-256 after this note was written:

```text
script  18f42dfc2ed67706cc61e84c1a03fab937d357a41afc7d8fa0c419402e3fd083
output  57821c6521a081943a5be377197f819b69570fdb611caacc007d12f45621240e
```

## Connection contract

```text
source:       labelled minority-anchor row {2h} union W, with W odd
target:       a finite directed interval graph on addressed tail teeth over
              every physical component address (j,epsilon)
map:          lift each G_h component to I_(j,epsilon), intersect it with all
              strict tail teeth, and orient proper crossings left-to-right
preserved:    exact 1/14 witness/nonwitness predicate; strict equality;
              MC2 owner bit; component, speed, and tooth addresses;
              endpoint and pairwise overlap slack
destroyed:    alternative cover chains after choosing the greedy path;
              geometry away from the selected anchor component; correlations
              between different component resets; exposed/nonexposed status
              in the global radius derivative
sidecar:      28h endpoint denominator, sheet bit epsilon, tooth integer n,
              q_i, gcd/lcm, and reset order around the 2h components
decisive test: directed reachability from the left anchor wall to beyond the
               right anchor wall; one failed address returns an exact witness
```

## Stopping reason and next use

The stochastic analogy has now paid out its deterministic content.  Treating
renewals probabilistically would immediately discard the address and gcd
coordinates that make `(6)--(15)` exact.  The remaining obstacle is reuse:
one high-frequency tail contributes many teeth, and the same owner pair can
pay transition obligations at many component addresses.  Therefore the local
collision debt `(14)` does not yet aggregate to a contradiction.

The next theorem-sized target is an exact reuse bound for transition masks:
for a fixed odd pair `(u,v)`, count the anchor addresses `k` on which a proper
`u`--`v` crossing can occur, retaining its `q` residue.  Comparing the union of
those pair masks with `P_h(W)` would turn the renewal lemma into a finite
component-address certificate on the `420|h` wall.  Until that reuse bound is
proved, the present result is a lossless local compiler and a sharp hostile
diagnostic, not a closure of the minority branch.
