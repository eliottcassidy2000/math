---
id: THM-4330
title: "LRC(14) affine two-adic root types and anchored fixed-pool entry sieve"
status: >
  PROVED RELATIVE TO CITED LRC THROUGH THIRTEEN TOTAL RUNNERS,
  THM-2060/4150/4156/4191/4326 + VERIFIED-EXACT.  The graph of pairwise
  differences of minimum two-adic valuation is the complete bipartite graph
  between the two primitive parity classes.  Its degree at a prescribed
  runner is exactly that runner's number of odd primitive relatives.  The
  singleton root type is safe for every prescribed runner, so the 12+2 split
  is the first unresolved dyadic type.  In its degree-two lane, midpoint
  collisions reduce to settled lower dimension and every collision-free row
  is literally an eleven-even/two-odd seam.  Re-referencing changes the
  distinguished runner and is not an entry proof.  A positive rational
  projective refactorization into at most two labels outside the fixed pool
  closes the seam; hence a counterexample in this branch has at least three
  outsiders after every such refactorization.  The actual primitive tail
  ratio gives the sharper sufficient test mu(G_H)>=mu(C_(p,q)).  LRC(14)
  and arbitrary-row entry remain open.
source: root + parity_entry + entry_corpus_audit / LRC14 continuation session, 2026-09-01
depends_on:
  - LRCUpTo13
  - "T. Sungkawichai and T. Trakulthongchai, Eleven, twelve, and thirteen lonely runners, arXiv:2604.23906v1 (preprint)"
  - THM-2060-crt-tail-coset-saturation
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4326-lrc14-rank-two-wall-graph-complete-typed-universe-closure
related:
  - THM-639-hamiltonian-path-frame-for-runner-families
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2888-eight-body-first-apex-global-pair-cap-atlas
  - THM-4203-fixed-pool-seventeen-body-depth-eight-haar-completion
script: 04-computation/lrc14_entry_parity_affine_classification_probe.py
output: 05-knowledge/results/lrc14_entry_parity_affine_classification_probe.out
script_sha256: 3e791e1a1a92028d8edd5bb70e90b3af248d97dab2a78991ff02dc7f4d4cdbfe
output_sha256: 6168aeccbb8c7bc5427aac49070948d61f5c1ceee4f268e785452e37a44a7226
audit: >
  PASS.  The graph theorem is proved symbolically below.  The exact probe
  checks all seven root-size profiles, all 3,060 fourteen-subsets of
  {-2,...,15}, all 858 configurations in a structured 12+2 universe, the
  AP and one-reference hostiles, literal and projective pool-entry controls,
  and 90 affine normalization images.  All 935,790 assertions pass.  Normal,
  optimized, and hash-seeded invocations produce the frozen transcript.
---

# THM-4330 -- affine two-adic root types and anchored pool entry

**PROVED RELATIVE TO CITED LRC THROUGH THIRTEEN TOTAL RUNNERS,
THM-2060/4150/4156/4191/4326 + VERIFIED-EXACT.  LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance pass

Let

```text
V={v_0,...,v_13} subset Z
```

be fourteen distinct labelled velocities.  Put

```text
g=gcd_(i<j)|v_i-v_j|,                 s=nu_2(g),          (1)
```

and form a labelled graph `Gamma(V)` on the runners by

```text
ij in E(Gamma(V))  iff  nu_2(v_i-v_j)=s.                 (2)
```

For a prescribed distinguished runner `r`, call

```text
u_i=(v_i-v_r)/g,                      i!=r,               (3)
```

its signed primitive relatives.  The values `|u_i|` form the primitive,
scale-equivalent positive-speed row after repetitions are identified.  If
`A={|u_i|:i!=r}`, the original relative row is `gA` and

```text
G_(gA)=m_g^(-1)(G_A),                    m_g(x)=gx.       (3a)
```

Thus the original and normalized rows have equivalent nonemptiness; the
prescribed anchor is unchanged.

> **Affine two-adic root theorem.**  The graph `Gamma(V)` is exactly the
> complete bipartite graph between the two parity classes of the normalized
> velocities.  Thus
>
> ```text
> Gamma(V)=K_(m,14-m)                                    (4)
> ```
>
> for some `1<=m<=13`, and
>
> ```text
> #{i!=r:u_i odd}=deg_Gamma(r).                          (5)
> ```
>
> Writing `k=min(m,14-m)`, the seven unordered root types and their edge
> counts are
>
> ```text
> k                 1   2   3   4   5   6   7
> |E(Gamma)|       13  24  33  40  45  48  49.           (6)
> ```

The root type `k=1` is safe at the LRC(14) threshold for **every prescribed
runner**.  Consequently a counterexample configuration must have `k>=2`.
This packages an old one-tail mechanism into the first intrinsic affine
classification; it is not a new proof of the one-tail lemma.

At `k=2`, the twelve majority-class runners have graph degree two and the two
minority-class runners have degree twelve.  For a prescribed degree-two
runner exactly one of the following holds:

1. two signed relatives have the same absolute value, equivalently the
   prescribed runner is the midpoint of two other velocities; settled
   lower-dimensional LRC closes the row;
2. the thirteen absolute relatives are distinct and have the literal form

   ```text
   2H union {a,b},       |H|=11,       a!=b positive odd. (7)
   ```

Globally, `k=2`, `|E(Gamma)|=24`, and the existence of a degree-two runner are
equivalent.  At least one degree-two runner is collision-free.  That last
existence statement is only diagnostic: changing to it changes the runner
whose loneliness is being proved.

The closest proved mechanism is THM-2060's one-tail sheet dodge and
THM-2061's exact two-tail dyadic seam.  The canonical hostile is
`V={0,1,...,13}`, whose graph is `K_(7,7)`.  The corrected near miss is the
forbidden inference "some reference has degree two, therefore the original
runner enters (7)"; THM-2888 already isolates this change-of-observer error.
The least-used relevant sidecar is the prescribed runner label itself.

## 2. Proof of the root theorem and the singleton closure

Choose any base velocity `v_*` and put `z_i=(v_i-v_*)/g`.  The gcd of all
`z_i-z_j` is one, so both parities occur.  Moreover

```text
nu_2(v_i-v_j)=s
 iff (v_i-v_j)/g is odd
 iff z_i and z_j have opposite parity.                  (8)
```

This proves `(4)`.  With `v_*=v_r`, equation `(8)` says that the neighbors of
`r` are precisely its odd primitive relatives, proving `(5)`.  Formula `(6)`
is `k(14-k)`.

Now suppose `k=1`.  If the prescribed runner is the unique minority vertex,
all thirteen primitive relatives are odd, and the clock `x=1/2` gives every
relative distance `1/2`.

If the runner lies in the thirteen-vertex majority class, its primitive row
has the form

```text
2C union {w},                    w odd,                 (9)
```

with at most twelve distinct body speeds.  Cited LRC through thirteen total
runners supplies a phase `y` at which every member of `C` has distance at
least `1/13>1/14`.  The two lifts

```text
x_0=y/2,                    x_1=(y+1)/2                 (10)
```

preserve every doubled-body distance.  Since `w` is odd, its two phases
differ by `1/2`; the two open danger arcs of radius `1/14` cannot contain
both.  One lift therefore handles the tail.  This is THM-2060's `q=2`
one-tail mechanism and proves the complete `k=1` closure.

For any prescribed runner, and distinct `i,j` with `i,j!=r`, repeated
absolute relatives occur exactly when

```text
|v_i-v_r|=|v_j-v_r|  iff  v_i+v_j=2v_r.                (11)
```

Such a row has at most twelve distinct positive constraints, so cited LRC
through thirteen total runners gives the stronger threshold `1/13`.  If a
degree-two runner has no collision, `(5)` leaves two distinct odd absolute
relatives and eleven distinct even ones, proving `(7)`.

Finally assume `Gamma=K_(12,2)`.  Take the least and greatest velocities in
the twelve-vertex parity class.  Neither can be the midpoint of two vertices
from that same class.  A collision can therefore use only the two minority
vertices, which have one midpoint; at most one of the two extremes equals
it.  Hence at least one extreme is collision-free.  The example

```text
V={-1,1,0,2,4,...,22}                                  (12)
```

has exactly one collision-free degree-two reference, namely `22`, so the
lower bound one is sharp.  Equation `(12)` does not license a change from a
different prescribed runner.

## 3. Exact projective fixed-pool entry

Retain the fixed pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (13)
```

Let `(7)` be an already anchored, collision-free row.  Suppose there is a
positive reduced rational

```text
lambda=p/q,       gcd(p,q)=1,                             (14)
```

such that

```text
H'=lambda H subset Z_(>0),       |H'\P|<=2.             (15)
```

Then `(7)` is safe.  Equivalently, a counterexample in this literal dyadic
branch must satisfy

```text
|lambda H\P|>=3                                           (16)
```

for every lawful `lambda` in `(14)` for which `lambda H` is integral.  In
particular, `|H\P|>=3` when `lambda=1`.  More generally, an already displayed
row `2cH_0 union {a,b}` cannot be a counterexample when `|H_0\P|<=2`.

### Proof

Integrality in `(15)` and `gcd(p,q)=1` imply `q|h` for every `h in H`; write
`H=qK`, so `H'=pK`.  Multiplication by any positive integer preserves Haar
measure on the circle.  Hence

```text
mu(G_H)=mu(G_K)=mu(G_(H')).                             (17)
```

If `H'` has zero outsiders, THM-4156 gives `mu(G_(H'))>4/63`; with one
outsider THM-4191 gives `mu(G_(H'))>=4/63`; with two outsiders THM-4326 gives
the same bound.  Equation `(17)` transfers that bound to the original body
without changing either tail, and THM-4150 closes `(7)`.

Conversely, every rational rescaling available to this Haar-invariance move
comes from an equality `pH=qH'` of positive integer multiples and therefore
has `H'=(p/q)H`.  Thus `(14)--(15)` exactly characterize **entry through this
projective fixed-pool mass template**, not safety itself.

The test is finite.  If `(15)` holds, at least nine image labels lie in `P`.
For one of them, `lambda=u/h` with `u in P` and `h in H`.  It is therefore
enough to inspect at most `30*11` candidate rational scales, reduce them,
check integrality, and count pool hits.  No tail coordinate is rescaled, so
the numerator parity is irrelevant.

## 4. Parity entry is strictly weaker than pool entry

Consider

```text
V={0,2,4,...,22} union {1,3}.                            (18)
```

The anchors `0` and `22` are collision-free degree-two references.  At either
one the half-body is

```text
H_0={1,2,...,11}.                                       (19)
```

Thus `(18)` passes the affine parity gate.  It does not pass the projective
pool gate.  Since `1 in H_0`, every integral reduced scale in `(14)` has
`q=1`; the exact finite candidate scan finds at most six labels in `P`, with
six attained only at scale `10`.  The required number is nine.  This is a
certificate failure, not an unsafe row.

The rational freedom is nevertheless real.  Put

```text
H_1={1,2,4,5,8,10,15,20,21,30,40}.                     (20)
```

There is no literal fixed-pool entry at the displayed anchor, but the even
numerator scale `lambda=2` gives

```text
2H_1={2,4,8,10,16,20,30,40,42,60,80},                  (21)
```

with nine pool labels and two outsiders.  This is the boundary control that
prevents replacing `(14)` by an integer body divisor test.

## 5. Pair-adaptive Haar transfer

For an eleven-body `H`, put

```text
G_H={y in R/Z:min_(h in H)||hy||>=1/14}.                (22)
```

Write the distinct odd tails as

```text
a=pt,       b=qt,       0<p<q,       gcd(p,q)=1,        (23)
```

so `p,q,t` are odd.  Let `C_(p,q)` be THM-4150's open two-sheet cross-comb.
Then

> **Pair-adaptive criterion.**  If `G_H` is nonempty and
>
> ```text
> mu(G_H)>=mu(C_(p,q)),                                 (24)
> ```
>
> the row `2H union {a,b}` is safe.  Hence every counterexample in this
> branch satisfies the strict, ratio-sensitive obstruction
>
> ```text
> mu(G_H)<mu(C_(p,q))<=4/63.                            (25)
> ```

For an eleven-body, the nonemptiness in `(24)` is automatic: cited
lower-dimensional LRC supplies a point with clearance at least `1/12`, and
continuity supplies a nonempty open subset at the weaker threshold `1/14`.

Indeed, failure gives exactly the THM-4150 containment

```text
G_H subset m_t^(-1)(C_(p,q)),             m_t(y)=ty.    (26)
```

The right side is a proper open set of measure `mu(C_(p,q))`, while `G_H` is
nonempty and compact.  Strict measure excess contradicts inclusion.  At
equality, the open difference would have measure zero and hence be empty, so
the two sets would form a nonempty proper clopen subset of the circle, again
impossible.  This proves `(24)--(25)`.

The exact threshold is executable from THM-4150:

```text
mu(C_(p,q))
 =2/49+2[B_2({1/2+(q-p)/14})-B_2({1/2+(q+p)/14})]/(pq), (27)
```

where `B_2(u)=u^2-u+1/6` on `[0,1)`.  Its maximum is `4/63`, uniquely at
`(p,q)=(1,9)`; it is zero at `(1,3)` and `(1,5)`.  Thus `(24)` is strictly
stronger than the uniform THM-4150 test on every submaximal ratio.

## 6. Connection contract and scope firewall

```text
source:       labelled affine velocity configuration with prescribed runner
target:       minimum-v2 cut graph, then a doubled body with odd tails
map:          divide affine content; take the prescribed anchor star;
              identify equal absolute relatives only after recording them
preserved:    anchor-relative LRC predicate, parity, common dilation class
destroyed:    magnitudes and signs after the parity quotient; the graph alone
              also forgets midpoint collisions and fixed-pool membership
sidecar:      anchor label, signed relatives, midpoint flag, H, and p/q
decisive test: anchor degree; midpoint equation; finite u/h pool scan;
               actual-ratio Haar comparison
```

Scope firewalls:

1. Re-referencing does not preserve the distinguished runner.  The existence
   of one collision-free degree-two anchor proves nothing for another anchor.
2. The `k=1` result is a complete root-type closure, but the theorem does not
   close `k=2,...,7`.  In type `k=2`, the two degree-twelve anchors are not in
   the two-tail lane.
3. Conditions `(14)--(15)` characterize entry through the stated rational
   fixed-pool mass template.  Their failure is not danger, and pool labels
   are not a projective invariant.
4. The scalar condition `(24)` is sufficient, not necessary.  When it fails,
   THM-2061's component locations and dyadic owner words remain load-bearing.
5. No arbitrary-row entry, owner/arrival transfer, or proof of LRC(14) follows.

The successful research move is to pay the affine rank first and then audit
the physical action: the quotient graph exposes the exact parity type, while
the anchor and projective-scale sidecars prevent a change-of-observer or
confusion between body-mass equivalence and physical tail rescaling.

## 7. Exact probe and replay

The probe is corroborative; the complete-bipartite proof above is
symbolic.  It freezes the following hostile and positive controls:

- all seven root sizes and edge counts in `(6)`;
- all `binom(18,14)=3,060` bounded configurations in `{-2,...,15}`;
- all `13*binom(12,2)=858` configurations in the structured split universe;
- the `K_(7,7)` AP hostile, the sharp one-reference example `(12)`, and the
  parity-without-pool example `(18)`;
- literal fixed-pool and even-numerator projective positive entries;
- translation and nonzero common-dilation invariance on 90 images.

Replay:

```bash
python3 04-computation/lrc14_entry_parity_affine_classification_probe.py \
  | diff -u \
  05-knowledge/results/lrc14_entry_parity_affine_classification_probe.out -

python3 -O 04-computation/lrc14_entry_parity_affine_classification_probe.py \
  | diff -u \
  05-knowledge/results/lrc14_entry_parity_affine_classification_probe.out -

PYTHONHASHSEED=1729 \
python3 04-computation/lrc14_entry_parity_affine_classification_probe.py \
  | diff -u \
  05-knowledge/results/lrc14_entry_parity_affine_classification_probe.out -
```

All three commands pass with `935,790` assertions.  The script and frozen
output hashes are recorded in the frontmatter.
