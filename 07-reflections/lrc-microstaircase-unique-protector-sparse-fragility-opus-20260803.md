# Unique-protector cells and boundary moments force scalar fragility

**Status: PROVED ELEMENTARY LEMMAS + FINITE-EXACT AUDIT; PROMOTED IN
[THM-3317](../01-canon/theorems/THM-3317-unique-protector-cells-and-weighted-scalar-fragility.md).**
For every `n>=3` and every coordinate, this note constructs
an explicit micro-staircase cell whose zero-ramp endpoint is protected by that
coordinate alone.  It then proves a weighted Hamming-ball certificate around
every scalar ramp.  At odd prime moduli `p>=7`, a finite-field moment argument
also excludes the first Hamming layer beyond that certificate.  The later
all-support interpolation argument in
[THM-3316](../01-canon/theorems/THM-3316-prime-right-boundary-interpolation-forces-scalar-rigidity.md)
subsumes that prime corollary and proves rigidity at every prime modulus.  This
is progress on the finite residue-cell problem, not a speed-to-residue lift
and not LRC(14).

The matching exact companion is
[`lrc_microstaircase_unique_protector_sparse_fragility_opus_20260803.py`](../04-computation/lrc_microstaircase_unique_protector_sparse_fragility_opus_20260803.py),
with frozen
[`output`](../05-knowledge/results/lrc_microstaircase_unique_protector_sparse_fragility_opus_20260803.out).

## Inheritance pass and board

- **Closest proved mechanism:** [THM-364](../01-canon/theorems/THM-364-lrc-scalar-ramp-cell-blocking.md)
  identifies scalar ramps with the Dirichlet equality spine, and
  [THM-363](../01-canon/theorems/THM-363-lrc-scalar-gauge-reindexing.md)
  makes scalar addition a genuine cell reindexing.
- **Canonical hostile:** the normalized `n=14` coordinate-six half-turn misses
  only `56/11368` candidates.  It shows that a one-coordinate deformation can
  be extremely close to a blocker.
- **Corrected near miss:** the right-adjacent rank-two affine cover in the
  independent HYP-1823 audit is necessary but not sufficient; the vector
  `(0,1,0,...)` passes every boundary test while failing on interior cells.
- **Least-used sidecar:** the complete set of zero-ramp endpoint owners of a
  cell.  Retaining the owner set, rather than only the fact that Dirichlet
  blocks the cell, exposes a cell on which every zero perturbation coordinate
  becomes irrelevant.

The live concepts were the scalar-gauge orbit, the endpoint-owner set, the
shift congruence fibres, and the physical lift.  The first three produce the
lemma below; the fourth remains absent.

## Explicit unique-protector cell

Fix `n>=3` and `1<=j<=n-1`.  Define

```text
                         1      1
t_(n,j) =              ----- + ----     if 2j<=n-1,
                         2     4n^2

t_(n,j) = 1/(2n)                         if 2j>n-1,

alpha_(n,j) = (1-t_(n,j)/n)/j.                         (1)
```

**Lemma 1.**  `alpha_(n,j)` is not a breakpoint of the full arrangement and

```text
floor(n*{i alpha_(n,j)}) in {0,n-1}  iff  i=j.          (2)
```

Thus every coordinate has an explicit open cell on which it is the unique
zero-ramp protector.

### Proof

Write `i=qj+r`, `0<=r<j`.  Then

```text
i alpha_(n,j)=q+r/j-i t_(n,j)/(nj).                    (3)
```

For `i=j`, the fractional part is `1-t/n`, strictly inside the last bin.

Suppose first that `2j<=n-1`.  For a larger multiple `i=qj`, `q>=2`, its
distance below the next integer is `qt/n>1/n`, while `qt<n-1`; hence it is
interior.  For `r>=1`, first put `t=1/2`.  The lower-bin margin has numerator

```text
2nr-i-2j >= 2n-(n-1)-(n-1)=2.                          (4)
```

The extra `1/(4n^2)` changes the fractional part by less than that strict
margin.  The upper margin follows from `r/j<=1-1/j<1-1/n`.

If `2j>n-1`, there are no larger multiples of `j`.  At `t=0`, every other
coordinate is at distance at least `1/j>1/n` from an integer.  The perturbation
`t=1/(2n)` moves it by less than half the smallest lower margin and only away
from the upper endpoint.  All inequalities are strict, so `(1)` is not a
breakpoint and `(2)` follows.  QED

## Weighted scalar-fragility certificate

For a nonzero residue `a mod n`, define its shift cost

```text
kappa_n(a) = 2                  if gcd(a,n)=1,
             gcd(a,n)           if gcd(a,n)>1.           (5)
```

Let `v=(v_1,...,v_(n-1))` be any residue vector and choose any scalar ramp
`mi`.  Put

```text
w_i=v_i-mi mod n,       J={i:w_i!=0}.                    (6)
```

**Lemma 2 (sparse fragility).**  If `J` is nonempty and

```text
sum_(i in J) kappa_n(w_i) < n,                          (7)
```

then `v` is not a full micro-staircase blocker.

### Proof

THM-363 permits replacing `v` by `w`.  Choose any `j in J` and use the cell
from Lemma 1.  If `i notin J`, then `w_i=0` and its bin is interior by `(2)`,
so it blocks no shift.

For `i in J`, its bad shifts solve one of the two congruences

```text
s w_i = -b_i       or       s w_i = -b_i-1  mod n.      (8)
```

Put `g=gcd(w_i,n)`.  If `g=1`, each congruence has one solution and the two
solutions are distinct.  If `g>1`, two consecutive right sides cannot both be
divisible by `g`; at most one congruence is soluble, and it then has exactly
`g` solutions.  Coordinate `i` therefore forbids at most `kappa_n(w_i)`
shifts.  The union bound and `(7)` leave a shift forbidden by no coordinate.
That shift and the unique-protector cell are a safe candidate.  QED

The statement is gauge invariant: it certifies a weighted Hamming
neighbourhood of **every** scalar ramp, not only the normalized zero ramp.

## Concrete consequences

- For prime `p`, every perturbation supported on at most `(p-1)/2`
  coordinates is safe.  In particular, the `n=13` classification is already
  proved on every normalized vector of support at most six before invoking an
  exact classifier.
- At `n=14`, every nonzero residue except the half-turn `7` has cost `2`,
  while `kappa_14(7)=7`.  Hence the lemma covers up to six ordinary defects,
  or one half-turn plus up to three ordinary defects.
- Every one-coordinate deformation of every scalar ramp is safe for every
  `n>=3`, including the historically hardest coordinate-six half-turn.  The
  theorem explains its `56` misses without needing their enumeration, though
  it does not recover the exact count.

After the exact gauge `w_1=0`, a dynamic-programming count gives the following
numbers of nonzero quotient classes certified by `(7)`:

```text
n=13: 1,501,621,044
n=14: 2,970,296,364
n=15: 7,588,647,950.                                   (9)
```

These large absolute counts are small fractions of the dense quotient spaces;
the point is the uniform mechanism and its near-scalar scope, not a global
classification claim.

## The first prime layer also fails

The next layer admits a second, quite different obstruction.

**Lemma 3 (prime critical-layer obstruction).**  Let `p>=7` be an odd prime.
If a full micro-staircase blocker is not a scalar ramp, then its Hamming
distance from every scalar ramp is at least

```text
(p+3)/2.                                                   (10)
```

### Proof

By Lemma 2, a nonzero scalar-gauge representative cannot have support at most
`(p-1)/2`.  Suppose its support `S` has the next possible size

```text
h=(p+1)/2.                                                 (11)
```

Use the open cell immediately to the right of `alpha=a/p`.  Its pattern is
`b_i=ai mod p`.  A zero coordinate blocks all shifts only when `a=0` or
`a=-i^(-1)`.  Therefore, for each `j in S`, setting `a=-j^(-1)` makes every
zero coordinate interior.  The nonzero coordinates must then cover all
`s in F_p` by the `2h=p+1` incidences

```text
s in { l/(j w_l), (l-j)/(j w_l) },   l in S.              (12)
```

After multiplying shifts by `j` and writing `q_l=w_l^(-1)`, these are the
pairs `{l q_l,(l-j)q_l}`.  Each pair has two distinct elements.  A cover by
`p+1` incidences consequently contains every field element once and one
element `delta_j` twice.

For `1<=m<=h-1`, the sum of the `m`-th powers of all elements of `F_p` is
zero.  Hence

```text
sum_(l in S) ((l q_l)^m+((l-j)q_l)^m) = delta_j^m.         (13)
```

At `m=1`, put `A=sum l q_l` and `B=sum q_l`; then
`delta_j=2A-jB`.  For general `m`, subtract `(2A-jB)^m` from the left side of
`(13)`.  This is a polynomial in `j` of degree at most `m`, vanishing at all
`h` elements of `S`.  Since `m<h`, it is identically zero.  Its leading
coefficient gives

```text
sum_(l in S) q_l^m = B^m,   1<=m<=h-1.                   (14)
```

If `B=0`, Newton's identities say that the first `h-1` elementary symmetric
functions of the nonzero `q_l` vanish.  Their root polynomial is therefore
`T^h-c`.  It has at most `gcd(h,p-1)<=2` nonzero roots, all simple, and cannot
split into `h>=4` linear factors.

If `B!=0`, replace every `q_l` by `q_l/B`.  Newton's identities applied to
`(14)` give root polynomial

```text
F(T)=T^h-T^(h-1)+(-1)^h E,                               (15)
```

where `E` is the product of the normalized roots.  The root recurrence from
`(15)`, together with `sum q_l^m=1` for `1<=m<h` and Fermat's
`sum q_l^(p-1)=h`, forces `(-1)^h E=h-1`.  Thus

```text
F(T)=T^h-T^(h-1)+h-1.                                   (16)
```

For nonzero `T`, Euler's criterion turns `(16)` into

```text
chi(T)(T-1)+h-1=0.
```

So a square root must be `T=2-h`, while a nonsquare root must be `T=h`.
There are at most two distinct field roots.  Moreover

```text
F'(T)=T^(h-2)(hT-(h-1)),
```

so at most one of them can be double and neither can have larger
multiplicity.  The total field-root multiplicity is at most three, strictly
less than `h`.  Again `(15)` cannot split into its `h` prescribed linear
factors.  This contradiction excludes `(11)`, proving `(10)`.  QED

At `p=13`, this rules out support seven: `(16)` is

```text
T^7-T^6-7,
```

whose only nonzero `F_13` root is `7`, and it is simple.  Combined with the
independent exact prime-unit classifier, scalar ramps are now the only full
blockers at `n=13`; the present human argument disposes of supports one
through seven, while the classifier closes the dense supports eight through
eleven.

THM-3316 subsequently removes the support split altogether: the same
right-boundary family forces `l -> l-yw_l` to permute the support off a small
collision set, and a degree-`<=p-2` finite-field sum rules out every nonzero
support.  Lemma 3 remains useful as an independent moment/Newton view of the
first critical layer, but it is no longer the strongest prime result.

## Equality boundary and hostile control

The strict inequality in `(7)` is sufficient, not necessary.  At `n=14`, put
half-turn defects at coordinates `2` and `3`.  Its cost is exactly

```text
kappa_14(7)+kappa_14(7)=14,                             (17)
```

so Lemma 2 makes no claim.  Direct exact checking nevertheless finds a safe
candidate at shift `1` in cell `254`, with pattern

```text
(4,8,13,3,7,12,2,7,11,1,6,10,1).                      (18)
```

Thus equality does not mean blocker.  At primes, Lemma 3 resolves the first
layer beyond equality by combining all right-boundary cells.  At composites,
covering all shifts at equality would still require the bad-shift cosets to
partition `Z/nZ` without overlap.  Varying the unique-protector cell can force
overlaps, suggesting a laminar/coset refinement of the union bound.

## Connection and loss ledger

```text
source:      scalar-ramp gauge orbit in the full open-cell system;
target:      one explicit endpoint-owner cell plus bad-shift cosets;
map:         subtract a scalar ramp, then choose alpha_(n,j);
preserved:   existence of a safe (shift,cell) candidate;
destroyed:   width, prime-grid realization, endpoint/gcd lift and physical time;
sidecar:     the unique endpoint owner j and gcd of each perturbation residue;
test:        sum kappa_n(w_i)<n.
```

This is a legitimate finite residue-cell theorem.  It does not turn a speed
tuple into `w`, guarantee a cell of usable width on a verifier grid, or handle
the endpoint/divisibility branches discarded by the model.

## Audit and next move

The companion checks all explicit cells for `3<=n<=30`, every congruence
cardinality for `3<=n<=60`, and every certified normalized vector for
`3<=n<=8`.  It audits the critical polynomial for every prime from `7` to
`97`, exhausts the `F_13` moment multisets, and checks the scalar positive
controls and equality hostile `(17)--(18)`.  Normal and optimized runs are
intended to have identical semantic output; the frozen transcript records the
exact hashes.

Reproduce from the repository root with

```powershell
python 04-computation/lrc_microstaircase_unique_protector_sparse_fragility_opus_20260803.py
python -O 04-computation/lrc_microstaircase_unique_protector_sparse_fragility_opus_20260803.py
```

Both modes emit the frozen transcript exactly.  SHA-256 (LF-normalized, which
equals the raw-file hash in this worktree) is

```text
script: 4ab6544d589f90263ecbeb4fc39b175d9cb788d406de911c76deeb45d87b96f0
output: f47f93daae5e9ebef79098e133fe00cc2ce6677ddcaa0d8a320166d01950eb9f
```

The prime degree-three excess layer is closed by THM-3316.  The next human
target is composite: classify when the cosets in `(8)` can partition all
shifts, retaining CRT valuations and the collision-set sidecar that the prime
field proof uses.  Every physical LRC lift remains a separate obligation.
