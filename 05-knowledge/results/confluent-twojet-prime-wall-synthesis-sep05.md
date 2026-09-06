# Full-residue two-jet Smith wall: proof candidate and exact audit

**Status: PROVED CANDIDATE / awaiting independent mathematical audit and
canon promotion.** This note closes the `s=p` formula left open by
[THM-4080, confluent two-jet single-scale Smith partition](../../01-canon/theorems/THM-4080-confluent-two-jet-single-scale-smith-partition.md).
The theorem identifier is deliberately unreserved here. Claims marked below
are proved in this note; the separate computational census is FINITE-EXACT.

**Concurrent-result integration (2026-09-05):** incoming commit `566677ae1d`
already proves precisely the `s=p` partition and `n<=p^2` corollary in
[the complete-residue Smith boundary](synthesis_20260905_wildcard_smith_boundary.md),
with a separate [analytical audit](synthesis_20260905_lrc_audit_smith_boundary.md).
Both incoming arguments were read and accepted against the present proof.
Sections 1--2 here are an independent reproduction of that same closure;
they must not be counted as another new theorem. Formula `(5)` makes the
incoming row-saturation step explicit as an integral divided trace. The
additional results of this note are Section 3's exact precision/kernel and
tensor-grid corollaries, and Section 4A's all-depth dyadic three-node law and
the failure of an unchanged Smith prefix under same-residue adjoining.
The two audit banks use different implementations and residue lifts, so both
remain useful reproducibility evidence.

## Inheritance and concept board

- Closest proved mechanism: THM-4080's weighted minors and residual
  confluent-evaluation matroid. Its scope is `s<p`.
- Canonical hostile: at `(p,e,s)=(2,1,2)` the Smith exponents are
  `(0,0,2,2)`, not the false extension `(0,0,1,3)`.
- Corrected near miss: [THM-4010, confluent consecutive Hasse observer](../../01-canon/theorems/THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall.md)
  separates the Hermite--Newton determinant from the false repeated-factorial
  Smith diagonal. Its cluster decomposition remains load-bearing.
- Least-used sidecar: THM-4010's higher Bockstein/saturation filtration.
  The residue matroid loses exactly one derivative direction, but this does
  not say how many powers of `p` that direction costs to recover.
- Live concepts: weighted minors; the derivative trace divided by `p`;
  residue-cluster CRT; finite-precision observers; tensor-product jet grids.

The relevant META-PATTERNS move is "Normalize repeated-response towers
before scalarizing": preserve the lost layer, rather than treating its
mod-`p` disappearance as total loss. A second useful move is "Divide
exceptional multiplicity before judging a wall." No new general method card
is proposed from this one thread.

## 1. Exact wall formula

Let `p` be any prime and let `e>=1`. Choose `p` integers `u_0,...,u_(p-1)`
whose residues are pairwise distinct modulo `p`, and any `a in Z`. Put

```text
x_i=a+p^e u_i.
```

The two-jet matrix on `Z[X]_<2p`, in monomial coordinates, has rows

```text
(P(x_i),P'(x_i)),             0<=i<p.
```

Here `P'=D^[1]P`; no higher ordinary derivative is substituted for a Hasse
jet. Its increasing `p`-Smith exponent list is

```text
(0,0,e,2e,...,(p-2)e,
 (p-1)e+1,
 (p+1)e,(p+2)e,...,(2p-2)e,
 (2p-1)e-1).                                         (1)
```

Empty ranges are omitted. In particular, `p=2` gives
`(0,0,e+1,3e-1)`, agreeing with THM-4010's exact two-node formula.
All inequalities in `(1)` hold even at `p=2,e=1`, where the last two
exponents agree. For `e=1` the largest two exponents agree for every prime.
The sum is `2e p(p-1)`, the confluent determinant valuation.

Equivalently, writing `d_h=v_p(Delta_h)` for determinantal divisors,

```text
d_0=0,
d_h=e(h-1)(h-2)/2,                    1<=h<=p,
d_h=e(h(h-1)/2-p)+1,                  p<h<2p,
d_(2p)=2e p(p-1).                                  (2)
```

Thus the wall changes the old tariff by exactly one on every middle rank.
It does not multiply that correction by the collision depth `e`.

### Weighted lower bounds

Translate `a` to zero by an integral binomial unitriangular coefficient
change. For an `h`-minor, let `q_1<...<q_h` be the column degrees and let
`d` be its number of derivative rows. Factoring scales gives

```text
v_p(minor)=e(sum q_j-d)+v_p(residual minor).          (3)
```

The residual entries are `u_i^q` and `q u_i^(q-1)`. THM-4080's elementary
column-degree bounds apply unchanged. For `h<=p`, its witness has one value
row and `h-1` derivative rows. Its residual determinant is `(h-1)!` times a
Vandermonde unit, so the first line of `(2)` is attained, including `h=p`.

For `p<h<2p`, the old lower tariff is

```text
T_h=e(h(h-1)/2-p).
```

Any minor not using both the first `h` columns and all `p` derivative rows
has scale tariff at least `T_h+e>=T_h+1`. For the remaining minors, the
residual derivative rows have rank `p-1` modulo `p`, so every such residual
minor is divisible by `p`. This proves `d_h>=T_h+1`.

To justify the rank assertion explicitly, summing the derivative rows gives
zero modulo `p` on every degree `q<2p-1`:

```text
sum_i q u_i^(q-1)=0 mod p.                          (4)
```

For `q=0,p` this is immediate. Otherwise the power sum over all `F_p`
vanishes: within `0<=q-1<=2p-3`, the only potentially nonzero positive
multiple of `p-1` is `q-1=p-1`, whose prefactor is `q=p`. The exponent-zero
sum is `p=0`. Meanwhile the first `p-1` derivative rows on columns
`1,...,p-1` have unit determinant `(p-1)!` times a Vandermonde, proving
rank exactly `p-1`.

### One divided trace attains every middle lower bound

Fix `p<h<2p`, and write `D_i` for the residual derivative row on columns
`0,...,h-1`. Replace the last row by the integral divided trace

```text
B=(D_0+...+D_(p-1))/p.                              (5)
```

Integrality follows from `(4)` for the actual integer lifts, not just for
canonical residues. The determinant of this row transformation is `1/p`.
The `p` transformed derivative rows on columns `1,...,p` have determinant

```text
(p!/p) product_(i<j)(u_j-u_i),                       (6)
```

a `p`-adic unit. They therefore have rank `p` modulo `p`.

All value rows together with the original derivative rows have rank `h`
on polynomials of degree below `h`: a common kernel polynomial has `p`
distinct double roots and hence is divisible by a degree-`2p` polynomial.
The original last derivative row is `pB-sum_(i<p-1)D_i`, so adjoining the
value rows to the transformed rows still gives rank `h` modulo `p`.
Select `h-p` value rows extending the transformed derivatives to a basis.

The resulting transformed `h`-minor is a unit; undoing `(5)` multiplies
its determinant by exactly `p`. We have therefore exhibited a literal
original residual minor of valuation exactly one. Restoring the scale in
`(3)` attains `T_h+1`.

Finally, the full confluent determinant is
`product_(i<j)(x_j-x_i)^4`, whose valuation is `2e p(p-1)`. Taking consecutive
differences in `(2)` proves `(1)`.

The mechanism can be seen on the missing Frobenius degree: modulo `p`,
`(X^p)'=0`, whereas `(1/p)sum_i (u_i^p)'=sum_i u_i^(p-1)` is a unit.
The actual lost direction survives in the integral lattice one level later.

## 2. Consecutive nodes through the complete quadratic band

For `n` consecutive nodes and `1<=n<=p^2`, write

```text
s_c=number of nodes congruent to c modulo p.
```

Each nonempty residue cluster has size at most `p` and exact internal
difference valuation one. THM-4010's local CRT makes the global exponent
partition the sorted multiset union of the cluster partitions. Use
THM-4080's formula when `s_c<p`, and `(1)` with `e=1` when `s_c=p`.
This extends the previously closed range `n<=p(p-1)` to `n<=p^2`.

The new boundary `n=p^2+1` really changes the object: one cluster has two
nodes whose difference has valuation two. It does not satisfy the exact
one-scale hypothesis. No extrapolation or closure for arbitrary larger
clusters is claimed.

## 3. Exact observation precision and a tensor transfer

Let `L=(2p-1)e-1`, the largest exponent in `(1)`, and work over `Z_p`.
For every `N>=1`, knowledge of all two-jet observations modulo `p^(N+L)`
determines every coefficient modulo `p^N`. This precision loss is optimal
uniformly over source polynomials. Indeed, Smith transformations are
invertible over `Z_p`, and a last-coordinate Smith basis vector gives an
integral primitive coefficient vector whose entire observed vector has
valuation exactly `L`. Multiplying that vector by `p^(N-1)` proves that
`N+L-1` digits do not suffice.

The observer reduced modulo `p^N` has kernel cardinality

```text
p^(sum_i min(N,alpha_i)),                            (7)
```

where the `alpha_i` are `(1)`. Thus `(1)` measures all finite precision
losses, and does more than give the determinant or first nonzero layer.

For `d` independent coordinate grids satisfying this theorem at depths
`e_1,...,e_d`, take polynomials of coordinate degree below `2p` and all
rectangular Hasse jets of orders zero or one at the product grid. The
matrix is the tensor product of the one-variable matrices. Tensoring their
invertible Smith transformations shows that the new exponents are exactly
all sums of one exponent from each coordinate. The sharp uniform recovery
loss is therefore

```text
(2p-1)(e_1+...+e_d)-d.                              (8)
```

This is an actual precision extension of
[THM-4255, specialization kernel and transverse Hasse-jet repair](../../01-canon/theorems/THM-4255-specialization-kernel-and-transverse-hasse-jet-repair.md),
Section 4.1, and of the independent-coordinate jet observer in
[THM-4263, moving multigraph filtered jets](../../01-canon/theorems/THM-4263-moving-multigraph-filtered-jet-and-finite-factor-density-transport.md),
Section 1. It is restricted to constant integral grids and the stated degree
box. Moving graphs or restricted arithmetic Hom lattices require their
actual source module and cannot inherit `(8)` from matching dimensions.

## 4. Connection and failure ledger

| Source -> target and map | Preserved predicate | Lost information / required sidecar | Decisive test |
|---|---|---|---|
| Integral two-jet lattice -> weighted residual rows, scale extraction | Determinantal valuations via `(3)` | Residue rank alone drops valuation depth; retain the divided trace `(5)` | The unit `(p-1)!` Vandermonde in `(6)` |
| Consecutive observer -> residue clusters, local CRT | Full `p`-Smith partition | Other primes and internal deeper collisions; retain all cluster sizes and depths | `n=p^2+1` leaves one-scale scope |
| Smith lattice -> finite-precision coefficient observer | Exact source congruence recovery | A single scalar specialization still has its original kernel; retain degree box and all jets | Primitive last Smith coordinate |
| Independent coordinate observers -> rectangular grid, tensor product | All exponent sums and sharp precision loss | Coupled node geometry is absent; retain independent product grid | Direct small product-grid Smith audit |

The old false extension loses the divided-trace coordinate. Its strongest
survivors are the total determinant, the zero layer, and all ranks `h<=p`.
Its repaired form is `(1)`, not a universal "scale the depth-one answer by
e" rule: e.g. `(p,e)=(3,2)` gives `(0,0,2,5,8,9)`, which is not twice
`(0,0,1,3,4,4)`.

No LRC(14), Rule-30, Jacobian, external p-adic-zeta, or Gaussian-moment
consequence is asserted. General jet order and arbitrary clusters with deeper
collisions remain open. The following first dyadic corner is explicit.

## 4A. A complete three-node dyadic collision formula

For nodes `(0,2^e,2^(e+1))`, with `e>=0`, the two-jet exponents are

```text
e=0: (0,0,0,0,2,2),
e=1: (0,0,2,2,5,7),
e>=2: (0,0,e+1,2e+1,4e,5e+2).                       (9)
```

This is an all-depth statement for one specified family with two collision
scales. It is independent of translating the three nodes by any integer.

To prove it, put `d=2^e` and clear the first node's value and derivative
identity block. The remaining four-by-four matrix has rows `(u,r)` equal to
`(1,0),(1,1),(2,0),(2,1)` and columns `q=2,3,4,5`, with entries

```text
q^r u^(q-r) d^(q-r).
```

For any `h`-minor, set `w=sum q-sum r`. Its determinant is `d^w` times
an integer residual determinant. The complete finite table below gives
`w:minimum v_2(residual determinant)` among nonzero minors of that weight:

```text
h=1: 1:1, 2:0, 3:0, 4:0, 5:0
h=2: 3:2, 4:0, 5:1, 6:0, 7:1, 8:0, 9:4
h=3: 7:2, 8:2, 9:2, 10:2, 11:3
h=4: 12:4.                                         (10)
```

This table covers all `16+36+16+1=69` residual minors; both the original
integer determinant path and the independent SymPy path reconstruct it.
For all real `e>=0`, its lower envelopes are

```text
delta_1=min(e+1,2e),
delta_2=min(3e+2,4e),
delta_3=7e+2,
delta_4=12e+4.                                     (11)
```

Each omitted tariff is dominated by one displayed tariff for every `e>=0`.
For integral `e>=0`, taking successive differences of `(11)` and prepending
the two unit exponents proves `(9)` without interpolation in `e`.

This also supplies a sharp warning for the next-scale recursion. The pair
`(0,8)` has exponents `(0,0,4,8)` by the proved wall formula. Adjoining `16`
gives `(0,0,4,7,12,17)`: an old exponent decreases. Decomposing node subsets
inside one residue class therefore cannot be treated as a direct sum or as
an unchanged Smith prefix. CRT requires unit resultants between the actual
factors. Its repaired carrier is the full weighted minor filtration, which
retains interactions between collision scales.

## 5. FINITE-EXACT evidence

Reproduce with:

```bash
python3 04-computation/confluent_twojet_prime_wall_synthesis_sep05.py
python3 -O 04-computation/confluent_twojet_prime_wall_synthesis_sep05.py
python3 04-computation/confluent_twojet_prime_wall_synthesis_sep05_independent.py
python3 -O 04-computation/confluent_twojet_prime_wall_synthesis_sep05_independent.py
```

The standalone standard-library script checks:

- 30 complete-wall DVR matrices for `p=2,3,5,7,11`, `e=1,2,3`, two residue
  lift systems, including integer translations;
- inherited strict-`s<p` positive controls;
- four independent complete determinantal-divisor reconstructions, using
  every minor, for `p=2,3`, `e=1,2`, with noncanonical lifts;
- 46 literal rank-by-rank residual witness minors after one divided trace;
- ten direct consecutive-node matrices in the newly added band;
- 1,027 consecutive formula rows through prime 19, checking rank and the
  independently computed factorial determinant valuation;
- the repeated-residue hostile `(nodes,p)=((0,2,4),2)`, with actual
  exponents `(0,0,2,2,5,7)`;
- three direct tensor-grid Smith matrices and complete finite-precision
  kernel censuses for `(p,e)=(2,1)` at precisions `N=1,2,3`;
- nine dyadic triple matrices at `e=0,...,8`, including the non-preserved
  prefix hostile `(0,8)->(0,8,16)`.

The independent script imports no primary code. It uses SymPy's integral
Smith normal form for twelve further translated noncanonical complete-wall
systems, all nine dyadic triples, and all 69 symbolic residual minors behind
the all-depth tariff table. Root independently reviewed the wall proof and
replayed twelve SymPy wall cases before this checker was added. Normal and
optimized output digests are frozen in the companion outputs at checkpoint.

Frozen SHA-256 values (raw LF bytes):

```text
1af12f22e5a64512910e547ed2dff3e543228fa6b718936682402d93c34f571f  primary script
cee4805d74339d6071dcf25bb3bd3f7fa9d03605e8c5e1bb6efc6c2fdaee367d  primary output
3eb56ee199fc8f0040e4f8ebec95497f64830436bbfbf3f80385370d24c3e3e6  independent script
dfecb937649be3fc56e3770f1fbcef8403f0f495d07827f69014e737dad5b264  independent output
```

Both normal/optimized output pairs byte-match. The final primary run has
441,971 exact gates. The semantic digests are
`9401f0259daa64e6a5d63c01e111d984717f42cb61e6c6bf18f7f83e8f4439c0`
and `396d73ba5764eb359f82b00c637295c8cf38a41b339b5fde817e1aca4d760444`.
