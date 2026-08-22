---
id: THM-3168
title: "Live-bank projected pole transition and sharp depth-two rank wall"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT THEOREM AUDIT.
  For the genuine support-(1,3), bank-I2 selector currents, adding any fixed
  physical pole value admits a linear parent-to-child map on the finite
  projected span of empty and one-pole parents in every degree 5 through 13.
  The extension fails sharply when two-pole parents are admitted: already in
  degree 5, for every pole value, the parent and child row spaces both have
  rank 6 and their stack has the maximal rank 12.  A seven-pair normalized
  determinant gives a direct exact obstruction for pole 1.  These are finite
  span statements, not positive, canonical, Markov, or stopping transitions.
source: root/multiscale-newton-flag/low-child-flag-extension-2026-08-02
depends_on:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
related:
  - THM-3158-sharp-depth-five-selector-resurrection-and-degree-thirteen-death-barcode
  - THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary
script: 04-computation/gmc_live_bank_projected_pole_transition_thm3168.py
output: 05-knowledge/results/gmc_live_bank_projected_pole_transition_thm3168.out
script_sha256: 59368b5e08f3654bc7f2da4849f1bad2b3e8d23d6023a206a59a6ec185fa4e16
output_sha256: 5f210a22561cd8d798f07264fdeb5a6a55645c6b24441d5f3cc6faf3011778f6
hash_basis: LF-normalized bytes
---

# THM-3168 -- live-bank projected pole transition and sharp depth-two rank wall

**PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT THEOREM
AUDIT.**

THM-3160 proves that same-degree selector currents do not carry a universal
pole-subtraction action.  Its counterexamples use valid endpoint functionals,
but not states of the live support-`(1,3)`, bank-`I2` selector bank.  The
present theorem closes that scope gap.  A finite-span transition exists on
the shallowest live prefixes, then fails as strongly as possible as soon as
two-pole parents are included.

The survival at depth one is only interpolation on a finite span.  It is not
a positive kernel or a local rule.  The depth-two failure is the substantive
statement: the projected child data add a second independent copy of the
entire degree-five zero-mass space.

## 1. Fixed-pole parent and child matrices

Fix support `(a,b)=(1,3)`, invariant `I=1`, and bank `I2`.  The physical pole
multiset is

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1).                    (1)
```

Include the empty prefix.  For a pole value `M in {1,...,8}` and depth cap
`d`, let

```text
E_d(M)={sigma: |sigma|<=d and sigma union {M} is a legal
                multiplicity-valid submultiset of P}.       (2)
```

For `sigma in E_d(M)`, write

```text
p_sigma=G_N^sigma,             c_sigma=G_N^(sigma union {M}) (3)
```

as columns indexed by the partitions of `N`.  Every current has total mass
zero, so delete the redundant `(1^N)` row.  Let `P_(M,N,d)` and `C_(M,N,d)`
be the resulting parent and child matrices.

The elementary transition criterion is

```text
there is L:span{p_sigma}->Q^(p(N)-1) with Lp_sigma=c_sigma

iff ker P_(M,N,d) subset ker C_(M,N,d)

iff rank(P_(M,N,d))=rank([P_(M,N,d);C_(M,N,d)]).            (4)
```

Indeed, the assignment on generators is well-defined exactly when every
parent relation remains a child relation.  This criterion concerns a linear
map on the displayed finite span; it says nothing about positivity or about
values away from that span.

## 2. Exact shallow survivor through degree thirteen

At parent depth at most one, the pair counts are

```text
|E_1(M)|=9 for M=1,...,5,          |E_1(M)|=8 for M=6,7,8. (5)
```

The difference occurs because values one through five have repeated physical
copies in `(1)`, whereas six through eight do not.  Exact rational rank and
kernel checks give:

| pole values | `N=5` | `N=6` | `N=7,...,13` |
|:---|:---:|:---:|:---:|
| `M=1,...,5` | `6=6` | `8=8` | `9=9` |
| `M=6,7,8` | `6=6` | `8=8` | `8=8` |

Each entry is

```text
rank(P_(M,N,1))=rank([P_(M,N,1);C_(M,N,1)]).               (6)
```

For `N=5,6`, the companion computes the exact rational parent kernels and
checks that every basis vector is killed by the child matrix.  For `N>=7`,
the parent columns are independent; a nonzero minor modulo the prime `65521`
is an exact certificate of full column rank over the rationals.  Criterion
`(4)` proves that a finite-span linear transition exists for every entry of
the table.

The high-degree entries are deliberately not overinterpreted: when the
parents are independent, any assigned child columns interpolate linearly on
their span.

## 3. The sharp live depth-two wall

Adjoin every legal two-pole parent and specialize to degree five.  The exact
census and ranks are

| `M` | `|E_2(M)|` | `rank P` | `rank C` | `rank[P;C]` |
|---:|---:|---:|---:|---:|
| 1 | 42 | 6 | 6 | 12 |
| 2 | 42 | 6 | 6 | 12 |
| 3 | 41 | 6 | 6 | 12 |
| 4 | 41 | 6 | 6 | 12 |
| 5 | 41 | 6 | 6 | 12 |
| 6 | 34 | 6 | 6 | 12 |
| 7 | 34 | 6 | 6 | 12 |
| 8 | 34 | 6 | 6 | 12 |

The value `6` is the full dimension `p(5)-1` of the degree-five zero-mass
space.  The value `12` is therefore maximal: the parent and child row spaces
inside the coefficient space of prefix-indexed functions meet trivially.
By `(4)`, no fixed-pole projected transition exists on `E_2(M)`, for any
physical pole value `M`.

## 4. A seven-pair direct determinant

For a compact exact witness, take `M=1` and the seven parent/child pairs

```text
()      -> (1),
(1)     -> (1,1),
(2)     -> (1,2),
(3)     -> (1,3),
(4)     -> (1,4),
(5)     -> (1,5),
(1,1)   -> (1,1,1).                                        (7)
```

Divide both columns in each pair by the gcd of all their parent and child
coordinates.  The seven divisors are

```text
(1440,720,1440,720,1440,720,720).                           (8)
```

Use the six parent rows

```text
(5),(4,1),(3,2),(3,1,1),(2,2,1),(2,1,1,1)
```

and append the child `(5)` row.  The resulting normalized matrix is

```text
[  61987, 123972,  61955, 123488,  60963, 117724, 123970 ]
[ 264982, 512220, 246780, 471920, 222210, 407844, 494480 ]
[ 212848, 422652, 205980, 391712, 180480, 320196, 419612 ]
[ 497718, 917412, 419790, 759408, 338238, 591756, 842154 ]
[ 431766, 816756, 375090, 670584, 292542, 500412, 770878 ]
[-500833,-986466,-478553,-908966,-420417,-753946,-970113 ]
[  61986, 123970,  61954, 123486,  60962, 117722, 123968 ]. (9)
```

Its exact determinant is

```text
-3821168719915360450560
 =-2^11 3^6 5 7 11 37 877 204869207 !=0.                  (10)
```

The first six rows span all parent coordinate rows because the parent
currents have total mass zero.  The last row is not in their span by `(10)`.
This directly contradicts the row-space form of `(4)` and certifies the live
depth-two wall without relying on the larger rank census.

The last row differs from the first by the small vector

```text
(-1,-2,-1,-2,-1,-2,-2),                                   (11)
```

so the lost depth datum remains visible after the large common arithmetic
scales are removed.

## 5. Exact verification

Run

```text
python 04-computation/gmc_live_bank_projected_pole_transition_thm3168.py
python -O 04-computation/gmc_live_bank_projected_pole_transition_thm3168.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_live_bank_projected_pole_transition_thm3168.out.
```

The companion reconstructs all live currents from the upstream exact signed
bank, checks physical multiplicities and total-mass cancellation, verifies
the shallow rational kernels and high-degree full-column minors, proves all
eight maximal depth-two rank jumps, and reconstructs `(8)--(11)` with exact
integer arithmetic.

## 6. Scope

The shallow maps in `(6)` are maps on finite projected spans.  They are not
asserted to be canonical, compatible across degrees or pole values,
commuting around prefix squares, positive on the Hasse cone, or induced by
the full Pluecker action of THM-3160.  For `N>=7` they are merely the unique
linear interpolation of independent displayed parent columns.

The depth-two wall excludes linear current-only response transport on this
live bank.  It does not exclude a nonlinear rule, nor a rule supplied with
cross-degree endpoint minors.  THM-3163 shows that an arbitrary terminal law
has an abstract posterior Markov realization; `(10)` concerns response
compatibility and does not contradict that theorem.  No selector positivity,
original-response decomposition, NC2, or Gaussian Moment conclusion follows.

QED (candidate pending independent theorem audit).
