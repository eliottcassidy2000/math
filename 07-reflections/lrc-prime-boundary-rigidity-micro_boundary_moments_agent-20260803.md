# Prime right-boundary interpolation proves all-prime micro-staircase rigidity

**Status: PROVED ELEMENTARY THEOREM + FINITE-EXACT CONTROLS; PROMOTED AS
[THM-3316](../01-canon/theorems/THM-3316-prime-right-boundary-interpolation-forces-scalar-rigidity.md).**
At every prime modulus, the cells immediately
to the right of `alpha=a/p` already force every full micro-staircase blocker to
be a scalar ramp.  In particular, the proposed `p=13` support-seven moment
obstruction is correct, support eight is also impossible, and in fact every
normalized support from one through eleven is impossible.

This closes the prime-modulus instances of the finite residue/cell rigidity
problem.  It is not a speed-to-residue lift, not LRC(13), and not LRC(14).
The composite `n=14` right-boundary hostile survives unchanged.

The matching exact companion is
[`lrc_prime_boundary_rigidity_micro_boundary_moments_agent_20260803.py`](../04-computation/lrc_prime_boundary_rigidity_micro_boundary_moments_agent_20260803.py),
with frozen
[`output`](../05-knowledge/results/lrc_prime_boundary_rigidity_micro_boundary_moments_agent_20260803.out).

## Inheritance pass and concept board

- **Closest proved mechanism:** [THM-363](../01-canon/theorems/THM-363-lrc-scalar-gauge-reindexing.md)
  identifies scalar addition with a shift of the cell parameter, while
  [THM-364](../01-canon/theorems/THM-364-lrc-scalar-ramp-cell-blocking.md)
  proves that every scalar ramp really is a full blocker.
- **Canonical hostile:** at `n=14`, the normalized vector
  `(0,1,0,...,0)` covers all `14^2` right-adjacent affine tests while missing
  interior micro-staircase candidates.  Boundary sufficiency is false at the
  target composite modulus.
- **Corrected near miss:** the previous `n=13` result was FINITE-EXACT by
  prime-unit orbit splitting and eleven UNSAT branches.  This note replaces
  the converse computation by a uniform proof; it does not retroactively turn
  solver evidence into proof.
- **Least-used sidecar:** for a fixed support coordinate `j`, retain the full
  multiset of the `2|S|` bad shifts, not merely its union.  Varying `j` makes
  its vanishing polynomial interpolate.

The live concepts were scalar gauge, right-boundary affine covers, support
interpolation, a moving root multiset, and the finite-field zero-sum
functional.  The last three collapse the whole prime case.

## Exact right-boundary convention

Let `p` be prime and index the coordinates by `1<=i<=p-1`.  Immediately to
the right of `alpha=a/p`, with the perturbation small enough for every
coordinate, the cell pattern is

```text
b_i = a i mod p.                                             (1)
```

Thus a residue vector `v` passes all right-boundary tests only if for every
`(a,s) in F_p^2` some coordinate satisfies

```text
a i+s v_i in {0,-1}.                                        (2)
```

Subtract the scalar ramp `v_1 i` and put

```text
w_i=v_i-v_1 i,        so w_1=0.                             (3)
```

Equation (2) is preserved: the expression for `w` at `(a,s)` is the
expression for `v` at `(a-sv_1,s)`.  This also gives a boundary-only proof of
the relevant part of the THM-363 gauge.

## The interpolation theorem

**Theorem.**  Let `p` be prime.  If `v in F_p^(p-1)` covers every
right-boundary test (2), then

```text
v_i=m i  for one m in F_p and every i.                       (4)
```

Consequently the full micro-staircase blockers at prime modulus are exactly
the `p` scalar ramps.

### Proof

Normalize as in (3).  Suppose `w` is nonzero and write

```text
S={l:w_l!=0},    h=|S|,    q_l=w_l^(-1).                    (5)
```

Because `w_1=0`, one has `1<=h<=p-2`.

Fix `j in S` and take `a=-j^(-1)` in (2).  A zero coordinate cannot block:
`-i/j` is neither `0` nor `-1` unless `i=j`, and `j` is in the support.  A
supported coordinate `l` blocks precisely when

```text
s = l q_l/j       or       s=(l-j)q_l/j.                   (6)
```

After scaling the shift by `y=js`, full boundary coverage says that, for
every `j in S`, the multiset

```text
{l q_l,(l-j)q_l : l in S}                                  (7)
```

covers all of `F_p`.

Fix `y in F_p` and define

```text
U(y)   = product_(l in S) (y-l q_l),
G_y(J) = product_(l in S) (y-(l-J)q_l),
Q      = product_(l in S) q_l,
Pi(J)  = product_(l in S) (J-l).                            (8)
```

By (7), `U(y)G_y(j)=0` at all `h` values `j in S`.  Let

```text
X={l q_l:l in S}.                                           (9)
```

If `y` is not in `X`, then `U(y)` is nonzero.  The degree-`h` polynomial
`G_y(J)` vanishes at every member of `S` and has leading coefficient `Q`.
Therefore

```text
G_y(J)=Q Pi(J).                                             (10)
```

The roots on the two sides give the multiset identity

```text
{l-y w_l:l in S}=S             for every y notin X.         (11)
```

Now set

```text
A(Y)=product_(l in S)(l-Yw_l),      L=product_(l in S)l.    (12)
```

Every label in `S` is nonzero, so `L!=0`.  Equation (11) says `A(y)=L` off
`X`, while the definition of `X` says `A(y)=0` on `X`.  If `r=|X|`, then

```text
sum_(y in F_p) A(y)=(p-r)L=-rL != 0,                        (13)
```

because `1<=r<=h<p`.

On the other hand, `deg A=h<=p-2`.  The sum over `F_p` of every monomial of
degree at most `p-2`, including the constant monomial, is zero.  Hence the
left side of (13) is zero, a contradiction.  Thus `w=0`, proving (4).
THM-364 supplies the converse for the full open-cell system.  QED

The proof is constructive.  For a normalized nonscalar vector, some
`y notin X` has `A(y)!=L`; otherwise (13) would hold.  Then (10) fails for
some `j in S`, and

```text
a=-j^(-1),       s=yj^(-1)                                 (14)
```

is an explicit missed right-boundary candidate.

## Finite-geometry interpretation and field tangent

For each coordinate `i`, equation (2) is the union of two parallel affine
lines in the dual `(a,s)`-plane,

```text
a i+s v_i=0,             a i+s v_i=-1.                    (15)
```

The theorem is therefore a structured affine-line-cover rigidity statement:
if the `2(p-1)` labelled lines cover `F_p^2`, all normals `(i,v_i)` have one
projective direction.  Indeed, for a scalar ramp `v_i=mi`, the first family
covers `a+ms=0`, and otherwise the unique label
`i=-(a+ms)^(-1)` supplies the second line.  The polynomial in (8) is a
Rédei-type direction polynomial; its unusual feature is that the coordinate
labels provide exactly the interpolation nodes needed to recover it.

The same calculation over a nonprime finite field `F_q` stops at a precise
torsion sidecar rather than giving an automatic contradiction.  If a
normalized cover with support `h` existed and `X` were defined by (9), then
for every `0<=k<=q-2-h`, summing `Y^k A(Y)` would force

```text
sum_(x in X) x^k=0.                                        (16)
```

In particular `|X|=0` in the ground characteristic.  Over `F_p` this is
impossible because `1<=|X|<p`; over a proper extension field it permits
characteristic-sized collision packets.  This cleanly identifies the
zero-sum step, rather than interpolation itself, as the prime-only closure.
It suggests that any extension-field or CRT analogue must retain the
collision packet `X` instead of quotienting it to its cardinality.

## Audit of the critical-support moment route

The originally proposed argument is sound, including its signs, but is now a
strictly weaker corollary.  For `h=(p+1)/2`, (7) has `p+1` incidences, so it is
`F_p` plus one repeated value.  Power moments in `J` of degrees
`1,...,h-1` and interpolation at the `h` support labels give, with

```text
B=sum_(l in S)q_l,
sum_(l in S) q_l^m = B^m       (1<=m<=h-1).                (17)
```

Newton identities then give the root polynomial.  If `B=0`, it has the form
`T^h-c`; it has at most

```text
gcd(h,p-1)<=2                                                  (18)
```

nonzero roots, all simple, and cannot contain the `h` values `q_l` with
multiplicity.

If `B!=0`, normalize `B=1`.  Fermat's relation
`sum q_l^(p-1)=h` forces the signed constant, giving

```text
F(T)=T^h-T^(h-1)+(h-1).                                    (19)
```

For a nonzero root, Euler's criterion reduces (19) to two candidates:

```text
T=h          on the quadratic-nonresidue branch,
T=2-h        on the quadratic-residue branch.              (20)
```

The derivative has only one nonzero critical point, and that root is simple
in the derivative.  Hence at most one candidate can be double and the total
`F_p` root multiplicity is at most three, less than `h` for `p>=7`.  At
`p=13`, this recovers the stated support-seven contradiction.  The main
interpolation theorem also disposes of support eight and every other
normalized support without a critical-cardinality split.

## Failure boundary and connection ledger

The field hypotheses enter in three precise places: every nonzero `w_l` has
an inverse, every support label `j` has an inverse, and the sum of low-degree
monomials over `F_p` vanishes.  They cannot simply be deleted.  At `n=14`,

```text
w=(0,1,0,...,0)                                             (21)
```

covers all `196` right-boundary tests, so the analogue of the theorem is
false even at support one.

```text
source:      prime right-boundary affine-cover shadow;
target:      a degree-h moving-root polynomial A(Y);
map:         choose a=-j^(-1), scale y=js, interpolate in j;
preserved:   a missed boundary candidate, hence failure of full blocking;
destroyed:   cell width, prime-grid realization, and physical speed data;
sidecar:     support S and collision set X={l/w_l};
test:        compare A(y) with product(S) off X.
```

The prime/composite split is therefore structural rather than merely
computational.  A composite successor would need a CRT/valuation sidecar that
tracks nonunits before any interpolation; the `n=14` hostile shows that the
right-boundary shadow alone cannot finish that case.

## Exact controls and reproduction

The companion verifies (1) at `3<=p<=31`, exhausts every normalized vector
for `p=3,5,7`, extracts and checks (14) for every nonzero vector in those
universes, tests a deterministic dense bank through `p=19`, and verifies
`5,780,507` scalar-gauge candidate equivalences.  It independently audits
(18)--(20) through `p=101` and rechecks the composite hostile (21).

Reproduce from the repository root:

```powershell
python 04-computation/lrc_prime_boundary_rigidity_micro_boundary_moments_agent_20260803.py
python -O 04-computation/lrc_prime_boundary_rigidity_micro_boundary_moments_agent_20260803.py
```

Both modes emit the frozen transcript exactly.  Their LF-normalized SHA-256
digests (equal to the raw-file hashes in this worktree) are

```text
script: 774288510f30ff36cc0b2f863f0f9d3af809c0ff958d8fc035b311ff7dea5cd8
output: 2eeff12d06d8588f94641877d1900388bfe2605214228f6bbc1b5ce85a31b83b
```

The exact n=13 SAT classification remains a valuable independent audit, but
is no longer the strongest truth source for the prime converse.
