# Independent referee report on THM-4009

**Final verdict: PASS.** The theorem's main implication passed on first audit.
The referee found one finite-atlas label/count defect; the primary theorem,
script, output, and navigation digest subsequently repaired it exactly as
specified below.

The audited working-tree artifacts were:

```text
THM-4009 theorem SHA-256:
  e4cb215eec944e146edcd36c4158abef299dee6db4e7057eb8faf6daf975e21d
primary script SHA-256:
  c1ebbb34e69102cac0af8c0b1e0b6fc6d97b86f99cb31da2b96aaaac7c3a608e
primary output SHA-256:
  01e3b6ce4ab22c5793dd1dfd2416bb15629def7c9309d6bcaf0c3f9558f942c4
```

The independent implementation is not an import or wrapper of the primary
script.  It uses a multiplicity generating function for the primitive
histogram/sign-orbit atlas, polynomial convolution for the labelled lattice
ball, a Moebius sieve for primitive vectors, and the balanced-squares formula
for every `l1` optimization.

```text
independent script:
  04-computation/lrc14_euclidean_covering_transference_independent_referee_thm4009.py
independent script SHA-256:
  3830131aaed502bde29e933d49c9c740c53718900ab55ae8bf736cac7c743511
independent frozen output:
  05-knowledge/results/lrc14_euclidean_covering_transference_independent_referee_thm4009.out
independent output SHA-256:
  0967cff265a8de3c8bafd3745eac946c6ce0fdfca42bcef9de890db907c7da77
```

Normal and `-O` runs match the frozen output byte-for-byte.

## 1. Geometry and the exact inradius

Let `K=pi([-1,1]^13)` in `V=n^perp`.  For `u in V`, orthogonality of the
projection gives

```text
h_K(u)=sum_i |u_i|.
```

Thus the centred Euclidean inradius is the minimum of `||u||_1` over unit
vectors in `n^perp`.  On a fixed sign orthant, normalize by `||u||_1=1`.
The feasible set is a nonnegative polytope with the two affine equations
`sum |u_i|=1` and `n dot u=0`.  A convex maximum of `||u||_2` occurs at an
extreme point, which has at most two positive absolute-coordinate variables.
One variable cannot balance a positive speed vector, so exactly two remain,
with opposite signs.  The pair `(i,j)` gives

```text
||u||_1/||u||_2=(n_i+n_j)/sqrt(n_i^2+n_j^2).
```

This quantity decreases with the ratio of the larger entry to the smaller.
The exact minimum is consequently attained at `(m,M)`, proving

```text
rho(n)=(m+M)/sqrt(m^2+M^2)>1.
```

The universal inclusion `B_V(0,1) subset K` is also immediate: a unit vector
in `V` has all coordinates at most one in absolute value and projects to
itself.  No quotient norm or projection factor is missing.  The independent
companion checks the extreme-pair minimizer on all `792` five-subsets of
`{1,...,12}`.

## 2. Closed boundary and transference normalization

A counterexample makes the **closed** zonotope disjoint from the projected
lattice.  Since a full Euclidean lattice is discrete, a nearest point to `c`
exists.  Therefore disjointness of the closed radius-`3/7` ball gives

```text
dist(c,Lambda)>3/7,
mu(Lambda)>3/7.
```

The strict sign is valid and load-bearing.  The one-dimensional hostile is
the open radius-`1/2` ball about `1/2` modulo `Z`: it is disjoint while the
distance is exactly `1/2`; the closed ball is not disjoint.

The cited Euclidean Banaszczyk inequality has the normalization

```text
mu(L) lambda_1(L*) <= d/2,
```

equivalently `2 mu(L) lambda_1(L*)<=d`.  Applying it to the rank-twelve
lattice `Lambda` has the correct primal/dual orientation.  Together with the
strict radius it gives

```text
lambda_1(Lambda*) < 6/(3/7)=14.
```

The bibliographic citation and DOI are correct.  The theorem would be even
easier to audit later if it pinned the result as the covering-radius
transference theorem inside the paper, but this is a citation-quality
suggestion, not a mathematical gap.

The quotient-dual identity also passes directly:
`<a,pi(z)>=<a,z>` for `a in V`, so integral pairing with all projected standard
basis vectors is exactly `a in Z^13`; imposing `a in V` is `a dot n=0`.

## 3. Integer caps and Graver scope

The strict Euclidean bound gives integer square norm at most `195` and
coordinate height at most `13`.  For thirteen nonnegative absolute
coordinates of total `L`, the minimum square mass is obtained by balancing
them.  It equals `194` at `L=50`, with sorted vector

```text
(3,3,4,4,4,4,4,4,4,4,4,4,4),
```

and equals `201` at `L=51`.  Hence the exact cap is `l1<=50`.  Parity of
square mass and `l1` also shows that no unbalanced second extremizer can fit
between square masses `194` and `195`; the sorted extremizer is unique.

Choosing a shortest nonzero Euclidean vector is sufficient for the Graver
claim.  In any nontrivial conformal integer decomposition, either summand is
coordinatewise dominated with a strict inequality somewhere, hence has
strictly smaller Euclidean norm.  There is no hidden assumption that the
transference vector itself already came Graver-minimal.

Every bounded-spread entry independently reproduces:

```text
(R, square cap, l1 cap)
(2,108,37), (3,122,39), (5,141,42),
(13,169,46), (21,178,47), (50,188,49).
```

The strict threshold is handled correctly when it is integral, notably
`<170` at `R=13`, hence `<=169`.

## 4. The support-two count repair

The mathematical ratio count is correct:

```text
#{(p,q): 1<=p<q, gcd(p,q)=1, p^2+q^2<=195}=47.
```

The maximum `p+q` is `19`, and the Pell packet leaves only `(1,5)`.  The
independent ordered ratio-list digest is

```text
78af8034cb48996c028cdc0fb5259542a3f042220d80e9022819a48deadbcda4.
```

However,

```text
47*C(13,2)=3,666
```

counts **unoriented ratio/support packets**.  For an arbitrary labelled speed
vector as stated in the theorem, each support has two assignments of `p<q`
to its two coordinate labels.  Modulo global sign there are therefore

```text
2*47*C(13,2)=7,332
```

oriented labelled coefficient assignments.  The assignment becomes
determined only after the two actual speeds are compared, or if a global
increasing speed order is explicitly imposed.  This is the same distinction
already kept in THM-3818 between unoriented supports and oriented
assignments.

Required repair, subsequently applied:

- change the status line's “3,666 labelled placements” to “3,666 unoriented
  ratio/support packets and 7,332 oriented labelled assignments”;
- make the same wording change in Section 3 and in the primary script/output;
- do not change the correct count `47` or any main theorem implication.

This is a finite-atlas bookkeeping defect, not a failure of the short-relation
theorem. The corrected final artefacts consistently use `3,666` for
unoriented ratio/support packets and `7,332` for oriented labelled assignments
modulo global sign.

## 5. Higher-support atlas and labelled ball

The independent multiplicity generating function exactly reproduces the
support histogram row

```text
2:47, 3:209, 4:566, 5:1177, 6:2057, 7:3180, 8:4490,
9:5911, 10:7374, 11:8805, 12:10188, 13:11455.
```

It also reproduces:

```text
55,459 primitive absolute histograms,
5,030,161 nonconstant sign-multiset orbits,
28,315 odd-l1 absolute histograms,
711,202,814,025,242 nonzero labelled ball vectors,
711,119,925,281,794 primitive labelled ball vectors.
```

For sign orbits the audit uses Burnside directly.  If the distinct absolute
values have multiplicities `m_v`, the raw sign-multiset count is
`prod_v(m_v+1)`; sign complement fixes a multiset exactly when every `m_v` is
even; the all-positive/all-negative orbit is then removed.  The special
support-two histogram `(1,1)` and its one mixed-sign orbit are correctly
removed because distinct positive speeds cannot support that relation.

## 6. Half-lattice character and rank-eleven join

The parity sidecar is correctly separated from the short row.  Since
`2c=pi(1) in Lambda` but `c notin Lambda`, bidual separation gives a dual
vector pairing by `1/2 mod Z`; the pairing is `(1/2)sum b_i`, so some relation
has odd coefficient sum.  Nothing in ordinary covering transference says
that a shortest vector lies in this character coset.

There is also an elementary hostile check.  A primitive speed row is either
all odd, in which case `t=1/2` makes every runner lonely and every kernel
relation has even coefficient sum, or it contains an even and an odd speed.
The primitive support-two relation on such a pair has odd coefficient sum.
The independent script checks all `8,191` nonzero mod-two speed patterns and
four integer pair controls.

The rank-eleven determinant cap also passes.  Eleven support-at-most-three
height-`Q` rows have norm at most `sqrt(3)Q`; the outside short row has norm at
most `sqrt(195)`.  Hadamard on the twelve maximal-minor rows gives

```text
floor(sqrt(195*3^11)*Q^11),  Q=91^6,
```

equal to the printed 134-digit integer.  Division of the cofactor vector by
its common gcd can only decrease coordinates, so primitiveness of the speed
row licenses the displayed speed cap.  A short row of support at most three
has height at most `13<Q`, hence is an actual generator of `W`; only support
at least four can be the outside-span twelfth row.

## 7. Final scope

After the support-placement wording repair, the result deserves
`CITED + PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED`:

```text
hypothetical LRC(14) counterexample
  => a nonzero Graver speed relation with
     ||a||_2<14, sum a_i^2<=195, ||a||_1<=50, |a_i|<=13.
```

The AP13 hostile confirms the stopping boundary: it is lonely at `t=1/14`
while carrying the square-norm-six Graver relation `(1,-2,1,0,...)`.
Short resonance remains a necessary atlas reduction, not a sufficient
loneliness obstruction and not a proof of LRC(14).
