---
id: THM-2954
title: "Opposite endpoint-link blocks and transitive-cycle obstruction"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  The opposite strict-danger endpoint relation
  14*gcd(u,v)|(u+v) is exactly equality of the 2- and 7-adic
  valuations together with opposition of the residual units modulo
  14.  Its graph is therefore a disjoint union of complete bipartite
  blocks.  A nonzero transitively C_7-equivariant supported matrix
  would contain an odd 7-cycle and is impossible.  This arithmetic
  support theorem does not govern a strict-open cover: an a.e.
  endpoint switch leaves a safe seam, while a pointwise cover of
  multiplicity one a.e. has no switches and is component-monochromatic.
  No LRC(14) exclusion is claimed.
source: root-lrc-endpoint-link-bipartite-2026-07-29
audit: Pending independent hostile audit.
depends_on: []
related:
  - THM-767-r7-collision-network-kcl-and-scale-free-witness
  - THM-2953-cyclic-charged-corank-parity-and-reynolds-averaging-boundary
  - MISTAKE-146
script: 04-computation/lrc14_opposite_endpoint_link_blocks_thm2954.py
output: 05-knowledge/results/lrc14_opposite_endpoint_link_blocks_thm2954.out
script_sha256: 5bbdeff06bd89c7076c2b9997d8381c4f12402db34de00fda0409e78e986d4f7
output_sha256: 5d758888d796a18738577d9a8df57bcd676ecc9bfba3761f4ed635ba7ae03375
hash_basis: LF-normalized bytes
---

# THM-2954 -- opposite endpoint-link blocks and transitive-cycle obstruction

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Exact opposite-endpoint relation

For positive integer speeds `u,v`, an exit endpoint of the strict
danger comb `D_u` can coincide with an entry endpoint of `D_v` exactly
when

```text
14 gcd(u,v) divides u+v.                               (1)
```

This is the arithmetic coincidence law retained by THM-767 after the
open-seam correction.  Write

```text
u=2^alpha 7^beta U,       v=2^gamma 7^delta V,         (2)
gcd(U,14)=gcd(V,14)=1.
```

Then `(1)` is equivalent to

```text
alpha=gamma,
beta=delta,
U+V=0 mod14.                                          (3)
```

Indeed, `(1)` first forces equal `2`- and `7`-adic valuations.  Under
that equality put `d=gcd(U,V)`.  Since `d` is a unit modulo `14`,

```text
14 gcd(u,v)|(u+v)
 iff 14d|(U+V)
 iff 14|(U+V).                                        (4)
```

The last equivalence uses both `d|(U+V)` and `gcd(d,14)=1`.

Consequently the link graph is the disjoint union, over

```text
(alpha,beta,{r,-r}),
r in (Z/14Z)^*/{+-1},                                 (5)
```

of complete bipartite graphs.  One side consists of speeds whose
unit part is `r mod14`; the other consists of speeds whose unit part
is `-r mod14`.  The three unordered residue pairs are

```text
{1,13},                     {3,11},                    {5,9}. (6)
```

There are no self-links and no odd cycles.

## 2. No nonzero transitive C7 link support

Let seven labeled speeds be indexed by `F_7`.  A support relation is
transitively `C_7`-equivariant if membership of `(i,j)` depends only on
`j-i`.  Suppose every supported off-diagonal pair is an arithmetic
link from `(1)`.

If displacement `d!=0` occurs, equivariance supplies

```text
(i,i+d),                 i in F_7.                     (7)
```

These seven edges form an odd `7`-cycle because `d` generates
`F_7`.  This contradicts the bipartiteness in Section 1.  Since
diagonal links do not exist, the supported matrix is zero:

```text
transitive C_7 support + opposite-endpoint links
                   ==> empty support.                 (8)
```

Thus a nonzero endpoint/link response cannot simultaneously retain the
literal arithmetic support and become a transitive cyclic response.
This explains why THM-2953's stable-kernel gate cannot be obtained by
simply declaring an endpoint matrix equivariant.

## 3. The strict-open seam boundary

Equation `(8)` is a theorem about an abstract support relation, not
about persistence of an LRC cover.  The distinction is load-bearing.

For the strict danger sets, take

```text
u=1,                         v=13,                     x_0=1/14. (9)
```

The two combs cover a punctured neighborhood of `x_0`, one on each
side, but

```text
x_0 notin D_1 union D_13.                              (10)
```

Thus an a.e. partition can switch labels at a linked seam while
leaving one literal safe point.  Such a partition is not a
counterexample to LRC(14).

Conversely let `R` be an open carrier and let

```text
E_i=R intersect D_(z_i)                                (11)
```

be finitely many open subsets which

```text
cover R pointwise,
sum_i 1_(E_i)=1 a.e. on R.                             (12)
```

For `i!=j`, the intersection `E_i intersect E_j` is open.  By `(12)`
it is null, hence empty.  Therefore the `E_i` form an actually
disjoint open cover of `R`.  Every connected component of `R` is
contained wholly in one `E_i`:

```text
pointwise strict-open cover + multiplicity one a.e.
                ==> component-monochromatic.           (13)
```

There are no internal endpoint switches at all.  If seven labels also
have equal mass `mu(R)/7`, the remaining condition is exactly a finite
assignment of whole component lengths into seven equal bins.  Link
blocks neither strengthen nor weaken that allocation.

This is the repaired critical-wall split:

```text
bare a.e. wall: linked switches allowed, but seams are safe;
actual pointwise cover: no switches, whole-component assignment.    (14)
```

The endpoint-capacity inequalities corrected in MISTAKE-146 must not
be recovered from the first line of `(14)`.

## 4. Exact companion

Run

```text
python 04-computation/lrc14_opposite_endpoint_link_blocks_thm2954.py
python -O 04-computation/lrc14_opposite_endpoint_link_blocks_thm2954.py
```

The companion uses integer and rational arithmetic with explicit
exception gates.  It checks:

1. equivalence `(1)<->(3)` on every ordered pair `1<=u,v<=200`;
2. the three complete bipartite residue blocks inside every valuation
   class occurring in that box;
3. absence of diagonal links and odd cycles;
4. all six nonzero `C_7` displacement orbits in `(7)`;
5. the exact seam hostile `(9)--(10)` on two rational side probes; and
6. positive controls in each residue-opposition block.

No finite scan is used in the proof of `(3)`, `(8)`, or `(13)`.
Promotion requires an independent audit of both divisibility
directions, the cyclic-support quantifier, and the pointwise-versus-a.e.
scope.
