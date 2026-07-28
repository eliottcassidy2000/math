---
id: THM-2783
title: "Weighted long-wall binary null avoidance and ternary state reconstruction"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For an arbitrary integral wall h and every independent
  (k-1)-row B_k root frame, the determinant is
  +/-2^c<epsilon,h>, where c counts unbalanced signed-cycle components and
  epsilon is the unique deficient-tree sign vector.  Hence h is transverse
  to every such frame iff its Boolean subset sums are distinct.  Among
  positive integral walls, the unique minimum-L1 universal null-avoider is,
  up to permutation, (1,2,...,2^(k-1)).  If the stronger requirement is to
  reconstruct every {-1,0,1} tree state, the unique minimum-L1 wall is
  (1,3,...,3^(k-1)).  Thus binary is the optimal zero detector and ternary
  the optimal full signed-state address.  These are coding theorems, not a
  PSL2(Z), graceful-tree, Keller, or LRC action.
source: root/weighted-long-wall-coding-2026-07-28
depends_on: []
related:
  - THM-2774-tree-path-smith-index-ladder-and-binary-ternary-lattice-defects
  - THM-2776-one-long-wall-signed-component-imbalance-determinant-classification
script: 04-computation/weighted_long_wall_binary_ternary_thm2783.py
output: 05-knowledge/results/weighted_long_wall_binary_ternary_thm2783.out
script_sha256: e6280567d7e895a54061fb1610ab9493f7509eb3a6e65a45dbd041194cf9bd6c
output_sha256: 9543c166f90d9ab3a7da6b09e1d2034f8e72f952fc53a24d65cd1afe1ceb684c
hash_basis: LF-normalized bytes
---

# THM-2783 -- binary avoids the wall; ternary remembers the state

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

There are two different reasons that the primes two and three appear in the
signed-tree arrangement.  They should not be conflated:

```text
binary:   cheapest positive integral wall avoiding every determinant zero;
ternary:  cheapest positive integral wall naming every absent/+/- tree state.
                                                                    (1)
```

The distinction is exact.  It comes from the two alphabets `{0,1}` and
`{-1,0,1}`, not from silently identifying them with the torsion generators
of the modular group.

## 1. Arbitrary-wall determinant formula

Fix `k>=1`.  Up to row sign, a type-`B_k` root is

```text
e_i,                    e_i-e_j,                    e_i+e_j.          (2)
```

Let `N` be a `(k-1) x k` matrix of roots `(2)` with rank `k-1`, and let

```text
h=(h_1,...,h_k) in Z^k.                                  (3)
```

Represent `e_i` by a half-edge and `e_i-s_ij e_j`, `s_ij in {+/-1}`, by a
signed ordinary edge.  The component count from the signed-graph proof is
forced: there is one deficient unanchored tree `U`; every other component is
an anchored tree or an unbalanced signed cycle.  Let `c` count the latter.
On `U`, choose either of the two sign vectors

```text
epsilon in {+/-1}^U,       epsilon_i=s_ij epsilon_j.      (4)
```

Then

```text
det [N;h] = +/- 2^c sum_(i in U) epsilon_i h_i.           (5)
```

Here is a self-contained proof.  In a connected component with `v` vertices,
`q` ordinary edges, and `a` half-edges, row independence gives
`q+a<=v`, connectivity gives `q>=v-1`, and the global row deficit is one.
Thus precisely one component has `q=v-1,a=0`; it is the tree `U`.  Every
square component is either `q=v-1,a=1`, with determinant `+/-1` by leaf
elimination, or `q=v,a=0`, with determinant `1-product(s)=+/-2`; the balanced
cycle is singular.  The maximal cofactors of the tree block are, up to a
common sign, the entries of `(4)`.  Expanding the appended row `h_U` gives
its factor `+/-<epsilon,h_U>`, and the square blocks give `2^c`.  This proves
`(5)`.

The all-ones case is the imbalance theorem isolated in candidate THM-2776.
Equation `(5)` is stronger in a different direction: the wall itself is now
allowed to carry arithmetic state.

## 2. Universal regularity is exactly Boolean dissociation

Call `h` **universally regular** if

```text
det[N;h]!=0
```

for every rank-`k-1` root frame `N` above.  The following are equivalent:

```text
(i)   h is universally regular;
(ii)  delta dot h !=0 for every nonzero delta in {-1,0,1}^k;
(iii) the 2^k subset sums sum_(i in A) h_i are pairwise distinct.    (6)
```

The equivalence `(ii)<->(iii)` is the disjoint-support decomposition of a
difference of two subset sums.  Equation `(5)` proves `(ii)->(i)`.

For the converse, fix a nonzero `delta` and put `U=supp(delta)`.  Join the
vertices of `U` by any tree and give an edge `ij` the row

```text
e_i-(delta_i delta_j)e_j.                                (7)
```

Its tree kernel is `delta|_U`.  Put a half-edge `e_j` on every vertex outside
`U`.  These are `k-1` independent rows, there are no cycle components, and

```text
det[N_delta;h]=+/- delta dot h.                           (8)
```

Thus every signed-sum zero is witnessed by an actual lawful root frame, and
`(i)->(ii)`.  This proves `(6)` with no genericity or asymptotic qualifier.

The simple wall

```text
h=(1,2,3)
```

already fails: `delta=(1,1,-1)` gives zero.  One explicit witness has rows
`e_1-e_2` and `e_2+e_3`.

## 3. Binary powers are the unique cheapest null avoider

Now suppose all `h_i` are positive integers and put

```text
H=sum_i h_i.                                               (9)
```

If `h` is universally regular, its `2^k` subset sums are distinct integers
in `[0,H]`.  Therefore

```text
H>=2^k-1.                                                 (10)
```

Equality is attained by

```text
h_bin=(1,2,4,...,2^(k-1)).                               (11)
```

Moreover `(11)` is the unique equality case up to coordinate permutation.
Indeed, equality in `(10)` forces the subset sums to be every integer from
zero through `2^k-1`, each once, so

```text
product_i (1+x^(h_i))=1+x+...+x^(2^k-1).                (12)
```

After sorting, the coefficient of `x` forces `h_1=1`.  Dividing `(12)` by
`1+x` gives

```text
product_(i>1)(1+x^(h_i))=1+x^2+...+x^(2^k-2).           (13)
```

Every coefficient is nonnegative, so `(13)` forces every remaining `h_i`
to be even.  Substitute `y=x^2` and induct.  This gives exactly
`1,2,...,2^(k-1)`.

Thus binary is not merely a convenient hostile-proof encoding.  It is the
unique positive integral wall of minimum total height that makes every
rank-`k-1` signed-root frame transverse.

## 4. Balanced ternary is the unique cheapest full-state address

Universal regularity remembers only whether the signed tree sum is zero.
To reconstruct the entire deficient-tree state, require the map

```text
{-1,0,1}^k -> Z,               delta |-> delta dot h      (14)
```

to be injective.  There are `3^k` values in `[-H,H]`, so

```text
H>=(3^k-1)/2.                                             (15)
```

Equality is attained by the balanced-ternary wall

```text
h_ter=(1,3,9,...,3^(k-1)).                               (16)
```

Again the equality case is unique up to permutation.  Shift every balanced
digit by one.  Equality in `(15)` says that all exponents from zero through
`3^k-1` occur exactly once:

```text
product_i (1+x^(h_i)+x^(2h_i))=1+x+...+x^(3^k-1).       (17)
```

The coefficient of `x` forces `h_1=1`.  Dividing by `1+x+x^2` gives

```text
product_(i>1)(1+x^(h_i)+x^(2h_i))
  =1+x^3+x^6+...+x^(3^k-3).                             (18)
```

Nonnegative coefficients force every remaining exponent `h_i` to be
divisible by three.  Put `y=x^3` and induct.  This proves `(16)`.

The stronger requirement is genuinely stronger.  The binary wall `(1,2)`
is universally regular, but the two signed states

```text
(1,0) and (-1,1)
```

both have address one.  Binary detects that neither is zero; ternary
distinguishes them.

There is one more necessary sidecar.  Equation `(5)` contains both the tree
address and the cycle count `c`.  Even balanced ternary does not promise to
recover `c` from the single product `2^c<epsilon,h>`; a complete frame state
must retain `c` separately.

## 5. What this does—and does not—say about two and three

The exact state alphabets are

```text
{absent,present}                  of size 2;
{absent,positive,negative}        of size 3.              (19)
```

Equations `(10)--(18)` show that radices two and three are the optimal
interval packings for those two tasks.  This gives a rigorous second answer,
complementary to THM-2774's first `Z/2` and `Z/3` lattice defects, for why
the same primes recur.

It is not a modular-group action.  In

```text
PSL_2(Z)=C_2*C_3
```

the numbers two and three are orders of torsion operations and reduced-word
history matters.  Here they are cardinalities of coding alphabets and the
maps are linear evaluations.  No action of either free factor on a wall,
tree frame, Keller field, or LRC packet has been constructed.

The same scope boundary matters elsewhere.

- Replacing the geometric all-ones long wall by `(11)` or `(16)` deforms the
  graceful arrangement.  It regularizes signed frames but does not select a
  graceful labeling or prevent cancellation in the original graceful
  coefficient.
- The weighted wall is a useful lexicographic owner/state sidecar, but its
  exponential height supplies no physical LRC relation, endpoint phase, or
  terminal current.
- The `D3`/quartic wall in THM-2775 is Weyl-geometric.  A radix wall breaks
  that symmetry and is not an affine operation on a Keller map.

## 6. Exact verification

Run

```bash
python 04-computation/weighted_long_wall_binary_ternary_thm2783.py
python -O 04-computation/weighted_long_wall_binary_ternary_thm2783.py
```

The dependency-free companion uses exact integer arithmetic, explicit
exceptions, and no truth-bearing Python assertions.  It checks `(5)` for
`43,880` weighted full-rank frames through `k=5`; constructs all `1,086`
nonzero signed-state witness frames through `k=6`; verifies that binary and
ternary walls never vanish there; finds the expected `145` linear-wall
zeros; exhausts the minimum-sum binary equality cases through `k=6` and the
ternary equality cases through `k=4`; and checks the radix polynomial
identities and sharp hostiles.  Normal and optimized runs byte-match the
stored transcript.

```text
PROVED HERE (candidate):  arbitrary-wall determinant formula (5);
                          universal regularity iff Boolean dissociation;
                          an explicit frame for every signed-sum zero;
                          sharp binary minimum and unique equality wall;
                          sharp ternary minimum and unique equality wall;
                          binary null-detection / ternary reconstruction fork;
                          cycle-count and modular-action stopping boundaries.

NOT PROVED:               a canonical radix deformation of the graceful wall;
                          graceful-label existence or Graceful Tree;
                          a PSL2(Z) action or free-factor realization;
                          a Keller, JC(2), or DC(2) consequence;
                          an LRC relation, endpoint, or current;
                          LRC(14).                                      (20)
```

QED (candidate).
