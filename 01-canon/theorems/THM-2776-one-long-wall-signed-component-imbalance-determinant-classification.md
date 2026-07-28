---
id: THM-2776
title: "One-long-wall signed-component imbalance determinant classification"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For k-1 independent B_k roots and the long wall h=(1,...,1),
  the signed support graph has one deficient unanchored tree U; every other
  component is a half-edge tree or an unbalanced signed cycle.  If c is the
  number of cycle components and epsilon is the +/-1 tree kernel, then the
  full determinant is +/-2^c sum_U epsilon.  This classifies the zero
  boundary and every attainable nonzero magnitude.  Odd primes are bounded
  by |U|<=k; a difference spanning tree attains k.  Binary cycle factors
  first create values beyond that spanning-tree range: 6 at k=5 and 8 at
  k=6.  The associated diagonal
  character is only total degree, so it cannot split homogeneous graceful
  coefficient cancellation.
source: >
  a4-resolvent-next-gate/one-long-wall-imbalance-2026-07-28;
  root/one-long-wall-classification-completion-2026-07-28
depends_on:
  - THM-2770-tree-incidence-a-d-weyl-clutch-and-four-vertex-fan-dichotomy
  - THM-2774-tree-path-smith-index-ladder-and-binary-ternary-lattice-defects
related:
  - THM-2777-marked-d3-six-root-determinant-tournament-and-binary-ternary-edge-spectrum
script: 04-computation/one_long_wall_signed_imbalance_thm2776.py
output: 05-knowledge/results/one_long_wall_signed_imbalance_thm2776.out
script_sha256: 0ee0950fda2a5f5c0da4f153fd6e663852545fecbec72d93de45f51419259ee6
output_sha256: 867e1e58d4082a15d57841b531a81c3aaea18d722aefb7a446c2801f88862ea7
hash_basis: LF-normalized bytes
---

# THM-2776 -- one long wall and the signed-component imbalance

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2774 isolates the canonical path frame and its `Z/k` quotient.  The
surrounding one-long-wall arithmetic admits a complete classification.
Binary cycle defects and long-wall imbalance factor independently:

```text
determinant = power of two * signed tree imbalance.       (1)
```

This explains both why the path gives exactly `k` and why larger composite
determinants can appear without introducing any new odd prime larger than
the path support.

## 1. Signed support graph

Fix `k>=2` and put

```text
h=(1,...,1) in Z^k.                                      (2)
```

A positive type-`B_k` hyperplane normal is, up to row sign,

```text
e_i,                       e_i-e_j,                e_i+e_j. (3)
```

Let `N` be a `(k-1) x k` matrix of such rows and assume

```text
rank_Z N=k-1.                                             (4)
```

Represent `e_i` by a half-edge at vertex `i`, and represent
`e_i-s e_j`, `s in {+/-1}`, by a signed ordinary edge.  Connected
components refer to the ordinary support graph; an isolated vertex is a
component, and half-edges do not join vertices.

For a component let

```text
v=#vertices,       q=#ordinary edges,       a=#half-edges,
r=q+a=#rows.                                             (5)
```

Because the ordinary graph is connected, `q>=v-1`.  Row independence gives
`r<=v`.  Summing over components and using `(4)` gives

```text
sum_components (v-r)=k-(k-1)=1.                          (6)
```

Hence there is exactly one component `U` with `v-r=1`; every other
component is square.

## 2. The three possible components

The inequalities in `(5)--(6)` leave only three row-independent types.

1. **Deficient tree.**  Here `r=v-1`, so `a=0`, `q=v-1`, and the ordinary
   graph is a tree.  This is the unique component `U`.
2. **Anchored tree.**  Here `r=v`, `a=1`, and `q=v-1`.  Leaf elimination
   against the half-edge gives determinant `+/-1`.
3. **Unbalanced cycle.**  Here `r=v`, `a=0`, and `q=v`.  The graph is
   unicyclic.  Switching signs along a spanning tree reduces its last row to
   `1-product_cycle(s)`.  Row independence forces the cycle unbalanced, so
   the determinant is `+/-2`.

A balanced unicyclic block has rank `v-1` and is excluded by `(4)`.  No
component with two half-edges or more than one cycle can occur in an
independent `(k-1)`-row selection.

Let `c` be the number of unbalanced-cycle components.  The determinant
contributed by all square components is therefore

```text
+/-2^c.                                                   (7)
```

## 3. The deficient-tree cofactor is its sign imbalance

Write a row of the tree `U` as `e_i-s_ij e_j`.  Since `U` is a tree there is,
up to global sign, a unique vector

```text
epsilon in {+/-1}^U,       epsilon_i=s_ij epsilon_j,     (8)
```

and `N_U epsilon=0`.  The maximal cofactors of the tree-incidence block are,
up to one common sign, the entries of `epsilon`.  Appending the restriction
`h_U=(1,...,1)` and expanding along its last row gives

```text
det [N_U;h_U]=+/- sum_(u in U) epsilon_u.                 (9)
```

After permuting components, every square block must use its own columns in
the full determinant expansion.  Combining `(7)` and `(9)` proves the exact
factorization

```text
det [N;h]=+/- 2^c sum_(u in U) epsilon_u.                (10)
```

In particular, under the rank hypothesis `(4)`,

```text
det[N;h]=0  iff  sum_U epsilon=0.                        (11)
```

This is the sharp zero boundary: `U` has even size and its two switched
sign classes have equal cardinality.  If `(4)` fails, the full determinant
is of course also zero, but for the earlier row-dependence reason.

## 4. Exact attainable magnitudes

Put `v=|U|` and `s=|sum_U epsilon|`.  Every nonzero determinant magnitude is

```text
2^c s,
1<=v<=k,        1<=s<=v,        s congruent_to v (mod 2),
0<=c<=floor((k-v)/2).                                  (12)
```

Conversely every value in `(12)` occurs.  Choose a tree on `v` vertices and
edge signs realizing any desired `+/-1` kernel with imbalance `s`; use `c`
disjoint two-vertex unbalanced cycles

```text
{e_i-e_j, e_i+e_j};                                     (13)
```

and fill the remaining vertices by singleton half-edges.  This realizes the
component data and hence `(12)`.

Two consequences are immediate.

- Every odd prime divisor of a nonzero determinant divides `s`, so

  ```text
  p_odd<=s<=v<=k.                                        (14)
  ```

- Taking `U` to be a spanning tree of ordinary differences gives
  `epsilon=(1,...,1)`, `c=0`, and determinant `k`.  THM-2774's path frame is
  the canonical geodesic realization.  It maximizes the imbalance, though
  not always the total determinant: binary cycle factors first produce
  values beyond the spanning-tree range, namely `6` at `k=5` and `8` at
  `k=6`.

The exact nonzero supports through six coordinates are

```text
k=2: {1,2};                 k=3: {1,2,3};
k=4: {1,2,3,4};             k=5: {1,2,3,4,5,6};
k=6: {1,2,3,4,5,6,8}.                                  (15)
```

Thus two is the primitive signed-cycle factor and three is the first odd
long-wall imbalance, but neither prime is exclusive at larger support.

## 5. The natural cyclic character does not repair coefficient signs

For the spanning difference-tree frame, THM-2774 identifies the quotient
by the coordinate-sum map

```text
sigma_k(z)=sum_i z_i mod k.                              (16)
```

This is valuable lattice data, but it is maximally coarse on a homogeneous
polynomial: every exponent vector of total degree `D` has the same image
`D mod k`.  Equivalently, the dual diagonal `mu_k` acts on the entire degree-
`D` piece by one scalar.

For the four-vertex path, THM-2770's graceful obstruction has degree `12`.
The ternary path quotient therefore sees

```text
12=0 mod 3.                                              (17)
```

It acts trivially on every monomial and cannot split the zero balanced
coefficient into noncancelling character sectors.  Any signed coefficient
repair needs a non-diagonal reference phase, factor-selection history, or
another sidecar; the Smith character alone cannot provide it.

## 6. Exact verification and scope

Run

```bash
python 04-computation/one_long_wall_signed_imbalance_thm2776.py
python -O 04-computation/one_long_wall_signed_imbalance_thm2776.py
```

The exact companion uses explicit exceptions, integer arithmetic, and no
truth-bearing Python assertions.  It exhausts all `390,242` selections of
`k-1` positive `B_k` roots for `2<=k<=6`.  It independently classifies the
signed components, checks `(10)` on all `277,730` rank-`k-1` selections,
checks `35,358` balanced full-rank zero frames, verifies the exact
determinant histograms and attainable sets `(15)`, and checks the difference-
tree extremizers.  Normal and optimized runs byte-match the stored
transcript.

```text
PROVED HERE (candidate):  unique deficient signed-tree component;
                          anchored-tree/unbalanced-cycle square blocks;
                          exact determinant factorization (10);
                          full-rank zero boundary (11);
                          complete attainable-magnitude formula (12);
                          odd-prime bound and spanning-tree extremizer;
                          homogeneous diagonal-character no-go.

NOT PROVED:               a global arrangement cokernel;
                          a signed graceful coefficient selector;
                          graceful-label existence;
                          a modular-group or quartic-monodromy action;
                          Graceful Tree, JC(2), DC(2), or LRC(14).       (18)
```

QED (candidate).
