---
id: THM-3269
title: "Scale-invariant clutch strength and canonical weighted bispanning polarization"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  The exact reset-link trap intervals canonically weight every edge of
  THM-3260's bispanning core by a positive row-scale- and orientation-invariant
  overlap ratio.  Among all 9,920 ordered complementary-tree charts, the
  product weight has one exact minimizer.  Hence the analytic clutch data
  selects an ordered integral polarization and removes THM-3260's external
  tree-pair sidecar.  It does not supply the remaining C12 vertex label.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-03
audit: >
  The assertion-independent exact companion pins promoted THM-3254 and
  THM-3260.  It recovers THM-3254's directly reconstructed response bank but
  trusts none of its chosen covers; rebuilds every endpoint trap and overlap
  strength; verifies orientation invariance; independently enumerates all
  complementary spanning trees and both orderings; proves the minimum unique
  by exact rational comparison; compares the intrinsic C2 image; and rebuilds
  the selected chart's two unimodular incidence matrices and GL11 transition.
  It also exhausts the selected tree's degree-compatible automorphisms and
  checks its unique center, radius, diameter and leaf set.
  Normal, optimized and stored replay plus LF hashes are required.
  Independent hostile audit is pending.
depends_on:
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
  - THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary
related:
  - THM-3255-twelve-balance-multiplicative-singer-rank-defect-and-phase-marker-boundary
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
script: 04-computation/gmc_scale_invariant_weighted_bispanning_thm3269.py
output: 05-knowledge/results/gmc_scale_invariant_weighted_bispanning_thm3269.out
script_sha256: ff4583490afbeefb6a507215604b697d412c1ebd7c547b32d05c9dac19f0c474
output_sha256: 721ba5ec05adbeeba6bf9cf6dceb1fe5e4c1585ad999cf2c240c04284263f645
hash_basis: LF-normalized bytes
---

# THM-3269 -- scale-invariant clutch strength and canonical weighted bispanning polarization

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-3260 proves that the 22-edge reset-link core has 4,960 complementary
spanning-tree decompositions, connected by symmetric exchanges.  Any ordered
decomposition gives the integral polarization needed by its conditional
bridge to the eleven missing `C_12` modes, but the bare graph distinguishes no
chart.

The analytic clutch intervals do distinguish one.  The selector is intrinsic
under the harmless choices that a response atlas actually permits: positive
rescaling of either row and reversal of an undirected edge.

## 1. A scale-free weight on a response edge

For a first-link-blocked pair `e={i,j}`, use the orientation `(i,j)` and let

```text
I_sigma(i,j) subset [0,infinity]                       (1)
```

be THM-3254's closed interval of positive ratios `lambda` for which the reset
link state `sigma` traps `lambda f_i+f_j`.  Let

```text
A_e={U : I_sigma(i,j)=[0,U] for some sigma},
B_e={L : I_tau(i,j)=[L,infinity] for some tau}.        (2)
```

THM-3254 proves that every one of the 23 reset-link edges has a sharp
two-interval cover of the positive line and no one-state cover.  Thus the two
finite sets in `(2)` are nonempty and some `U>=L>0`.  Define the **clutch
strength**

```text
kappa_e=max {U/L : U in A_e, L in B_e}.                (3)
```

This is a positive rational number at least one.

The definition is independent of row normalization.  Replacing
`f_i,f_j` by `alpha f_i,beta f_j`, with `alpha,beta>0`, multiplies every ratio
endpoint by the common factor `beta/alpha`; every quotient `U/L` is unchanged.

It is also independent of edge orientation.  For the reversed pair,

```text
mu f_j+f_i = mu F_(1/mu).
```

Hence `[0,U]` becomes `[1/U,infinity]` and `[L,infinity]` becomes
`[0,1/L]`.  The corresponding overlap is

```text
(1/L)/(1/U)=U/L,                                      (4)
```

and inversion bijects all endpoint pairs.  Therefore `kappa_e` is an
invariant of the undirected, positively projectivized response edge.

On THM-3260's 22-edge leafless core, exact reconstruction gives

```text
12654135326210/11849988963879
  <= kappa_e <=
66479967673052544/306432210196733.                    (5)
```

The complete strength bank, including its maximizing state pairs, has digest

```text
c7b8441f54bd893f33f1767db913da0d40f392f20d74f6d177782a0490dc911a. (6)
```

## 2. The unique weighted chart

For an ordered complementary-tree chart `(T,T')`, put

```text
W(T)=product_(e in T) kappa_e.                         (7)
```

This avoids logarithms: comparing `W` is an exact rational version of
minimizing the sum of the logarithmic clutch strengths.  Exhaustion of all
4,960 unordered decompositions and both orderings gives a unique minimizer,

```text
T_*={(2,13),(2,17),(2,19),(2,21),(3,11),(3,16),
     (7,22),(10,11),(16,17),(18,19),(21,22)}.          (8)
```

Its canonical complement is

```text
T_*'={(2,10),(3,22),(7,18),(10,16),(10,22),(11,13),
      (13,18),(13,22),(16,21),(17,22),(19,22)}.        (9)
```

Both are spanning trees.  The ratio between the second-smallest ordered
weight and the minimum is exactly

```text
111771229679330735594490655
-------------------------------------------  > 1.     (10)
106891942612503954628585203
```

Thus uniqueness is an exact separation, not a floating-point tie break.
Since the product over all 22 edges is fixed, `(8)` simultaneously orders the
decomposition: `T_*` is the weak-clutch tree and `T_*'` the strong-clutch
cotree.

The nontrivial automorphism of the unweighted core swaps rows 17 and 21.  Its
image of `T_*` is not the minimizer; its exact weight ratio is

```text
W(C_2 T_*)/W(T_*)
 =719636186849558063/684759823142197440 > 1,           (11)
```

and it is the fourth ordered chart in strict weight order.  Hence the analytic
weights break the whole intrinsic `C_2`; the weighted response graph has no
residual chart ambiguity.

The selected weak tree is itself rigid:

```text
Aut(T_*)=1,     center(T_*)={17},
radius(T_*)=4,  diameter(T_*)=8,
leaves(T_*)={7,10,13,18}.                              (11a)
```

Thus the analytic selector canonically supplies row 17 as a root as well as
the ordered polarization.  Rigidity does **not** by itself provide a cyclic
ordering: imposing a generic graph-canonicalization algorithm and numbering
its output would be a convention, not an intertwiner with the physical phase
module.

## 3. Canonical integral polarization

Root at the canonical center row 17 and orient every edge by increasing
label.  The two reduced
incidence matrices for `(8)--(9)` satisfy

```text
det B_(T_*)=1,       det B_(T_*')=-1.                 (12)
```

Consequently

```text
U_*=B_(T_*)^(-1) B_(T_*') in GL_11(Z),
det U_*=-1.                                            (13)
```

The exact matrix digest is

```text
93f9ae1643ad9606f15826d8f1d3809fd0b04ce1d74c11d92d354edbdd737d4f. (14)
```

As in THM-3260, the columns `(-U_*z,z)` give an integral cycle frame and the
two incidence matrices give dual cut coordinates.  Equations `(3)` and `(8)`
make this polarization canonical from the normalized response problem itself;
no external tree-pair choice remains.

## 4. Consequence for the `C_12` bridge

THM-3260's integral bridge

```text
Aug(Z[C_12]) --> H_1(G_0;Z)                           (15)
```

was conditional on two sidecars: a cyclic vertex label with a selected root,
and an ordered tree-pair chart.  The chart `(T_*,T_*')` and its center 17
supply the second sidecar and the root canonically.  Therefore `(15)` now
needs only an external bijection

```text
V(G_0) <--> C_12,       17 |--> 0.                    (16)
```

The load-bearing cyclic-difference map `Z[C_12]/Z N -> Aug(Z[C_12])` from
THM-3260 remains unchanged.

This is a genuine positive holotopy move: the unweighted exchange atlas has
no preferred origin, while the exact analytic clutch supplies a discrete
potential with one global minimum.  It selects a chart; it does not create
phase labels.

## 5. Failure boundary and scope

The selector uses only THM-3254's eleven reset-link states and the 22
first-link-blocked core edges.  It is special to the complete support-(1,3),
bank-I2 response normalization up to positive projective row scaling.

The argument would fail if the minimum in `(7)` were tied.  Merely knowing
the graph, its Betti number, its critical group, or the unordered interval
cover does not imply uniqueness; the exact rational endpoint values are
load-bearing.

The theorem does not provide the cyclic vertex label `(16)`, an owner/phase
map, a `C_7` carrier on the delayed relative layer, a positive physical
current, or a scalar Gaussian-moment response.  It proves no row exclusion,
no arbitrary-radial NC2 theorem, no Gaussian Moment Conjecture consequence,
and no `LRC(14)` decrement.

## 6. Exact verification

Run

```text
python 04-computation/gmc_scale_invariant_weighted_bispanning_thm3269.py
python -O 04-computation/gmc_scale_invariant_weighted_bispanning_thm3269.py
```

and compare LF-normalized bytes with

```text
05-knowledge/results/gmc_scale_invariant_weighted_bispanning_thm3269.out.
```

The companion uses exact integer and rational arithmetic only.  It has no
floating point, randomness, discovery cache, graph-library dependency, or
optimization-sensitive assertion.  It independently reconstructs all
complementary tree pairs and the selected integral transition.

QED, conditional only on the pending independent audit required for status
promotion.
