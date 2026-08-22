# Exact-six divisor lift and effective-kernel sidecars

## Status

**PROVED finite graph reduction + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This note proves an elementary quotient-lifting lemma and freezes its complete
application to the current refined exact-six ledgers.  It is not a theorem-ID
dependency, physical mutation, row exclusion, or LRC(14) terminal.

## Inheritance and corrected question

The raw exact-six relation has vertices `F in C({1,...,14},6)` and an edge
`F -> C` whenever `C` covers the strict unsupported target of at least one
retained divisor row `(F,D)`.  The refined SCC sidecar showed one body-level
three-cycle at `k=2` and no body-level cycle at `k=3`, but it still projected
away `D`.  The corrected question is therefore not whether a body arrow
exists.  It is whether consecutive arrows have compatible row and occurrence
sidecars.

The canonical hostile is the `k=2` cycle

```text
A=(1,3,7,8,9,10) -> C=(4,5,7,8,9,10)
 -> B=(2,4,5,7,9,12) -> A.
```

Its three unique source divisors are `17640,17640,4410`.

## 1. Exact current universes

The companion reconstructs the frozen support rows, current generating-
function occurrence weights, all strict endpoints and open atoms, and every
pool-14 six-subset.  The `k=3` post-screen universe has exactly `1897`
body/divisor keys and `2,548,893,834` occurrences:

| exact completion type | keys | occurrences |
|---|---:|---:|
| unique self completion | 1823 | 2,547,058,578 |
| unique nonself completion | 20 | 313,120 |
| no pool-14 six-completion | 54 | 1,522,136 |

Every self key is full-period, `D=L(F)`.  Every nonself key is half-period,
`D=L(F)/2`.  The twenty nonself arrows involve `39` bodies.

The exact Boolean potential

```text
Phi(F)=1_[1 in F]+2*1_[5 in F]-2*1_[10 in F]             (1)
```

strictly increases on all twenty arrows, with increment histogram

```text
1:3, 2:8, 3:4, 4:1, 5:4.                                (2)
```

Hence the body graph is a DAG.  Its unique longest path has two arrows:

```text
(1,2,3,7,10,12)
 -> (1,5,6,7,8,10)
 -> (3,4,5,7,8,12),                                     (3)
```

with potentials `-1,1,2`.  Even `(3)` is divisor-spliced: its two source
divisors are `2940` and `5880`.

## 2. Divisor and feature lifts

Retain `D` and require it to remain fixed across a transition.  For the twenty
`k=3` arrows, exactly eight target keys are absent from the current ledger and
twelve are full-period self keys.  Thus every live same-divisor nonself move
has height one: it lands at a self state.

The next refinement retains the exact generating-function feature class:

```text
(k,D,pattern,c,capacity)
```

and, on the one-spike slice, also its spike divisor and allowance.  Twelve
arrows transport `38` common feature classes representing exactly `298227`
occurrences.  All twelve again land at self states.  The other eight arrows
transport no feature because the target body/divisor key is absent.

This is a feature-class join, not an occurrence-identity theorem.  A feature
class still aggregates literal denominator multisets and discards phase.

For `k=2`, the same-divisor lift of the body cycle is

```text
(A,17640) -> (C,17640) -> (B,17640) -> self,
(B,4410)  -> (A,4410),
```

where `(A,4410)` has zero current occurrence weight.  The body `C3` therefore
arises exactly by splicing the second chain back to the first after forgetting
`D`.  One bit distinguishing `4410` from `17640` is necessary and sufficient
for this named cycle.  The literal divisor is sufficient but is not proved
globally minimal; even the coindices `L/D=(2,2,4)` detect this splice.

## 3. Kernel-pair lifting lemma

Let `E => X` be a directed multigraph and `q:X -> Y` a quotient.  A sequence
of edges can be composed after quotienting when

```text
q(t(e_i))=q(s(e_(i+1))).                                  (4)
```

Its seam defects lie in the kernel pair

```text
K_q=X times_Y X,
d_i=(t(e_i),s(e_(i+1))).                                  (5)
```

The quotient word is an actual lifted walk if and only if every `d_i` is on
the diagonal of `X`.  This is immediate because composition upstairs means
literal endpoint equality.  Thus every quotient-created cycle is assembled
from at least one off-diagonal seam in `K_q`.

A sidecar `c:X -> S` reflects all relevant length-two compositions exactly
when

```text
q(x)=q(x'), c(x)=c(x'), and (x,x') is a possible seam
 imply x=x'.                                               (6)
```

For finite `X`, make the conflict graph whose edges are the off-diagonal
possible seams in one `q`-fibre.  Condition `(6)` says precisely that `c` is a
proper coloring of this conflict graph.  Consequently its chromatic number is
the minimum number of labels needed to reflect the specified path family.
This proves why “minimal sidecar” is meaningless until the observable and
allowed seams are fixed.

## 4. Effective forgotten kernels across three frontiers

The same proof strategy, not the same algebraic object, appears elsewhere.
For a quotient `q:G -> A` and retained action `rho:G -> H`, the coarsest group
quotient carrying both is

```text
E_min=G/(ker q intersect ker rho),
K_eff=ker q/(ker q intersect ker rho) ~= rho(ker q).       (7)
```

The action descends to `A` exactly when `K_eff` is trivial.

- Here, the defect is an off-diagonal divisor seam in a finite transition
  kernel pair.
- In THM-3450, rational double centring keeps `Q(zeta_91)` but discards the
  nonzero `Q(zeta_13)` and `Q(zeta_7)` margin modules seen at orders eight and
  fourteen.  Through order fourteen those two irreducible modules are the
  minimal rational representation sidecar for the displayed obstruction.
- In THM-3449, abelian exponent coordinates forget the first visible
  commutator.  Through `3c`, the effective kernel is the central Heisenberg
  digit; at `3c+1` a higher bracket proves that repair no longer sufficient.

The divisor, margin modules, and commutator digit are not one object.  The
proved commonality is only failure of descent through a quotient and recovery
of the effective forgotten kernel.  The exact-six seam is not a cocycle
without additional torsor structure, and the D5 margin obstruction is not
graph holonomy.

## 5. Loss boundary and reproduction

The feature lift still destroys literal denominator tuples, common phase,
next-sector alignment, quotient-tail labels, inherited-clock distinctness,
ancestry, endpoint owners, and physical time.  Seven inherited clocks plus
six complements still make thirteen clocks, the open LRC(14) boundary.

Run

```bash
python -B 04-computation/lrc14_exact_six_divisor_feature_lift_20260815.py
python -B -O 04-computation/lrc14_exact_six_divisor_feature_lift_20260815.py
```

The companion is exact `Fraction`/integer/bitset/GF arithmetic and compiles
the targeted `k=2` engine independently under `-O2` and `-O3`.  Ordinary and
optimized Python runs reproduce the stored transcript.  LF-normalized hashes:

```text
script    e20dd7dda94afef1d8a044fd130832018f8b5c8b6859bb8a078b6f6dbaa16943
output    eddf99812bceda20fddef6d5cad649232b1ab4d3796ef2a96d1e1329ae17c0a1
semantic  c423bca65a43fdbfddc25b972a373002890f2c0a6f42c6aea96a8c2e2dacb8e5
```

An independent audit reconstructed the universes, all completions, divisor
rows, occurrence weights, potential, longest path, and stopping boundary.
