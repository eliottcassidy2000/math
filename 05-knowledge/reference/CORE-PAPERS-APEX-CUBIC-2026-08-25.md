# Apex cubic three-edge-coloring: source audit and transfer boundary

> **Audit date: 2026-08-25. Status: CITED / EXTERNAL PREPRINT UNDER
> AUDIT.** The source is a one-day-old, computer-assisted arXiv v1.  This
> repository has audited the statement and proof interfaces but has not
> replayed its three large computational obligations.  Do not promote the
> external theorem to internally `PROVED` on this intake alone.

## Source identity

- **Primary:** Yuta Inoue, Ken-ichi Kawarabayashi, Ritarou Matsuo, Atsuyuki
  Miyashita, Bojan Mohar, and Tomohiro Sonobe,
  [*Three-edge-coloring apex cubic graphs*](https://arxiv.org/abs/2608.22870),
  arXiv:2608.22870v1 [math.CO], submitted 2026-08-24.
- **Full text:** [official arXiv HTML](https://arxiv.org/html/2608.22870v1).
- **Code/data:** the authors' public
  [GitHub organization](https://github.com/three-edge-coloring-apex-cubic-graphs)
  separates configurations, discharging rules, computer checks, and the
  semi-reducibility checker.
- **Metadata discrepancy:** the arXiv abstract page currently renders the
  third author as “Rintaro Matsuo,” while the paper body renders “Ritarou
  Matsuo.”  This intake follows the PDF/body spelling and records the mismatch
  rather than silently choosing an identity claim.

## Exact imported claim

Theorem 1.4 states:

```text
Every 2-connected apex cubic graph is three-edge-colorable.              (1)
```

Here `apex` means that deleting one vertex makes the graph planar.  Section 8
claims a constructive `O(n^2)` algorithm: find an apex by planar deletion
tests, locate one member of the fixed unavoidable configuration family,
reduce recursively, and extend the coloring using a stored reducibility
certificate.

The paper combines `(1)` with the Robertson--Seymour--Thomas reduction and
the published doublecross case to claim the final piece of Tutte's **cubic
Petersen-minor-free three-edge-coloring conjecture**.  It does not settle
Tutte's general 4-flow conjecture for arbitrary noncubic graphs.

Under the paper's definition of a snark—cyclically 4-edge-connected, cubic,
girth at least five, and not three-edge-colorable—`(1)` implies:

```text
no snark is apex; equivalently, a snark needs at least two vertex deletions
to become planar.                                                         (2)
```

This does not say that every apex cubic graph is Petersen-minor-free, nor
does it classify snarks, resistance, oddness, or criticality.

## Proof architecture

### 1. Bounded-defect planarization

From a smallest noncolorable apex cubic graph `G_mc`, Section 2 deletes either
one apex-incident edge when that planarizes, or the apex itself.  The planar
subcubic core `G` has exactly two or three degree-two vertices; its dual
`G*` has exactly two or three digons.  Lemma 2.2 records minimum degree five,
at most one digon incident at a vertex, and explicit separation bounds for
cycles of length at most five.

The bounded defect ledger is load-bearing.  Deleting an arbitrary number of
exceptional vertices would not preserve the later constant-size charge
budget or state space.

### 2. Boundary-state greatest fixed point

Sections 3--4 use multi-boundary islands, with at most three boundary faces.
An exterior two-color Kempe component acts on the ordered boundary colors.
Planarity makes its paired endpoints a noncrossing semi-matching; a degree-two
defect can terminate a component, creating a singleton half-chain.

A set of bad boundary colorings is `semi-consistent` when it is closed under
all subset switches for one admissible semi-matching.  Semi-D-reducibility
means that the greatest move-closed set of nonextendable words is empty.
Semi-C-reducibility uses a small deletable edge set and asks the bad fixed
point to avoid the corresponding restricted extension set.  Lemmas 3.7 and
3.10 exclude either kind of island from a minimal counterexample.

Lemma 3.15 proves a monotonicity property needed by the compiler: deleting
boundary leaves preserves both reducibility notions.  This is an actual
restriction/lift statement, not the heuristic that a smaller pattern should
be easier.

### 3. Defect-compensated discharging

The initial charge is `10(6-d(v))`.  Lemma 6.2 gives total charge `60` with
three digons and `80` with two.  The unavoidable family has `915` normal
configurations.  A local allowance `c(v)` compensates for nearby digons; the
five degree cases make one digon's total contribution at most

```text
18, 12, 16, 8, or 12,                                                   (3)
```

always strictly below `20`.  At most three digons therefore consume less
than the global surplus `60`, forcing a configuration homomorphism in
Theorem 6.7.

The portable principle is `global conserved surplus > sum local defect
allowances`.  Merely ignoring exceptional neighborhoods would not give the
contradiction.

### 4. Rotation-preserving quotient compiler

Weak connectivity and digons let a configuration embed non-inducedly or with
identified darts.  Section 7 enumerates free homomorphic images, prunes by
planarity and forbidden short separators, converts survivors to
multi-boundary islands, and checks their semi-reducibility.  Theorem 7.1
turns the unavoidable configuration from Theorem 6.7 into an excluded
reducible island.

These are homomorphisms of a dart/rotation representation: they preserve
heads, edge reversal, non-null successor/predecessor incidences, and specified
degrees.  An adjacency-only graph quotient loses cyclic order and boundary
nil pointers and is not a faithful replacement.

## Computer-assisted boundary

The mathematical proof isolates the finite obligations as Lemmas B.1--B.3.
The authors' [reproducibility
instructions](https://github.com/three-edge-coloring-apex-cubic-graphs/instructions-for-checking-reproducibility)
give these audit metrics:

```text
unblocked combined rules                         747
maximum combined-rule transfer                     5
cartwheels at center degrees 7,8,9,10,11
                                      4438,4939,2409,567,38
multi-boundary islands with 0,1,2,3 boundaries
                                      254,88393,20836,18.                 (4)
```

The paper reports independent implementations reconstructed from pseudocode
with AI assistance.  It explicitly treats those reconstructions as added
reproducibility evidence, not as part of the mathematical justification.
This repository has not rerun `(4)` and therefore records the whole external
claim as under audit.

## Repo consumer and exact sharpening

[THM-4116](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md)
imports no unverified coloring theorem.  It extracts the paper's native
interface operation and proves it independently:

```text
# three-edge-colorings of a glued graph
  = dot product of its ordered boundary extension vectors.               (5)
```

Its exact Petersen four-edge boundary has disjoint supports and its K4
control has dot product six.  Cubic parity plus a Kempe involution compresses
every adjacent-edge four-pole of an uncolorable finite simple cubic graph to
one even multiplicity.  Exact atlases distinguish the Petersen, `J_5`, and
two Blanusa graph families, but that scalar does not in general classify
their edge orbits.  Transporting
`(5)` as an interface-state operation to the
already-proved THM-4110 AP13 torsion quotient yields the new odd-shell
component law `2^(c(F)-1)` and classifies all minimal phase synchronizers as
spanning trees.  That sharpening is internally `PROVED`; it does not validate
the external computer-assisted theorem `(1)`.

## Transfer rules and stopping conditions

The mechanisms worth reusing are:

1. **bounded-defect reduction:** tractable core plus an explicit finite defect
   sidecar and a proved lift back;
2. **boundary action carrier:** compute the greatest native-move-closed bad
   set, including half-moves when a conservation law fails;
3. **multi-boundary quotient closure:** enumerate folds while retaining
   rotation and attachment data;
4. **defect subsidy:** dominate every exceptional neighborhood by a local
   allowance whose total is below a conserved surplus; and
5. **unavoidable plus reducible compiler:** store a decreasing rewrite and
   its inverse extension certificate, yielding construction as well as
   existence.

The state dynamics is a reversible action hypergraph/groupoid, not a
tournament.  A Kempe move flips any subset of components from one
noncrossing semi-matching; orienting pairwise comparisons loses the
simultaneous-subset operation and is not sound.

## Does not prove

This intake and THM-4113 prove none of:

- the external computer checks or the full theorem `(1)` internally;
- Tutte's general 4-flow conjecture;
- a classification or enumeration of snarks;
- a `k`-apex extension of the paper;
- a new Euclidean Kakeya estimate;
- a planar Jacobian reduction;
- physical LRC safety from phase synchronization; or
- LRC(14).
