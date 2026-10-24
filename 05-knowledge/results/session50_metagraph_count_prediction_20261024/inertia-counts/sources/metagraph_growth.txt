# Additive growth, canonical inflation, and the missing complement mask

**Status:** PROVED here for the canonical inflation and complement-mask
identities; inherited PROVED results are named below. FINITE-EXACT companion
checks cover the explicitly stated finite universes. No general recurrence
for A000568, no classification of all metagraph components, and no asymptotic
complexity improvement is claimed.

## Inheritance and the exact object

- **Closest proved mechanism:** [THM-796, three-sorted recursive incidence](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md),
  especially sections 2, 3, 9, 10: additive recursion already exists on
  tilings and complement lines, with boundary phase retained.
- **Canonical hostile example:** endpoint deletion does not descend to a
  function on unmarked tournament classes. An explicit order-four witness
  appears below; THM-796 has the stronger node-kernel non-lumpability audit.
- **Corrected near miss:** uniform transitive inflation embeds isomorphism
  classes and respects converse, but ordinary all-free-bit complement leaves
  its image already at order four.
- **Relevant sidecar recovered:** [THM-4463, homogeneous tournament contacts](../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md)
  section 2 already proves that homogeneous-pair components are the unique
  maximal transitive modules. We use that canonical decomposition; we do not
  claim to have discovered it here.
- **Existing multiplicative response theory:** [THM-3235, transitive tournament
  blowup saddle and observer twins](../../../01-canon/theorems/THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation.md)
  already studies this inflation's determinant and Hamiltonian response.
  Its observer twins warn that a few scalar responses do not determine the
  tournament as a substitution operator.

Use path order `1 -> 2 -> ... -> n`, the relabelled reverse of THM-796's
convention. Let `X_n` be its fixed-path tournament cube; its free coordinates
are pairs `(i,j)` with `j>=i+2`. Write

```text
m_n=binom(n-1,2),           |X_n|=2^m_n,
tau_n = complement all free coordinates,
sigma_n(i,j)=(n+1-j,n+1-i).
```

The converse-merged node set is `I_n/(T~T^op)`, where `I_n` consists of
ordinary isomorphism classes. A line is a `tau_n` orbit; it is blue exactly
when its tilings are fixed by `sigma_n`, and black otherwise. These are the
literal definitions of [THM-643, blue/black line parity](../../../01-canon/theorems/THM-643-blue-black-line-parity-structure.md)
and [THM-1430, tiling/class/metagraph dictionary](../../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md).

The board is: **boundary phase; canonical transitive modules; primitive
scale; complement defect; full Hamiltonian-path fibre**. Each map below
specifies which of these survives.

## 1. The additive laws already proved, and their boundary

Right-end growth from `X_n` to `X_(n+1)` freely adds `n-1` coordinates.
Consequently each parent has `2^(n-1)` right-end lifts, and two successive
right-end additions have `2^(2n-1)` lifts. Right deletion commutes with
`tau`; `sigma` exchanges right and left deletion. This is part of the
stronger compatible-two-face pullback in THM-796 section 2.

Adding one vertex at each end gives the middle restriction
`X_(n+2) -> X_n`, which commutes with **both** involutions. Its `2n-1`
new free coordinates consist under `sigma` of `n-1` exchanged pairs and
one fixed outer chord. Thus a blue parent tiling has `2^n` blue lifts;
a black parent has none.

For `n>=3`, the same statement on complement-line instances is

```text
                     upper blue        upper black
lower blue             2^n            2^(2n-1)-2^n
lower black             0                 2^(2n-1).
```

Entries count lifts **per individual lower line**, with line instances
retained. This is precisely THM-796 section 10 in parent-order notation.
It yields the total blue recurrence `B_(n+2)=2^n B_n`, but does not tell
which isomorphism classes contain those upper lines.

For a small witness to that missing class information, take the tournament
with arcs `12,23,34,13,24,41`. It has Hamiltonian paths `(1,2,3,4)` and
`(3,4,1,2)`. Deleting the last vertex of the first path leaves a transitive
triple. Deleting the last vertex of the second leaves the cycle `3->4->1->3`.
Thus even one fixed order-four class has nonisomorphic endpoint children.
Retain a marked Hamiltonian path or the full deletion correspondence with
multiplicities; an unmarked node alone is insufficient.

## 2. Unique uniform transitive roots

For an integer `k>=1`, define

```text
D_k(T)=T[K_k],
```

where `K_k` here means the **transitive tournament** of order `k`:
replace every vertex by an ordered chain of length `k`, and every base arc
by all `k^2` arcs between its two chains. This is not the unrelated skew
Hadamard doubling used elsewhere in the repository.

Let `C_1,...,C_b` be the canonical maximal transitive modules from THM-4463,
with lengths `l_1,...,l_b`. Their quotient `Q` has no homogeneous pair.
Indeed, two homogeneous quotient vertices would correspond to two modules
whose ordered union is a larger transitive module.

**Proposition 1 (unique root criterion).** The maximal transitive module
lengths of `D_k(T)` are exactly `k*l_1,...,k*l_b`, with the same quotient
`Q`. A tournament `S` is isomorphic to some `D_k(T)` if and only if every
one of its canonical module lengths is divisible by `k`. Its root is then
unique up to isomorphism.

**Proof.** Inflating a canonical chain produces a homogeneous transitive
chain of multiplied length. A homogeneous pair joining two distinct such
chains would force their quotient vertices to have identical incidence to
every other quotient vertex, contradicting the absence of quotient
homogeneous pairs. Since canonical chains are the components of the
homogeneous-pair graph, no merger is possible. Conversely,
given divisible canonical lengths, divide them all by `k` while retaining
`Q`; reinflation reconstructs `S`. The quotient and ordered module lengths
are canonical under isomorphism, so this operation is unique. This also
shows that every automorphism acts on the weighted quotient and is unique
on each transitive chain. Multiplying all weights by `k` changes neither
the allowed quotient automorphisms nor their action, so
`Aut(D_k(T))` is naturally isomorphic to `Aut(T)`. ∎

Converse reverses the quotient and the order inside each chain without
changing chain lengths. It commutes with inflation up to the evident
isomorphism. Therefore `D_k` is injective on ordinary isomorphism classes,
on self-converse classes, and on converse-merged classes; self-converseness
is preserved **in both directions**.

Define the intrinsic scale

```text
g(T)=gcd(l_1,...,l_b).
```

Every nonempty tournament has a unique expression `D_g(P)` with `g(P)=1`.
If `A(n)` counts ordinary tournament classes of positive order `n`, and
`P(n)` counts those with scale one, then

```text
A(n)=sum_(d|n) P(d),
P(n)=sum_(d|n) mu(d) A(n/d).                              (1)
```

This is a canonical divisor stratification, not a closed recurrence for
the unknown primitive counts. The identical divisor law holds separately
for self-converse and converse-merged counts, using the corresponding
primitive strata. In particular, the number of order-`kn` classes admitting
a uniform `k`th root is **exactly `A(n)`**.

Small exact controls are

| n | A(n) | P(n) | Nonprimitive scale counts |
|---:|---:|---:|---|
| 1 | 1 | 1 | none |
| 2 | 1 | 0 | scale 2: 1 |
| 3 | 2 | 1 | scale 3: 1 |
| 4 | 4 | 3 | scale 4: 1 |
| 5 | 12 | 11 | scale 5: 1 |
| 6 | 56 | 54 | scale 2: 1; scale 6: 1 |

Scale one does not mean that the tournament has no nontrivial modules.
It only means their **length gcd** is one. It is weaker than being
indecomposable, and must not be substituted for that notion.

## 3. The exact complement defect of multiplication

In the fixed-path presentation, place the `k` clones of vertex 1 first,
then the `k` clones of vertex 2, and so forth. This yields an injection
`D_k:X_n -> X_(kn)`. It repeats each old free coordinate on a `k` by `k`
rectangle. All other target coordinates are fixed to the forward direction.
It satisfies

```text
sigma_(kn) D_k = D_k sigma_n.                            (2)
```

Let `E_(n,k)` be the following fixed set of target **free** coordinates:

1. nonadjacent edges internal to each clone block;
2. edges between two adjacent base blocks, except their single connecting
   Hamiltonian-path edge.

Then

```text
tau_(kn) D_k(t) = D_k(tau_n(t)) XOR E_(n,k),              (3)

|E_(n,k)| = n*binom(k-1,2)+(n-1)*(k^2-1)
           = (k-1)*(3*n*k-2*k-2)/2.                     (4)
```

**Proof.** Edges between nonadjacent base blocks are exactly the repeated
old free coordinates; both sides of (3) complement them. Every other
target free coordinate belongs to the displayed set `E`: the left side
complements it while the inflation on the right fixes it forward.
There are `binom(k-1,2)` nonpath edges in each length-`k` block and
`k^2-1` nonpath edges between each adjacent pair of blocks. This proves
both identities. Reflection preserves this partition, so `sigma(E)=E`. ∎

The corrected involution

```text
hat_tau_(n,k)(s)=tau_(kn)(s) XOR E_(n,k)
```

preserves the inflated subcube and satisfies
`hat_tau_(n,k) D_k = D_k tau_n`. For `n>=3`, it is fixed-point-free there,
and (2) preserves blue versus black tilings in both directions. Therefore
inflation is an exact embedding of the tiling/line incidence structure
**when the corrected complement is retained**.

For doubling, the defect is particularly economical:
`|E_(n,2)|=3(n-1)`. There are three defective free edges at every boundary
between consecutive clone pairs. This gives a precise multiplicative
boundary sidecar instead of an unspecified failure of self-similarity.

The correction is necessary. `D_2(K_2)` is the transitive order-four
tournament, with one canonical chain of length four. Its ordinary `tau`
mate has arcs `12,23,34,31,41,42` and chain lengths `(2,1,1)`. It has no
uniform square root. Thus ordinary black/blue complement lines do not
stay in the doubling image. The corrected edges are generalized flips,
not ordinary edges of the original metagraph.

## 4. A second missing coordinate: non-block Hamiltonian paths

Even the corrected operation does not identify the **entire** metagraph
induced on inflated isomorphism nodes. The injection above retains only
Hamiltonian paths whose clone blocks remain consecutive.

For the directed triangle `C3`,

```text
                    H(T)     |Aut(T)|       tiling fibre
C3                    3          3               1
D_2(C3)              45          3              15.
```

The fifteen target presentations include fourteen outside the fixed
block-respecting image of the one source presentation. They arise from
Hamiltonian paths that interleave clone blocks. These finite counts are
checked directly over all `6!` vertex orders. The automorphism equality
also follows generally from Proposition 1.

Consequently the faithful object for multiplicative reconstruction must
retain **both** the complement defect and the full Hamiltonian-path fibre.
The canonical root alone gives exact node strata; it does not determine
all line incidence at those nodes.

## 5. The precise remaining research target

The useful new map is

```text
source: (T, Hamiltonian path, canonical chain lengths)
target: inflated tournament and its block path
preserved: isomorphism/converse, canonical quotient, scaled lengths,
           grid symmetry, corrected-complement incidence
lost by bare-node projection: interleaving HPs and endpoint phase
required sidecars: E_(n,k), full HP-orbit fibre
cheapest next test: classify the additional HP interleavings by
                   the sequence of visited clone blocks.
```

Such a sequence is a walk in the base tournament **with a loop added at
each vertex**, with each vertex appearing exactly `k` times. Each maximal
consecutive run in one block must use its clone labels in increasing
order. Separated visits to that block need not occur in label order.
If `W_k(T)` is the set of these block words, then the elementary count is

```text
H(D_k(T)) = sum_(w in W_k(T)) (k!)^n / product_(runs R of w) |R|! . (5)
```

Indeed, split the `k` clone labels at each vertex among its ordered runs;
the multinomial coefficient counts the choices, and each run has one
allowed ascending order. Multiplying these coefficients gives (5).
The formula counts Hamiltonian paths, not their automorphism orbits;
division by the unchanged `|Aut(T)|` gives the target tiling fibre.

This is a direct bridge to the existing walk-content substitution kernel
in THM-3235 and its dependencies. However, the cheap hostile test already
shows that the block word alone, even with all run lengths, is insufficient
for ordinary complement incidence.

In `D_2(C3)`, call the clone chains `(A1,A2)`, `(B1,B2)`, `(C1,C2)`, with
base arcs `A->B->C->A`. Both

```text
P = (A1,B1,B2,C1,C2,A2),
Q = (A2,B1,B2,C1,C2,A1)
```

are Hamiltonian paths with the same block word `ABBCCA`. Complementing
all nonpath arcs in the respective path coordinates gives degree multisets

```text
tau_P: (2,2,2,3,3,3),
tau_Q: (1,2,2,3,3,4).
```

They are nonisomorphic even after allowing converse, since both displayed
multisets are invariant under `d -> 5-d` and are distinct. The first failed
implication is **same block word implies same complement target class**.
The strongest survivor is the exact weighted path count (5). The missing
coordinate is which ordered clone is assigned to each separated visit.

The next question is whether clone allocations admit an operation-stable
quotient smaller than the complete path, while still determining ordinary
complement's target class. This is a concrete finite-state minimization
target; the assertion of a compact sufficient state remains OPEN.

Uniform inflation gives a genuine arithmetic scale decomposition of the
nodes. The difficult primitive classes remain; equation (1) isolates them
instead of counting them for free.

## Reproduction and finite scope

```text
python 05-knowledge/results/metagraph_growth_20261020/growth_check.py
python -O 05-knowledge/results/metagraph_growth_20261020/growth_check.py
```

The companion has no `assert` statements. Its controls include all 1,099
labelled tournaments of orders 1 through 5, comparing homogeneous-pair
components against an independent all-subset maximal-transitive-module
test; complete fixed-path cubes through order 6 for the canonical-root
census; 152 inflation cases at multipliers 2 and 3; an independent full
permutation check of 76 class representatives; three hostile witnesses;
16 block-run multinomial counts checked against independent Hamiltonian
subset DP; and all 33,856 target tilings of orders 5 through 7 for the
inherited two-ended colour channel. Normal and optimized Python passed
with identical displayed output. The recorded `.out` is literal program
output.
