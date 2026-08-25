> **RESEARCH SYNTHESIS (2026-08-25).**  Internal claims route to proved
> THM-261, THM-4113, and THM-4116.  The apex-cubic theorem of
> arXiv:2608.22870v1 is `CITED / EXTERNAL PREPRINT UNDER AUDIT`; its large
> computer obligations were not replayed here.  LRC(14), JC(2), and the
> Euclidean four-dimensional Kakeya conjecture remain OPEN.

# Apex Snarks and Boundary Compilers Across Frontiers

## Portfolio and inheritance pass

**Anchor -- LRC(14).**  THM-4110 had isolated 64 AP13 circle-phase sheets.
THM-4116 now identifies seven odd-shell vertices carrying six free bits
(speed one is anchored) and proves that a graph `F` of same-shell constraints
leaves exactly `2^(c(F)-1)` sheets.  This
closes phase synchronization, not physical safety.

**Niche -- snark interfaces.**  THM-261 gives the closest old exact carrier:
Petersen is the disjoint-support graph on positive `A_4` roots.  The new
interface carrier is stronger for coloring: an ordered boundary word is
retained with its extension multiplicity, and gluing is tensor contraction.

**Wildcard -- Jacobian and Kakeya boundary data.**  The response profiles of
THM-4103 and the tube/cell interfaces of four-dimensional Kakeya both invite
boundary-state descriptions, but neither currently has the lift lemma or
bounded defect budget needed to import the apex proof.

The canonical hostile is Petersen itself.  It has no three-edge-coloring,
but it is not cycle-poor:

```text
(c_5,c_6,c_8,c_9)=(12,10,15,20).
```

On the THM-4116 four-edge cut, both sides admit many local boundary states,
yet their supports are disjoint.  This distinguishes “no local solution”
from the genuine global failure “no compatible pair of local solutions.”

The corrected near misses are MISTAKE-069 (first-found versus smallest/full
cycle profile) and MISTAKE-507 (global Ω-girth, Lie-carrier, flow, and
snark/prime conflations).  The least-used sidecar is the full cycle-witness
incidence with the chosen cut; scalar cycle counts alone cannot locate the
coloring obstruction.

## Live concept board

| Concept | Representation and invariant | Native operation | Status / hostile |
|---|---|---|---|
| Boundary extension tensor | `f_X(sigma)` with ordered colors and multiplicity | tensor contraction along shared words | PROVED, THM-4116; Petersen dot zero, `K_4` dot six |
| Half-Kempe boundary atlas | maximal noncrossing partial matchings with singleton half-chains | simultaneous subset switches / greatest fixed point | PROVED combinatorics, THM-4113; coarse Petersen compatibility has 15 states versus 10 planar states |
| Bounded defect ledger | planar core, two or three digons, local allowances | discharging plus inverse reducibility rewrite | CITED external preprint; computation not replayed |
| AP odd-shell sheet graph | `F_2` bit on each odd speed, anchored at speed one | equality constraints and edge deletion | PROVED, THM-4116; a six-edge cyclic hostile leaves two sheets |
| Tournament deletion charge | `chi_x=2[(H-disc)(C)-(H-disc)(C-x)]` with rooted sign and insertion fibers | delete one vertex / sum local charges | identity PROVED, local unavoidability OPEN, HYP-9080; source-over-`C_3` pointwise hostile |
| Jacobian response packet | labelled infinity indices and Newton-face data modulo target shear | response-fiber intersection / chart gluing | OPEN; coarse degree packet need not realize |
| Kakeya tube-cell packet | scale, direction cap, tangency type, multiplicity | polynomial partition and cross-scale induction | OPEN; THM-4035 all-broad/all-narrow collision is hostile |

Every pull below is compared against these six concepts: does it produce a
finite interface, retain the target predicate, admit a native composition,
and come with a hostile object that detects quotient loss?

## What arXiv:2608.22870 contributes

The paper claims that every 2-connected apex cubic graph is three-edge-
colorable and gives a constructive quadratic-time algorithm.  The exact
proof architecture, rather than the headline alone, is the useful import:

1. delete the apex or one incident edge to obtain a planar subcubic core with
   exactly two or three degree-two defects;
2. encode exterior Kempe components as noncrossing boundary semi-matchings,
   allowing singleton half-chains at defects;
3. define reducibility as emptiness of the greatest move-closed set of bad
   boundary states;
4. use total charge `60` or `80` and a strict allowance below `20` per digon
   to force one of 915 configurations; and
5. enumerate rotation- and dart-preserving homomorphic images, reduce, and
   store the inverse coloring extension.

The boundedness is load-bearing.  “Delete some exceptional set” is not a
transfer unless the number of defects, the total subsidy, the boundary-state
size, and the lift back are all controlled.  Likewise, adjacency-only
quotients are too coarse: cyclic order, dart reversal, boundary face, and
nil attachments are part of the object.

The paper's move system is a reversible action hypergraph or groupoid.  A
single move can flip any subset of Kempe components belonging to one
noncrossing semi-matching.  There is no intrinsic complete antisymmetric
binary relation, so Tournament Analysis is rejected on this carrier.

## Exact sharpening obtained in this session

THM-4116 proves three layers independently of the external theorem.

### 1. Coloring is a boundary tensor contraction

For a vertex cut `X`,

```text
#Col_3(G)=sum_sigma f_X(sigma)f_(V-X)(sigma).
```

For any partition `V=V_1 disjoint_union ... disjoint_union V_k`, the same
restriction/union bijection gives

```text
#Col_3(G)=sum_tau product_i f_i(tau|delta(V_i)).
```

This is an associative local-to-global compiler.  It preserves extension
counts and exact compatibility, while discarding interior geometry.  Boundary
order, repeated interface-edge incidence, color gauge, and automorphism
action are required sidecars.

For the displayed Petersen four-edge boundary the left support has 12 multiplicity-one
words, the right support has nine multiplicity-two words, and the supports
are disjoint.  The analogous `K_4` cut has dot product six, one coloring
modulo the free global `S_3` action.  An independent referee additionally
checked the contraction identity on every simple graph through five vertices
and every vertex cut.

For an adjacent edge of any finite simple cubic graph, parity reduces the
complement four-pole tensor to four `S_3` word orbits `(S,P,X,Y)`, and
`#Col_3(G)=6(x+y)`.  If `G` is uncolorable, a Kempe-path involution further
forces the signature to `(m_e,m_e,0,0)` with even `m_e`.  The resulting exact
edge histograms are

```text
Petersen:   {2:15};
J_5:        {4:5, 6:10, 10:10, 12:5};
Blanusa I:  {4:23, 6:4};
Blanusa II: {4:25, 6:2}.
```

Thus the smallest interface already separates the three snark families by
extension multiplicity even though their allowable boundary support is the
same forced nine-word set.

### 2. AP odd-shell phase synchronization is connectivity

For `v=(1,...,2q-1)`, the complete unequal-`v_2` reciprocal graph leaves one
half-toggle bit on every odd speed other than the anchored speed one.  Adding
an odd-shell graph `F` imposes equality of bits along its edges, so

```text
K_(Gamma_F)(v)/P_v ~= F_2^(c(F)-1).
```

Therefore saturation is equivalent to connectivity, minimum synchronizers
are the `q^(q-2)` labelled trees, and AP13 has exactly `7^5=16,807` minimum
six-edge synchronizers.  Among all 54,264 six-edge sets the sheet histogram
is

```text
1:16807,  2:32417,  4:5005,  8:35.
```

Six constraints are necessary; their incidence, not their count, is the
missing coordinate.

### 3. Robust synchronization is edge connectivity

After deleting a constraint set `D`, the number of sheets is exactly
`2^(c(F-D)-1)`.  Thus survival under every deletion of fewer than `k` edges
is equivalent to `lambda(F)>=k`.  Trees are optimally sparse and maximally
fragile.  For `q>=3`, the minimum one-edge-fault-tolerant synchronizers are
the labelled `q`-cycles, counted by `(q-1)!/2`; AP13 has 360.

This is the precise survivor of the snark intuition that “bridges matter.”
It transfers edge robustness, not noncolorability, irreducibility, or prime
factorization.

### 4. Tournament unavoidability separates from local reducibility

Concurrent HYP-9080 extracts the other half of the apex architecture without
claiming a graph map.  With `S(C)=H(C)-disc(C)`, its exact deletion charge is

```text
chi_x(C)=2[S(C)-S(C-x)].
```

Thus `chi_x>=0` is a native size-reducing certificate, while the conjecture
that every tournament has such a vertex is the missing unavoidability step.
It would imply the open inequality `H>=disc` by induction.  The stronger
pointwise claim is already false for a source over `C_3`, whose source charge
is `-2`; exact checks only establish the existential and averaged proposals
through labelled order six and all order-seven isomorphism representatives.

The boundary-tensor lesson suggests a sharper target than another scalar
average: factor `chi_x` through the rooted sign and complete insertion-fiber
response, then seek a local positive atom whose sum is the global charge.
This retains precisely the data that an Euler-style scalar analogy would
discard.  Local unavoidability, its averaged strengthening, and `H>=disc`
all remain **OPEN**.

## Past snark themes: survivors and withdrawals

### Exact survivors

- THM-261: `{i,j}->e_i-e_j` identifies `K(n,2)` with the orthogonality graph
  of positive `A_(n-1)` roots.  At `n=5` this is Petersen.
- THM-4113: maximal noncrossing half-Kempe states are exactly a rooted
  combination of noncrossing partitions into blocks of sizes two and three.
  On five boundary points, cyclic planarity deletes a crossing `C_5` from
  the 15 Petersen compatibility edges and leaves ten states.
- Full cycle-length profiles and witness cycles distinguish displaced
  structure from genuine sparsity.
- Resistance, oddness, criticality, extension-support rank, and deletion
  distance are distinct obstruction coordinates worth tabulating together.

### Withdrawn transfers

- `E_ij-E_ji in so(n)` and `E_ij in sl(n)` share a pair label but are not the
  same vector.  The root-support graph forgets orientation and sign.
- Cubic three-edge-coloring is equivalent to a nowhere-zero `F_2^2` flow,
  hence a 4-flow, not a `Z/3Z` flow.  Cubic graphs cannot be two-edge-colored.
- THM-264's Ω-girth dichotomy is finite-exact only through `n=6`.  THM-343's
  order-seven `alpha_1=3`, `Omega=K_1 disjoint_union K_2` witness refutes the
  old global threshold; the global absence of intermediate girth is OPEN.
- Snark sums and dot products do not have a proved unique factorization.
  Petersen/Blanusa/flower labels therefore do not support a prime-number
  dictionary.

## Cross-frontier connection contracts

### LRC(14): synchronization needs a safety sidecar

```text
source:       AP odd-shell half-toggle packets;
target:       reciprocal-commutator phase kernel;
map:          odd-edge constraint op -> epsilon_o=epsilon_p;
preserved:    phase-sheet identity, component count, deletion robustness;
destroyed:    physical time, safe interval, clearance, owner, ancestry;
sidecar:      labelled physical orbit t*v, interval walls, owner packet;
cheap test:   star/path positive, isolated-13 cyclic hostile, t=0 unsafe.
```

The next useful object is not another unlabelled constraint graph.  It is a
**safe response tensor** indexed simultaneously by a synchronized physical
phase, the AP collar interval, and endpoint owners.  Connectivity should be
used only as the entry gate.  A saturated graph cannot make an unsafe
physical point safe; `t=0` is the decisive hostile.

### Planar Jacobian: response states modulo target shear

THM-2230 proves that for a fixed `P` possessing a constant-Jacobian mate,
the fiber of `Q -> Jac(P,Q)` is exactly the target-shear orbit
`Q+k[P]`.  THM-4103 reduces the live theta-only boundary to the degree packet
`{7,12,21}` and only five fully labelled target-infinity response profiles.
This suggests the following finite compiler:

```text
source:       shear-normalized Newton-face chart plus puncture labels;
target:       one of the five labelled infinity-response profiles;
map:          chart -> locally realizable response relation;
preserved:    Jacobian response and shared puncture labels;
destroyed:    representative Q, leading-face coefficients, valuations;
sidecar:      a proved shear normal form, cyclic face order, coefficient
              valuation, and Galois orbit labels;
cheap test:   all five THM-4103 profiles on each of the three live walls.
```

The first required lemma is a restriction/lift theorem: response relations
computed after shear normalization must lift to actual polynomial charts and
glue without changing the retained leading face.  Without it, intersecting
five finite profile sets is only a necessary sieve.  JC(2) remains OPEN.

### Four-dimensional Kakeya: bounded tangency debt

The apex discharging mechanism suggests a conditional template after
polynomial partitioning: ordinary cells contribute a conserved broad surplus,
while tangent or multiply incident tube packets consume local allowances.
To become mathematics, a Kakeya version would need

```text
global broad surplus > sum_(defect packets) local tangency allowance
```

uniformly in the partition degree and at every induction scale.  The present
obstruction is exact: unlike two or three digons, tangency packets can grow
with degree and proliferate across scales.  A four-dimensional boundary has
no cyclic noncrossing-matching carrier, so link topology, direction cap,
scale, and tube multiplicity must replace the planar rotation sidecar.

The cheapest hostile is THM-4035: the same coarse `C_60` phase admits both an
all-broad and an all-narrow finite embedding.  Any proposed boundary tensor
that cannot separate those two has already forgotten the Kakeya predicate.
No Euclidean Kakeya estimate follows from the apex theorem.

### Sixty-phase clocks and arithmetic sequences

There is an exact elementary common clock modulo ten.  Fibonacci residues
`F_n mod 10` have Pisano period 60.  Triangular residues
`T_n=n(n+1)/2 mod 10` have period 20: for a period `p`, the identity

```text
T_(n+p)-T_n=pn+p(p+1)/2
```

forces `p=0 mod 10`, while `p=10` leaves the constant residue five and
`p=20` works.  Hence the paired residue state `(F_n,T_n) mod 10` has period
`lcm(60,20)=60`.  This is a clock identity, not a conjugacy with a Sturmian
or AP tail state.  The apex/AP lesson says to retain

```text
clock position + boundary state + native transition + target predicate.
```

The exact `C_60` collision in THM-4035 shows why the period alone is blind.
A Sturmian/AP “60-phase tail law” can transfer only after its phase-state
map is specified and shown to preserve the event of interest; numerical
period agreement is not such a map.

## Precise next problems

1. **Snark boundary atlas, next layer.**  The adjacent-edge tensors for
   `K_4`, Petersen, `J_5`, and both Blanusa snarks are now exact.  Determine
   whether the normalized `k_e=m_e/2` controls critical edge deletion or
   resistance, then add six-pole cuts together with the cycle witnesses that
   meet them.  The hostile is two graphs with the same `k_e` histogram but
   different criticality.
2. **Kempe fixed-point compression.**  Starting from THM-4113's exact
   2--3-block partition grammar, compute the orbit algebra of simultaneous
   subset switches.  Determine the smallest sidecar beyond the noncrossing
   semi-matching that decides semi-D-reducibility.
3. **Fault-tolerant LRC entry.**  Classify minimum odd-shell constraint graphs
   robust under two edge failures, then pair each class with the actual AP13
   safe-interval/owner tensor.  Connectivity is the proved entry gate;
   clearance is the unsolved consumer.
4. **Jacobian face-response relation.**  For each live THM-4103 wall, build
   the five-profile local extension table before and after target shear.
   A single profile that passes all coarse degree tests but fails labelled
   face gluing would be the most informative hostile.
5. **Kakeya defect-subsidy test.**  On a finite polynomial-partition model,
   measure defect count, tangency allowance, broad surplus, and scale loss.
   Reject the compiler unless it distinguishes THM-4035's all-broad and
   all-narrow realizations of the same clock.
6. **Tournament response-charge lift.**  Expand HYP-9080's `chi_x` into the
   rooted ear/deletion boundary packet and test whether an orbit-summed local
   atom is nonnegative.  The source-over-`C_3` example is the pointwise
   hostile; an all-negative charge row would refute local unavoidability
   without by itself refuting `H>=disc`.

The common research move is now precise: replace a scalar obstruction bit by
the smallest boundary relation that composes, audit what the quotient loses,
and demand a lift theorem before importing a local certificate globally.
