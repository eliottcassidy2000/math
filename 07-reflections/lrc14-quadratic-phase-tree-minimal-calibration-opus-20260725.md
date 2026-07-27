---
source: codex-2026-07-25-quadratic-phase-tree-calibration
status: >
  PROVED ALGEBRAIC REFINEMENT + VERIFIED-EXACT C_13 CENSUS +
  INDEPENDENTLY HOSTILE-AUDITED. Complete untwisted pair-union energies need
  only the total singleton energy, not every singleton separately, to
  decide grouped-current vanishing. For p>=3 pair-twist orbits, a
  connected probe graph and one anchored singleton magnitude reconstruct
  every component current up to common phase; a connected nonbipartite
  graph needs no anchor. A single global cyclic twist bank has the same
  power when its unique-difference graph is connected. In C_13 this holds
  for every two-support, exactly the 208 non-AP three-supports, and 364
  four-supports, but for no support of size at least five. These are
  algebraic phase-retrieval statements. Current LRC theorems do not yet
  supply the required component-labelled lawful twists, so no LRC(14)
  profile is excluded.
depends_on:
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2355-quadratic-grouped-current-repair
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2312-sparse-root-bispectrum-positive-word-current
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
  - lrc14-two-frequency-orientation-current-leakage-opus-20260725
script: 04-computation/lrc14_quadratic_phase_tree_minimal_calibration_probe.py
output: 05-knowledge/results/lrc14_quadratic_phase_tree_minimal_calibration_probe.out
script_sha256: 1cdeed147301dbc28ff416fe1d91af44f7ac376259bd98a362e8958056d09133
output_sha256: 239aad210fe947c5592802a61ff204017a0042956069acd5e069769a3f03f381
hash_basis: working-tree bytes (LF)
---

# One calibration is enough for a quadratic phase tree

THM-2303 identifies the missing LRC coordinate as relative phase between
terminal-component currents. THM-2355 gives the exact quadratic repair:
singleton component energies plus complete pair-union energies decide
whether the grouped current vanishes, while lawful pair twists recover the
complex phase transports.

This note makes that repair measurement-economic and connects it to the
global shifted-frequency banks already present in the LRC route. The
algebraic gain is unconditional. Its realization by genuine LRC operations
remains open.

## 1. Complete pair unions need only one diagonal total

Let

```text
z_1,...,z_r in C,

A=sum_i |z_i|^2,

E_ij=|z_i+z_j|^2.                                  (1)
```

Then

```text
sum_(i<j) E_ij
 =(r-1)A+2 Re sum_(i<j) z_i conjugate(z_j),         (2)

|sum_i z_i|^2
 =A+2 Re sum_(i<j) z_i conjugate(z_j).              (3)
```

Subtracting gives the exact compression

```text
|sum_i z_i|^2
 =sum_(i<j) E_ij-(r-2)A.                            (4)
```

Thus THM-2355's complete-pair criterion does not require the singleton
energies separately. Their one total `A` is enough. Formula (4) is an
identity, so it detects vanishing with no positivity margin and no genericity
hypothesis.

## 2. Pair twists and the signless-incidence rank

Fix a primitive `p`-th root `omega`, with `p>=3`. For an oriented probe edge
`e->f`, put

```text
Q_ef(q)=|z_e+omega^q z_f|^2,          q in Z/pZ.    (5)
```

With normalized Fourier transform

```text
Qhat_ef(k)
 =(1/p)sum_q Q_ef(q)omega^(-kq),                    (6)
```

orthogonality and `p>=3` give

```text
Qhat_ef(0)=|z_e|^2+|z_f|^2,

Qhat_ef(1)=conjugate(z_e)z_f.                       (7)
```

The first mode is exactly the relative phase transport requested by
THM-2303.

Let `Gamma` be a connected probe graph on the nonzero component support.
Choose one root component and fix its irrelevant common phase. If one
anchored magnitude `|z_root|` is known, the second equation in (7)
recovers every neighboring current and induction along any spanning tree
recovers every `z_e`. Hence:

> **One-anchor phase-tree lemma.** A connected lawful pair-twist graph and
> one anchored singleton magnitude reconstruct all nonzero component
> currents up to a common phase. In particular, they decide whether their
> grouped sum vanishes.

There is an exact anchor-free alternative. Put `x_e=|z_e|^2`. The zero modes
in (7) form the signless-incidence system

```text
x_e+x_f=Qhat_ef(0),                  {e,f} in E(Gamma).  (8)
```

For a connected graph, the kernel of this system is zero exactly when
`Gamma` is nonbipartite. If `Gamma` is bipartite with parts `L,R`, its
kernel is the line which is `+t` on `L` and `-t` on `R`. This follows by
propagating `x_f=-x_e` along paths; an odd cycle forces `x_e=-x_e`, while
on a bipartite graph the alternating assignment survives.

Therefore a connected nonbipartite twist graph reconstructs all magnitudes
from (8), then all phases from (7), with no singleton input. A connected
bipartite graph needs one scalar anchor.

The anchor is genuinely necessary. On the labelled path `1-2-3`,

```text
(-1,2,-1)       and       (-2,1,-2)                 (9)
```

give on both edges

```text
Qhat(0)=5,             Qhat(1)=-2,                  (10)
```

so their complete pair-twist responses agree, but their grouped totals are
respectively `0` and `-3`. This is the repaired missing-singleton hostile
from THM-2355, now identified as the one-dimensional bipartite kernel.

## 3. One global twist bank can replace pairwise probes

Pair-specific masks may be unavailable in LRC. A global cyclic action can
nevertheless expose the same phase tree.

Give each component a distinct label `a_e in Z/pZ` and define

```text
H(q)=|sum_e omega^(q a_e) z_e|^2.                   (11)
```

Its normalized Fourier coefficients are

```text
Hhat(d)
 =sum_(a_e-a_f=d) z_e conjugate(z_f).               (12)
```

Call `{e,f}` a **unique-difference edge** when the ordered difference
`a_e-a_f` is realized by no other ordered component pair. On such an edge,
(12) is the individual phase transport

```text
Hhat(a_e-a_f)=z_e conjugate(z_f).                   (13)
```

Consequently:

> **Global-twist phase-tree lemma.** If the unique-difference graph of the
> labelled nonzero support is connected, one global twist bank (11), plus
> one anchored component magnitude, reconstructs every current up to common
> phase. If that graph contains an odd cycle, the edge-product magnitudes
> determine all component magnitudes and the anchor is unnecessary.

For the last sentence, take absolute squares in (13):

```text
|Hhat(a_e-a_f)|^2=x_e x_f.                          (14)
```

After logarithms on the positive support, (14) has the same signless
incidence matrix as (8); equivalently, multiply and divide successively
around an odd cycle.

This is the exact global surrogate for THM-2355's pair-specific twists.
Its failure coordinate is also exact: collisions of equal label differences
aggregate several products in (12). A full-support CAZAC is the extreme
case, while the pending THM-2353 hostile realizes the same problem as a
target-residue collision multiplicity.

## 4. Exact `C_13` boundary

The companion exhausts every subset of `Z/13Z`. The number whose
unique-difference graph is connected is

```text
support size 1:       13 / 13
support size 2:       78 / 78
support size 3:      208 / 286
support size 4:      364 / 715
support size >=5:      0.                            (15)
```

The `78` bad triples are exactly the cyclic three-term arithmetic
progressions. For a non-AP triple all three pair differences are unique,
so the graph is a triangle. For an AP triple only the outer pair is unique
and the midpoint is isolated.

At size four, the connected cases split into `312` four-edge graphs and
`52` six-edge graphs. At size five the maximum number of unique-difference
edges is three, already too few to connect five vertices. The complete
histogram in the stored transcript records every higher support.

This is simultaneously encouraging and restrictive for LRC. THM-2305's
root fibres have support at most two, so their *root-address* twist is in
the favorable line of (15). But the desired terminal current is a sum over
continuous base components. Existing rooted energy has type

```text
integral |M_k(y)|^2 dy,                              (16)
```

whereas THM-2355 requires component-current products of type

```text
conjugate(integral_(E_e) e(-nx)dx)
          integral_(E_f) e(-nx)dx.                  (17)
```

For disjoint terminal components, the fibrewise cross term in (16) can be
zero while (17) changes under relative translation. This is exactly the
THM-2303 `F_0/F_1` hostile. Root support two therefore does not by itself
give a terminal phase-tree edge.

## 5. Connection contract and decisive test

```text
source:
  THM-2305's exact terminal word and THM-2303's component currents;

target:
  the grouped word/full current at the prescribed shell frequency;

map:
  either take lawful pair-twist DFTs edge by edge, or take one global
  affine twist DFT and select its unique-difference edges;

preserved algebraically:
  component label, complex cross product, grouped-current vanishing, and
  common-phase gauge;

destroyed by current rooted energies:
  cross-component continuous base phase;

destroyed by a coarse global twist:
  individual pairs inside one residue-difference collision class;

minimum calibration:
  one singleton magnitude on a connected bipartite phase tree, or no
  singleton on a connected nonbipartite probe graph;

needed LRC sidecar:
  prove that the twist operation preserves owner, exact word, time, and
  primitive shell colour, and that some resulting unique-difference
  forest connects the nonzero terminal components;

cheapest decisive test:
  on the actual THM-2349/2305 marked word, build the endpoint-address
  difference multigraph at successive thirteen-adic depths. Stop at the
  first depth whose unique-difference edges connect the active component
  support, and then check whether the corresponding affine frequency
  shifts remain inside the forced-current bank.                       (18)
```

If the graph never connects, the collision classes are the exact ancestry
state that must be retained. THM-2352 warns that no fixed depth is uniform;
a coefficient-dependent stopping depth remains compatible with its dense
finite-cylinder hostile.

Current canon supplies neither componentwise phase masks nor the required
global unique-difference forest with shell legality. Thus (18) is a sharp
service target, not an LRC proof, and no scalar profile is removed.

## 6. Exact verification

Run

```text
python3 04-computation/lrc14_quadratic_phase_tree_minimal_calibration_probe.py
python3 -O 04-computation/lrc14_quadratic_phase_tree_minimal_calibration_probe.py
```

The companion checks (4) on `2,400` deterministic Gaussian-integer packets;
checks the hostile (9)--(10); verifies the signless-incidence rank criterion
on all `27,475` connected labelled graphs through six vertices; exhausts all
`8,191` subsets of `C_13`; and checks `2,964` unique-edge coefficient
selectors. Normal, optimized, and stored transcripts must agree exactly.

An independent proof/code audit rederived (4), both Fourier modes in (7),
the signless-incidence kernel, the one-anchor and odd-cycle reconstructions,
and the global unique-difference selector. It independently replayed normal
and optimized Python against the stored transcript, reproduced the complete
`C_13` census, and confirmed the fixed-frequency/integrate-after-square scope
boundary.

PROVED THM-2356 (promoted at `59c933aae`) supplies a different, stronger
measurement bank: quadratic chirps over `F_(p^d)` eliminate every difference
collision by pairing the target address with a Bockstein jet. That theorem does
not supersede the sparse criterion here. The present global *linear* twist is
cheaper exactly when a difference-Sidon tree already exists; the chirp route
is designed for the collision classes where it does not.

## 7. Independent hostile audit and promoted form of THM-2356

The early algebraic candidate at `79f63619f` passed an independent
sign, normalization, boundary, and finite-field audit. With unnormalized
characters and

```text
F_z(eta,chi)=sum_x z(x) eta(phi(x)) chi(x),
```

expanding the intensity and averaging first in `chi` forces `x-u=h`.
Averaging in `eta` then forces

```text
D_h phi(u)=D_h phi(y).
```

For `h!=0`, planarity gives `u=y`; both character sums contribute `q`,
so the exact normalization is `q^(-2)` and the signs in formula (5) are
correct. At `h=0`, however, every diagonal term survives and the expression
is the total norm, independent of `y`. The off-diagonal qualifier is
load-bearing.

The specialization

```text
K=F_13[theta]/(theta^2-2),

Tr(u+v theta)=2u,

phi(x)=x^2/2
```

also passes. The element `2` is a nonsquare modulo `13`, and

```text
D_h phi(y)=h y+h^2/2
```

is a permutation for each of the `168` nonzero `h`. The off-diagonal Gram
matrix recovers its diagonal without extra data only when the support has
at least three nonzero points:

```text
|z(x)|^2=Gamma_(x,y) Gamma_(w,x) / Gamma_(w,y).
```

Singleton location and the fixed-two-site `(2,3)<->(3,2)` magnitude swap
are exact uniform obstructions, so labelled singleton energies are the
correct sharp sidecar.

The first scope draft was too pessimistic. PROVED THM-2337
(`THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar`)
already supplies the joint Abel array

```text
A(q,z),          (q,z) in G x B,
```

with all linear target/jet character transforms. After chosen additive
identifications `G,B~=K=F_169`, fix a nontrivial additive character
`Psi_K` of `K`. The `169` resulting graphs

```text
Gamma_c={ (q,z): z=q^2/2+c }
```

partition `G x B`. The joint table is canonical; the graphs depend on the
chosen field identifications. Thus

```text
Z_c(q)=A(q,q^2/2+c)
```

is a lawful **coefficient grouping**. If the joint Fourier convention is

```text
Ahat(beta,alpha)
 =sum_(q,z) A(q,z) Psi_K(beta q+alpha z),
```

then its graph-chirp amplitude is the derived finite functional

```text
F_c(a,b)
 =1/169^2 sum_(alpha,beta) Ahat(beta,alpha) Psi_K(-alpha c)
    sum_q Psi_K((b-beta)q+(a-alpha)q^2/2).          (21)
```

The denominator, both minus signs, and the order `(beta,alpha)` follow
from two-dimensional Fourier inversion. Equivalently, on graph `c` the
linear character `Psi_K(bq+az)` restricts to

```text
Psi_K(ac) Psi_K(bq+a q^2/2).
```

The constant factor disappears from intensity. Because the graphs
partition the joint array and THM-2337 proves a nonzero joint Abel limit,
some graph signal is nonzero.

This is an algebraic grouping, not an independently supplied physical
pair-probe measurement. THM-2337 gives every coefficient by finite
inversion, so `F_c` and `|F_c|^2` are exact derived linear and quadratic
functionals. It does not say that the scalar LRC dynamics separately
measure that intensity, nor provide a positive lower bound for it. The
promoted theorem repairs its older “no planar section” scope without
overstating a new physical observable.

The exact residual is **graph-singleton localization**. A nonzero graph
signal of support at least two already lands at some `q!=0`; support at
least three is fully reconstructed from off-diagonal Gram data. But a
singleton graph signal has the same complete chirp-intensity table at
every target location. THM-2337 does not yet force an aligned singleton
energy away from `q=0`, and derived chirp intensities cannot manufacture
that location information graphwise. Coherent cross-frequency
interference can locate it, so the precise next sidecar is a
phase-sensitive target anchor or aligned singleton energy. Thus the
promoted theorem provides a coefficient-level tomography identity but
still removes no LRC profile.

The singleton boundary has an exact table-level normal form. If every
chosen graph signal `Z_c` is at most one-sparse, then either one of its
singletons lies at `q!=0`, which is the desired target landing, or

```text
A(q,z)=1_(q=0) B(z).                                (22)
```

Indeed, under failure of target landing every joint coefficient with
`q!=0` is zero, and at `q=0` graph `c` reads exactly `A(0,c)`.
Conversely (22) gives

```text
Z_c(q)=B(c)1_(q=0)
```

on every graph. This vertical tensor is sharp as a formal coefficient
hostile: tensor THM-2333's full-term target convolution `delta_0(q)` with
any nonzero jet profile `B(z)`. Taking `B=1` makes every graph survive but
keeps every survivor at the zero target. It also saturates the finite
Fourier uncertainty product; over `F_13^2`,

```text
delta_0(x,y)=(1-x^12)(1-y^12)
```

is the `13 x 13` interpolation-footprint extremizer, not a claim about its
monomial-support size. This is a group-algebra boundary,
not yet a canonical interval-weight LRC row.

The promoted theorem has an exact scalar detector for this residual. On the
unchirped target line write

```text
F_c(b)=sum_q Z_c(q) Psi_K(bq)
```

and put

```text
D_c=1/169 sum_b |F_c(b)-F_c(0)|^2.                 (23)
```

Character orthogonality gives

```text
D_c
 =sum_(q!=0)|Z_c(q)|^2
   +|sum_(q!=0)Z_c(q)|^2.                          (24)
```

Thus `D_c=0` exactly when graph `c` is supported at the zero target. This
does not follow from the separate intensities `|F_c(b)|^2`: the difference
intensity retains the phase ratio between `F_c(b)` and `F_c(0)`. At the
coefficient level (23) is another exact derived quadratic functional of
THM-2337's table. As a physical LRC service it names the missing operation
precisely: one lawful graph-channel pair twist or phase-coherent
cross-frequency anchor. It is not a profile exclusion by itself.

An independent exact companion is

```text
04-computation/thm2356_planar_chirp_independent_audit.py
```

with stored output

```text
05-knowledge/results/thm2356_planar_chirp_independent_audit.out.
```

It checks `68` prime-field off-diagonal selector cells, `134,003`
prime-field graph-restriction selector cells (including the signs and
`1/169^2` analogue of (21)), the `h=0` failure, all `168` `F_169` planar
derivatives, all `168` nondegenerate trace pairings, the full `28,561`-cell
graph partition, the diagonal hostiles, `35` exact Gaussian-integer
instances of (24), `373` graph-singleton/vertical-tensor cells, and `80`
joint-to-word energy domination controls including the sharp constant.
Normal, optimized, and stored transcripts agree after LF normalization.
The working-tree SHA-256 hashes are respectively

```text
04107d2d41dd899efd3cbd963400c921a87f6afa8ea39530092ea53a5f44f6ef

92ffc132e1827cef8fa2714e64e03be4a69ff5e5578a5382311eeb9f4a0f319f.
```

One nonfatal verification-coverage caveat in the original `79f63619f`
script is that its two-site hostile constructs `first` and `second` from literally
identical hand-written cyclotomic terms, so that assertion is tautological
rather than an independent execution of the two swapped signals. The
underlying hostile is nevertheless proved by the displayed norm/cross-product
calculation and independently checked here. This is a coverage note about the
early script, not a defect in the promoted statement.

The promotion-safe LRC composition wording is:

> After chosen additive identifications `G,B~=F_169`, THM-2337's canonical
> joint coefficient table partitions into the resulting graphs
> `z=q^2/2+c`. Graph restriction is a lawful finite coefficient grouping,
> and formula (21) expresses each quadratic-chirp amplitude with the exact
> `1/169^2` scalar and negative Fourier phases. Its squared modulus is a
> derived quadratic functional, not an independently supplied physical
> pair-probe observable. The graph bank is nonzero, but a one-sparse
> survivor may lie at `q=0`; failure of all nonzero-target landings has the
> sharp vertical form `A(q,z)=delta_0(q)B(z)`. Hence no scalar profile is
> removed without an aligned singleton energy or phase-coherent
> cross-frequency target anchor.
