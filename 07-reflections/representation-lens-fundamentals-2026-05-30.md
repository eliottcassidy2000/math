# Representation Lens: What Is Fundamental Here

**Session:** codex-2026-05-30-representation-lens  
**Status:** abstract synthesis, not canon  
**Prompt:** think from a representation-theory lens and identify what is most
fundamental across even the smallest concepts in the repo.

## Executive Thesis

The most fundamental object in the repo is not a tournament, not a Hamiltonian
path count, and not the conflict graph. It is:

```text
a sign assignment on the positive roots of type A, modulo the Weyl group,
then probed by functors that forget, split, pack, and evaluate.
```

Concretely:

```text
vertices                         basis vectors e_i
arc {i,j}                         positive root e_i - e_j
tournament                       orientation/sign choice on every positive root
vertex relabeling                 Weyl group S_n action
opposite tournament               Chevalley/complement involution
score sequence                    Cartan/standard projection
3-cycles                          first root relation/curl
Hamiltonian paths                 compatible Weyl-chamber chains
Omega(T)                          overlap algebra of odd-cycle root packets
H(T)=I(Omega,2)                   graded trace/evaluation after packing
```

The representation-theoretic core is therefore:

```text
Root sign cube -> S_n quotient -> residue/phase/incidence probes -> scalar invariants
```

The repo keeps rediscovering this same process in different dialects.

## The Base Module

Let

```text
E_n = Q-span of 2-subsets {i,j}, dimension C(n,2).
```

This is the permutation representation of `S_n` on unordered arcs. It
decomposes as:

```text
E_n = S^(n) + S^(n-1,1) + S^(n-2,2).
```

Interpretation:

```text
S^(n)        trivial component       total mass / orbit count / spine
S^(n-1,1)    standard component      score, hierarchy, Cartan shadow
S^(n-2,2)    two-row component       residual pair geometry, sea, complement bulk
```

This is the cleanest rigorous triple in the repo. Many older triples are
shadows of it: hypotenuse/vertical/horizontal, spine/ribs/sea,
fixed/boundary/free, inert/ramified/split, scalar/tournament/cooperation.

The Schur-Weyl/Burnside script confirms that the transition-orbit count splits
into exactly these three multiplicity channels. For small `n`:

```text
n   V_n      m_triv   m_standard   m_two_row   total
3   2        2        2            0           4
4   4        4        8            4           16
5   12       12       36           40          88
6   56       56       240          408         704
7   456      456      2584         5872        8912
8   6880     6880     47376        134032      188288
```

This table is more than a count. It says the smallest global decomposition of
tournament space has only three visible isotypic directions.

## The Smallest Concepts, Re-read

### Vertex

A vertex is not merely an element of `[n]`. Under the representation lens it is
a weight coordinate. Deleting a vertex is restriction from `S_n` to the
parabolic subgroup `S_{n-1}`. This is why so many deletion results look like
branching rules.

Examples:

- Claim A / OCF deletion formulas;
- projection-kill and near-kill vertices;
- good-cut/SCC defect;
- endpoint-transfer witnesses.

### Arc

An arc is a signed positive root. Flipping an arc is multiplication by the
corresponding coordinate sign. In Coxeter language it resembles a root
reflection, but in the tournament cube it is better treated as a coordinate
flip in the sign representation of the root space.

The arc is the minimal place where three structures meet:

```text
root coordinate       Lie/Weyl language
binary bit            Hamming/Krawtchouk language
comparison outcome    application/ranking language
```

### Score

The score sequence is the projection onto the standard representation. It is
the Cartan shadow of the full orientation. This explains both its power and
its limits:

- it predicts much of `H`;
- it controls `c3` through Rao-style formulas;
- it cannot see the `S^(n-2,2)` residual where many failures live.

When the repo says "score captures 90 percent but not all," representation
theory says:

```text
the standard isotypic component is large enough to control coarse behavior,
but the two-row component carries the irreducible residual.
```

### 3-Cycle

A directed 3-cycle is the first nonzero root relation:

```text
(e_i-e_j) + (e_j-e_k) + (e_k-e_i) = 0.
```

It is therefore the first obstruction to a global chamber order. This is why
3-cycles keep appearing as atoms:

- first non-transitive pattern;
- first OCF term `alpha_1`;
- first path-homology obstruction;
- first "curl" after the score gradient;
- first ramified/tournament sector in the Cartan trichotomy.

The 3-cycle is not just a small cycle. It is the minimal relation among roots.

### Odd Cycle

An odd directed cycle is a higher root packet with parity that survives the
Redei involution. OCF says Hamiltonian path count can be computed by selecting
vertex-disjoint odd-cycle packets and weighting each packet by `2`.

Representation reading:

```text
odd cycles are packets of root relations;
disjoint packets commute;
Omega records noncommutation by shared support;
I(Omega,2) is the graded trace over commuting packets.
```

This is the most promising route to a representation-refined OCF: replace the
scalar identity by a graded character identity for the module spanned by
independent odd-cycle packets.

### Hamiltonian Path

A Hamiltonian path is a Weyl chamber chain: an ordering of weights whose
adjacent root signs agree with the orientation. Counting Hamiltonian paths is
therefore counting compatible chambers.

This makes `H(T)` less mysterious:

```text
H(T) = number of Weyl chambers compatible with the sign pattern T.
```

The OCF identity is then a chamber-counting theorem expressed in the packet
algebra of odd root relations.

### Good Cut

Good cuts looked like base-path coordinates. THM-354 says they are really:

```text
goodCutCount(T,P) = n - #SCC(T).
```

Representation reading: a cut is a parabolic boundary along a chosen chamber
chain. Bad cuts are boundaries between composition factors in the condensation
poset. The base path is a coordinate choice; the SCC defect is the invariant.

This is a recurring lesson:

```text
coordinate statistic upstairs -> invariant defect downstairs.
```

### Self-Complementary And Dark Objects

Complement is the `Z_2` character attached to root-sign reversal. SC objects
are fixed or self-dual under that involution. Dark/even graph phenomena are
the sign-representation complement: the same Burnside sum split by whether an
automorphism acts trivially or nontrivially on the edge sign.

The smallest principle:

```text
light/dark = trivial/sign character decomposition.
```

This is one of the cleanest representation-theoretic bridges in the repo.

## The Four Fundamental Operations

Every recurring thread I read can be classified as one of four representation
operations.

### 1. Invariants

Take fixed vectors or trivial-isotypic projection.

Examples:

- Burnside orbit counts;
- `V_n = A000568(n)`;
- Hamming/Krawtchouk spherical averages;
- self-complementary fixed strata;
- quotient buckets and class counts.

Mathematical form:

```text
dim M^G = <chi_M, 1_G>.
```

### 2. Restriction And Branching

Forget a vertex, endpoint, or coordinate. This restricts a module from `S_n`
to a subgroup.

Examples:

- vertex deletion;
- source/sink recursion;
- endpoint transfer;
- fractal codec;
- single-core H=63 signatures;
- Paley orbit subcomplexes.

Mathematical form:

```text
Res^{S_n}_{S_{n-1}} M.
```

Most "residue" phenomena are failures of a scalar invariant to commute cleanly
with restriction.

### 3. Character/Eigenspace Splitting

Decompose into characters, Walsh levels, Fourier modes, or Cartan sectors.

Examples:

- complement kills odd Walsh levels;
- Krawtchouk coordinates;
- Paley `C_p` and QR eigenspaces;
- trace alternation;
- Cartan scalar/antisymmetric/symmetric split;
- path-homology symbol matrices.

Mathematical form:

```text
M = direct sum_chi M_chi.
```

Most "phase" phenomena are statements that one channel vanishes, flips sign,
or dominates.

### 4. Intertwiners And Incidence

Study maps between modules rather than the modules alone.

Examples:

- GLMY boundary maps;
- endpoint-transfer matrices;
- bucket transport matrices;
- collision hypergraphs;
- Smith/torsion profiles;
- private witness row-rank conditions.

Mathematical form:

```text
d : C_k -> C_{k-1},     rank(d), ker(d), coker(d).
```

Most "incidence" phenomena are failures of a support-level picture to remember
rank, parity, or torsion.

## The Representation Meaning Of Residue, Phase, Incidence

The recent residue/phase/incidence split becomes sharper:

```text
residue   = cokernel or fiber left by restriction/projection
phase     = character channel after eigenspace decomposition
incidence = rank/torsion of an intertwining map
```

So the triad is not just a useful taxonomy. It is the standard triad:

```text
objects, decompositions, maps.
```

This may be the most fundamental abstraction in the current repo.

## Paley As The Clean Test Case

Paley tournaments are special because the symmetry group is large enough to
replace brute force by character theory:

```text
C_p translations              circulant diagonalization
QR subgroup                   signed permutation symmetry
Legendre character            distinguished eigenvector
Gauss sums                    flat spectral profile
F_2 rowspace                  QR/Hamming/Golay code shadows
```

Paley is not "the best tournament" in every sense. It is the tournament where
the representation theory has the fewest moving parts. It is the lab where
phase channels become visible.

For `p=7`, the existing representation script shows:

```text
Paley orientation sigma_P = (1,1,-1)
H(sigma_P) = 189
quadratic H Hessian eigenvalues = 7, -3.5, -3.5
Paley direction is the trivial QR character channel
```

This is the small model for every later Paley/circulant question:

```text
trivial character channel       coherent/Paley direction
nontrivial character channels   oscillatory alternatives
H comparison                    competition among channels
```

The Paley-versus-interval story is therefore not primarily residue. It is
phase competition: flat QR channels versus concentrated interval channels.

## OCF As The Missing Character Theorem

The scalar OCF identity is:

```text
H(T) = I(Omega(T), 2).
```

A representation-theoretic refinement would ask for an identity before taking
dimensions:

```text
Chamber module of T
    ?= graded packet module generated by vertex-disjoint odd cycles,
       evaluated with a 2-dimensional local fiber.
```

Informally:

```text
Hamiltonian paths are chambers.
Odd-cycle packets are commuting root-relation defects.
The factor 2 is the local two-state fiber created by resolving each defect.
```

If such a module identity exists, it would explain at once:

- Redei parity;
- forbidden small H values;
- Walsh evenness;
- the 2-adic tower;
- why independent odd cycles are the right packets;
- why path homology sees exactness in low degree.

This is the best abstract target from this session.

## What The Forbidden Values Look Like In This Lens

`H=7` would require:

```text
7 = 1 + 2*3.
```

So the packet module would need exactly three rank-1 odd-cycle packets and no
higher independent packets. But the root system of a complete graph does not
permit that shape: once three such relations exist in the wrong support size,
another relation or longer odd cycle is forced.

Representation reading:

```text
the requested small character vector is not in the image of the root-sign
orientation module.
```

The H=63 single-core examples show the positive version: a parabolically
induced single-core module can realize:

```text
H = 1 + 2*31.
```

The immediate theorem target remains the single-core image gap:

```text
r_core(signature) notin {3,10}.
```

Representation translation: classify the character image of the induced
binary signature module for complete-Omega single-core tournaments.

## Path Homology As Exactness Of A Root Complex

The theorem `beta_2(T)=0` for all tournaments wants a representation proof.
The likely shape:

```text
C_3 -> C_2 -> C_1
```

is an `S_n`-equivariant chain segment built from allowed directed paths. The
claim says the relevant `C_2` homology has no surviving isotypic component.

This is stronger than a rank computation. It predicts that each irreducible
appearing in `ker d_2` also appears in `im d_3` with the right multiplicity.

Concrete next step:

```text
compute the S_n character of C_1, C_2, C_3 for the universal tournament path
complex and decompose ker/im for n <= 7.
```

If the cancellation is isotypic, the proof route is representation-theoretic.
If cancellation mixes isotypic components, the proof route is incidence/torsion.

## What Is Most Fundamental

After this pass, the deepest recurring units are:

1. **Sign on a root.**  
   The atomic choice behind every arc, bit, comparison, and flip.

2. **Support of a relation.**  
   The atomic obstruction behind cycles, Omega vertices, deletion residues,
   and path-homology chains.

3. **Group action on supports.**  
   The atomic symmetry behind isomorphism, Burnside, Paley orbits, and
   quotient transport.

4. **Character channel.**  
   The atomic phase behind Walsh parity, Krawtchouk modes, trace alternation,
   and circulant eigenspaces.

5. **Intertwiner rank.**  
   The atomic incidence fact behind boundary maps, endpoint transfer, bucket
   matrices, and Smith/torsion warnings.

Everything else is a construction from these five units.

## The Cleanest Dictionary

```text
tiny concept              representation meaning
---------------------------------------------------------------
arc                       signed positive root
flip                      coordinate sign action
vertex deletion           restriction to parabolic subgroup
score                     standard/Cartan projection
3-cycle                   minimal root relation
odd cycle                 parity-surviving root-relation packet
disjoint cycles           commuting packets
Omega                     noncommutation graph of packets
H                         compatible chamber count
OCF                       graded trace over packet algebra
good cut                  parabolic boundary
SCC defect                invariant composition-length defect
self-complement           Chevalley-involution fixedness
dark graph                sign-character sector
Krawtchouk coordinate      spherical function of Hamming scheme
Paley                     Legendre-character fixed orientation
path homology             homology of an equivariant incidence complex
residue                   cokernel/fiber after projection
phase                     eigenspace/character channel
incidence                 rank/torsion of an intertwiner
```

## What To Do Next

1. **Representation-refined OCF.**  
   Try to define the graded packet module whose dimension is
   `I(Omega(T),2)`, then ask whether Hamiltonian paths give a basis or only a
   filtered shadow.

2. **Irrep audit of path homology.**  
   For `n <= 7`, compute characters of the chain groups and boundary ranks by
   isotypic component. Test whether `beta_2=0` is componentwise exactness.

3. **Single-core signature module.**  
   Treat H=63/H=7/H=21 complete-core signatures as a parabolic induction
   problem. Prove image gaps `3` and `10` if possible.

4. **Paley character table notebook.**  
   For `p=7,11,19,23`, write the same data in one language: `C_p` characters,
   QR orbit characters, H channels, path-homology channels, and code channels.

5. **Feature taxonomy update.**  
   In `tournament_tda.py`, make `residue_*`, `phase_*`, and `incidence_*`
   explicitly correspond to restriction residues, character channels, and
   intertwiners.

## Source Trail

Local sources read or re-run:

- `07-reflections/the-grand-synthesis-s21.md`
- `07-reflections/the-triple-representation.md`
- `07-reflections/krawtchouk-coordinates-of-tournament-space.md`
- `07-reflections/everything-is-sl-n.md`
- `07-reflections/waggly-as-schurian-coherent-configuration.md`
- `07-reflections/burnside-tower-as-lie-algebra.md`
- `07-reflections/hidden-orthogonality-everywhere.md`
- `07-reflections/residue-phase-incidence-synthesis.md`
- `07-reflections/cartan-ghosts-synthesis.md`
- `07-reflections/dark-tournaments-and-the-sign-representation.md`
- `07-reflections/paley-gives-dual-codes.md`
- `04-computation/gn_schur_weyl_s217.py`
- `04-computation/deep_synthesis_representations.py`
- `04-computation/representation_chi_bridge.py`

