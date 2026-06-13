# Residue, Phase, And Incidence: Three Ways Threads Recur

**Session:** opus-2026-05-30-S4
**Mode:** exploratory formalization after random-thread sampling
**Status:** synthesis, not canon
**Related:** HYP-1779, HYP-1780, HYP-1785, HYP-1795, HYP-1796, HYP-1797

The recent repo spine says:

```text
choose supports -> project/forget -> measure the survivor.
```

That is real, but a deliberately noisy file sample shows it is not the only
primitive.  Some threads are best described as residues after projection.
Others are better described as phase or orthogonality phenomena.  A third
family lives in incidence matrices where neither scalar residues nor adjacency
graphs remember the relation.

The useful split is:

```text
residue     = what survives a projection
phase       = what survives an orthogonal decomposition
incidence   = what survives the operation interface
```

These three are not competing theories.  They are three failure modes of
compression.

## Sampled Threads

The intentionally mixed sample included:

- good-cut/SCC residue and bucket gaps;
- H=63 exact kill versus THM-025 near-kill;
- Lonely Runner quotient gaps;
- `F_k(T) mod 2` tournament independence;
- trace alternation between Paley and interval circulants;
- super-orthogonality / Walsh parity / OCF;
- endpoint-transfer private witnesses and collision hypergraphs;
- Petersen and Coxeter side quests;
- lex/Fourier/information synthesis;
- older "what each piece represents" notes around completeness and the
  parabolic frustration law.

The sample is useful because several of these threads should not be forced
into the same vocabulary too early.

## Primitive 1: Residue

Residue begins with a quotient or projection:

```text
pi : X -> Y.
```

The visible object is a shadow in `Y`.  The information is the part of the
fiber, complement, or boundary that is not determined by the shadow.

Examples:

- THM-354: good-cut count is not a base-path statistic but the SCC defect
  `n - #SCC(T)`.
- THM-355: empty quotient fibers are zero transport rows and columns.
- H=63: deleting the core vertex kills every odd cycle, so the residue is
  empty and rigid.
- THM-025: deleting the high-loss vertex leaves only two old cycles, but they
  are disjoint, so the tiny residue has rank 2 and remains dangerous.
- Paley versus interval T7: the odd-cycle support shadow can be shared while
  multiplicity and disjointness fibers differ.
- Lonely Runner: forbidden arcs pull back to an interval cover, and a witness
  is either a positive quotient gap or a boundary residue.

Formal search rule:

```text
If an invariant is pure on a quotient but was defined upstairs,
find the defect downstairs that it is secretly counting.
```

Residue is the right language for support loss, missing buckets, boundary
witnesses, deletion profiles, and fiber multiplicity.

## Primitive 2: Phase

Phase begins with a decomposition into characters, eigenspaces, or orthogonal
channels.  The important event is not that a support is forgotten.  It is that
some channel vanishes or changes sign.

Examples:

- Complement symmetry `H(T)=H(T^op)` forces odd Walsh coefficients to vanish.
- THM-094 says `F(T,x)=(1+x)^(n-1) mod 2`; tournament-dependent channels
  disappear modulo 2.
- Super-orthogonality links complement parity, Walsh levels, OCF amplitudes,
  and path-homology balances.
- THM-136 trace alternation is a phase statement: Paley and interval
  eigenvalue phases approach `pi/2` from opposite sides, and the sign pattern
  is controlled by `k mod 4`.
- Circulant Paley/interval crossover is not a near-kill.  It is a competition
  among phase-aligned trace or Walsh channels.

Formal search rule:

```text
If a pattern is periodic in degree, parity, k mod 4, or character orbit,
do not first ask what projection forgot.  Ask which channel is orthogonal.
```

Phase is the right language for Walsh degree, mod-p Taylor zeros, trace signs,
Krawtchouk coordinates, circulant eigenspaces, and the Paley/interval
arithmetic dichotomy.

## Primitive 3: Incidence

Incidence begins with an operation:

```text
parent object --operation--> child object.
```

The carrier is a matrix or hypergraph.  The mistake is to replace that matrix
too early by a scalar invariant or an ordinary adjacency graph.

Examples:

- THM-356: private odd child witnesses imply full endpoint-transfer row rank,
  but support matching alone does not imply rank over `F_2`.
- Merged endpoint transfer loses private witnesses before it loses rank.
- Even-graph endpoint transfer can have full support matching and still fail
  rank because mod-2 cancellation occurs.
- Endpoint SC collision triples are not necessarily triangles in the parent
  merged metagraph.  They are 3-uniform transfer-incidence hyperedges.

Formal search rule:

```text
When a support shadow has the right shape but rank fails,
the missing structure is incidence parity or torsion, not adjacency.
```

Incidence is the right language for endpoint recursion, transfer matrices,
Smith factors, private pivots, support matchings, and collision hypergraphs.

## What The Three Primitives Represent

The older "what each piece represents" reflection argued that tournaments are
complete binary decision systems and that the parabolic law

```text
E[Delta c3 | score=s] = s(n-1-s)/2
```

is completeness made quantitative.  The three primitives can be read as three
ways completeness resists compression:

```text
residue   : completeness creates forced survivors after forgetting;
phase     : completeness creates forced vanishings after orthogonal splitting;
incidence : completeness creates forced parities across operations.
```

This gives a cleaner feedback loop:

1. Start with a scalar pattern or failed conjecture.
2. Ask whether it is about a projection, a phase channel, or an operation
   interface.
3. Add the missing data type: residue vector, phase vector, or incidence
   matrix.
4. Re-test the examples that originally looked unrelated.

## Thread Comparisons

| Thread | Best primitive | What is really happening |
|---|---|---|
| Good-cut count | residue | SCC defect survives base-path projection |
| Bucket gap `1` | residue | empty fiber creates zero transport row |
| H=63 classes | residue | exact kill leaves empty odd-cycle survivor |
| THM-025 | residue + phase | tiny deletion residue has enough disjointness to break root shape |
| Lonely Runner | residue | forbidden interval cover leaves boundary/gap residue |
| `F_k mod 2` | phase | tournament-dependent Taylor channels vanish |
| Paley/interval trace alternation | phase | eigenvalue phases straddle `pi/2` |
| Super-orthogonality | phase | multiple decompositions force compatible vanishings |
| Endpoint transfer | incidence | parity injectivity lives in matrix rank |
| SC collision triples | incidence | ternary collision is hyperedge incidence, not metagraph triangle |
| Petersen/Coxeter side quests | phase/incidence | orthogonality graphs and root-incidence shadows, not direct residues |

## New Hypotheses

### HYP-1795: Residue-Phase Bifurcation

The next useful tournament features should be separated into a residue vector
and a phase vector.  Real-root failures and H-spectrum gaps should lean on
residue rank; Paley/interval and circulant maximizer questions should lean on
phase-channel dominance.  Mixed phenomena, such as THM-025, should require
both a nontrivial residue and a bad phase/root channel.

### HYP-1796: Incidence-Layer Necessity

Whenever a quotient support graph has the expected matching or adjacency shape
but the `F_2` rank statement fails, the missing object is an incidence layer:
a transfer matrix, collision hypergraph, or Smith torsion profile.  Endpoint
transfer is the model case.

### HYP-1797: Completeness-Defect Principle

The parabolic law is the zeroth-order checksum of tournament completeness.
For generalized or applied settings with missing, tied, weighted, or noisy
pairwise comparisons, the first feature should be a completeness-defect vector:
how much the observed structure deviates from the forced tournament parabola.
This should predict when tournament-native invariants such as `H`, residue
rank, and phase vectors remain stable.

## Practical Next Steps

1. Extend `tournament_tda.py` into three blocks:

```text
residue_*   deletion, SCC, support, fiber, gap features
phase_*     Walsh/Krawtchouk, trace, complement-character features
incidence_* endpoint/transport rank and torsion summaries
```

2. Re-run the known contrast set:

```text
transitive, Paley, interval, H=63, THM-025, n=6 H=37 trap,
beta_3 anomalies, endpoint-transfer even-graph rows.
```

3. For any new conjecture, classify it first as residue, phase, incidence, or
mixed.  This should prevent the recurring mistake of proving a scalar shadow
when the theorem lives one layer higher.
