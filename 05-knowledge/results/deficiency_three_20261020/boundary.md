# Fifteen probes certify the junction target; golden compression has a physical collision

**PROVED, using the explicit finite cover certified below:** fourteen zero
closure counts plus one normalization determine the full extension tensor
of a physical cubic seven-pole when the intended target is one junction
and two wires. **PROVED:** Kempe uniformity, a positive near-perfect-matching
formula, and the aggregate closure inequality. **FINITE-EXACT:** the cover,
the displayed physical controls, and the golden Gram ranks and witnesses.
The general rank-84 blind space remains seven dimensional. There is no
general rooted-minor conclusion or Collatz consequence.

## 1. Inheritance, types, and the live board

The nearest proved mechanism is
[seven-port junction tomography](../blue_quotient_20261020b/junctions.md):
the 105 junction probes have rank 84 on 91 boundary-colour orbits, and
their seven missing coordinates count particular Kempe paths. The
[105/17 incidence note](../black_descent_20261020c/105_17.md) supplies the
separate eight-port matching carrier and explains why another occurrence
of 105, or rank 91, does not identify these spaces. The
[connected realization](../blue_quotient_20261020b/snarks.md) supplies a
17-vertex pole with tensor \(4J\) and additional rooted-minor information.

The canonical hostile is an ambient kernel vector that has no declared
physical realization. The correction is to retain nonnegativity, actual
Kempe moves, and port labels before drawing graph conclusions. The
least-used sidecar here is the *singleton port and the parity of the three
remaining two-colour paths*. Its near-perfect-matching description below
gives a second exact computational engine.

The wildcard is the earlier
[ternary ports / golden defect](../ternary_ports_golden_defect_20261018j.md).
Its three-terminal elimination coefficient contains
\(q^2-3q+1\), which vanishes at \(q=\phi^2,\phi^{-2}\). This is a
parameter-dependent identity, not an identity of ordinary three-colour
tensors. We test an actual golden diagram quotient below, retaining its
planarity and loop parameter.

| Lane | Source and map | Preserved predicate; loss or sidecar |
|---|---|---|
| Anchor | Physical seven-pole to 14 zero probes and one count | Exact \(mJ\); loses internal graph and rooted minor paths |
| Niche | Singleton-colour Kempe paths to near-perfect matchings | All seven missing counts; retains path parity and cycle weight |
| Niche | Sum over all path-tree closures | Detects simultaneous vanishing of the seven counts; loses their locations |
| Wildcard | Seven defect counts to a golden moment | Zero versus nonzero survives; individual counts do not |
| Wildcard | Planar matchings to the golden loop Gram form | Planar quotient only; crossing probes and a new loop parameter change the kernel |

No tournament black-degree deficit, symmetry slack, and boundary-kernel
dimension are identified in this note. In particular, the separate
\(D(z)=z^4\) center result has its own objects and proof.

## 2. The physical target needs only fifteen specified measurements

Let \(P\) be a finite loopless cubic multipole with seven distinct labelled ports,
allowing different ports at the same internal vertex and isolated free
wires. Let \(f_P(s)\) count proper edge colourings extending the boundary
word \(s\in\{0,1,2\}^7\). Closed vertex-free wire circles, if created by
closure or already present as components, have weight three. Counts are
nonnegative integers.

The parity lemma restricts support to the 546 words in which all three
colour multiplicities are odd. Global colour permutations partition them
into 91 free orbits. Let
\[
 J(s)=[s_0,s_1,s_2\text{ all distinct}]\,[s_3=s_4]\,[s_5=s_6].
\]
It supports 54 words, or nine orbits.

A closure label \((abc;de,fg)\) means one all-distinct junction on its
triple and equality wires on its pairs. Its measurement is
\[
 C_L(P)=\sum_s f_P(s)J_L(s).
\]
These are labelled-colour counts, not orbit-normalized inner products.

**Finite cover certificate.** The following fourteen closures have support
disjoint from \(J\), and their supports together contain every admissible
word outside \(J\):

~~~text
(034;12,56)  (034;15,26)  (034;16,25)  (056;12,34)
(134;02,56)  (134;05,26)  (134;06,25)  (156;02,34)
(234;01,56)  (256;01,34)  (256;03,14)  (256;04,13)
(345;01,26)  (346;01,25)
~~~

The checker verifies this statement directly on all \(3^7=2187\) words,
with the parity filter explicit, and checks its transport to all 105
target labels. There are 492 forbidden admissible words, or 82 orbits.
No minimality assertion about fourteen is needed or made.

Nonnegativity therefore gives the exact implication
\[
 C_L(P)=0\text{ for these fourteen }L
 \quad\Longleftrightarrow\quad \operatorname{supp}f_P\subseteq
 \operatorname{supp}J.                                      \tag{1}
\]
This support certificate does not use rank 84.

**Kempe uniformity lemma.** If a physical cubic seven-pole satisfies the
right side of (1), then \(f_P=mJ\) for an integer \(m\geq0\).

Proof. If it has no colouring, take \(m=0\). Otherwise fix a colouring
and one of the two target wire pairs, whose common colour is \(c\).
For another colour \(d\), the \(c/d\) subgraph consists of cycles and
boundary paths. Follow the path from one port of the selected wire.
Its other endpoint must be the other port of that same wire. If it
ended anywhere else, switching \(c,d\) on the path would change just one
port of this wire pair. That proper colouring would lie outside \(J\),
contrary to the global support hypothesis.

The switch thus changes only this wire's colour from \(c\) to \(d\),
leaving the junction palette and the other wire fixed. It is bijective:
the underlying two-colour component is unchanged, and switching again
reverses the operation. Consequently all nine assignments of colours
to the two wires have the same extension count for a fixed ordered
junction palette. Global \(S_3\) colour symmetry gives the same count
for all six ordered distinct junction palettes. This proves the claim.
The argument includes a free-wire component with no internal vertex. \(\square\)

**Target certification theorem.** The fourteen zero measurements in (1)
and the target self-closure measurement \(C_J(P)=54m\) are equivalent
to \(f_P=mJ\), for physical cubic seven-poles.

On the 91 orbit coordinates the normalization is
\(\langle f_P,J\rangle=9m\); the factor six must not be omitted.
These fifteen specified tests can certify a context-independent
replacement tensor. They do not provide disjoint minor branch sets,
simplicity, bridgelessness, or a general Petersen reduction. Nor is
their evaluation asserted to take polynomial time for arbitrary poles.

### Why this does not contradict the seven-dimensional kernel

Arbitrary tensor reconstruction and certification of one nonnegative
target are different tasks. There is also the following stronger
target-specific statement with a larger probe bank:

**Nonnegative uniqueness.** If an arbitrary nonnegative invariant parity
tensor \(f\) has all 105 measurements equal to those of \(mJ\), \(m\geq0\),
then \(f=mJ\), without assuming physical realizability.

By inherited rank 84, write \(f-mJ=\sum_i c_iw_i\).
On a \((5,1,1)\) orbit with singleton ports \(i,j\) its value is
\(c_i+c_j\); on a \((3,3,1)\) orbit with singleton \(i\), it is \(-c_i\).
For every \(i\) there is a latter orbit outside \(J\), so \(c_i\leq0\).
A \((5,1,1)\) orbit is supported by \(J\) exactly when its singleton
pair lies within \(\{0,1,2\}\). Every port belongs to an unsupported
pair. Nonnegativity on such pairs gives \(c_i+c_j\geq0\), forcing both
coefficients zero. All coefficients vanish. This proof uses all 105
measurements; the fifteen-probe theorem instead uses the physical
Kempe lemma after its support test.

The physical assumption is necessary for those fifteen probes. Put
weight nine on just one supported colour orbit of \(J\), and zero on
the others. This nonnegative invariant tensor has all fourteen zero
readouts and self-closure count 54, exactly as \(J\) does, but is not
\(J\). The checker retains this nonphysical positive hostile.

## 3. A positive matching formula for the missing seven coordinates

Use the inherited orbit vectors
\[
 w_i=+1\text{ on }(5,1,1)\text{ with singleton }i,\qquad
 w_i=-1\text{ on }(3,3,1)\text{ with singleton }i,
\]
and zero otherwise. The inherited exact identity is
\[
 6\langle f_P,w_i\rangle=-N_i,                               \tag{2}
\]
where \(N_i\) counts full colourings whose port-\(i\) colour is unique
on the boundary and whose other-colour paths all have differently
coloured endpoints. In particular \(N_i\geq0\).

Suppose now every port is incident to an internal vertex and there are
no vertex-free wire components, open or closed. Write \(v_i\) for the
incident vertex of port \(i\). Let \(M\) run over perfect matchings of
the internal graph with \(v_i\) deleted. Retain individual edge choices
if the graph has parallel edges. Remove \(M\) and port \(i\) from the
pole. The residual graph has degree two at every internal vertex and
six boundary ends, hence three paths and some cycles. Call \(M\)
admissible if every cycle has even length and each of the three paths
contains an odd number of internal vertices. If \(c(M)\) counts its
cycles, then
\[
                 N_i=24\sum_{M\text{ admissible}}2^{c(M)}.    \tag{3}
\]

Proof. Fix the unique boundary colour at port \(i\). Its internal edges
form exactly a perfect matching on the other internal vertices. The
remaining two colours must alternate on the residual components.
Closed cycles are colourable precisely when even, each with two choices.
Each path has two choices, and its endpoint edges differ precisely when
the number of internal vertices is odd. There are three choices of
unique colour and \(2^3\) choices for the three paths. This gives (3),
with no cancellation. \(\square\)

In particular \(24\mid N_i\). Formula (3) is a genuine dimer connection,
but the fixed monomer, path parity, and cycle weight are essential.
It is not the unconstrained count of all matchings of a path. If
vertex-free closed circles are added, each multiplies the answer by
three; it must not be treated as a residual two-colour cycle of weight two.

### One aggregate inequality detects the entire equality face

Let \(J_\Sigma\) sum the 105 junction kernels. Let \(K_\Sigma\) sum all
630 labelings of the cubic path-tree with five internal vertices,
quotiented by its two endpoint cherry swaps and path reversal.
Each such tree kernel supports 96 words, or sixteen orbits.
The orbitwise counts are

| Boundary multiplicity | \(J_\Sigma\) | \(K_\Sigma\) | \(\sum_iw_i\) |
|---|---:|---:|---:|
| \((5,1,1)\) | 15 | 120 | 2 |
| \((3,3,1)\) | 9 | 108 | -1 |

For a fixed tree, a majority colour appearing on five ports must occur
at all five internal vertices. Each endpoint cherry then contains one
majority port and one singleton: there are \(2\cdot2=4\) such orbits,
and alternating the two other colours along the internal path verifies
all four. The remaining twelve orbits have type \((3,3,1)\).
Port-transitivity gives \(630\cdot4/21=120\) and
\(630\cdot12/70=108\); the checker also directly enumerates the table.
Thus
\[
 11K_\Sigma-112J_\Sigma=-180\sum_iw_i.
\]
For labelled-colour closure counts this gives
\[
 11C_{\rm path630}(P)-112C_{\rm junction105}(P)
       =180\sum_iN_i\ \geq 0.                               \tag{4}
\]
Equality holds exactly when all seven \(N_i\) vanish, equivalently
when \(f_P\) belongs to the rank-84 row space. Equality alone does not
identify the smaller target family \(mJ\).

There is also a positive physical-ensemble hostile: eleven copies
of every one of the 630 tree tensors and 112 copies of every one of
the 105 junction tensors have identical 105 junction measurements.
Their difference is the nonzero vector \(-180\sum_iw_i\).
These are integer ensembles of physical poles, not a claim that either
sum is realized by one connected pole or a selector gadget.

## 4. Golden weighting preserves zero but loses real physical information

Let \(\phi=(1+\sqrt5)/2\), and define
\[
              G(N)=\sum_{i=0}^{6}N_i\phi^i.
\]
Because all \(N_i\geq0\), \(G(N)=0\) iff all \(N_i=0\).
This is an exact useful zero test, but any strictly positive weights
have that property. Golden weights do not reconstruct \(N\).

Here are two **connected single-pole** witnesses, not merely formal
kernel vectors. The first has internal path \(0-1-2-3-4\), with ports
incident in label order to
\[
                 (0,0,2,1,3,4,4).
\]
Its defect vector is \(N=24e_2\). For the second use the Petersen graph
with edge set
~~~text
01 04 05 12 16 23 27 34 38 49 57 58 68 69 79
~~~
delete vertices \(0,1,3\), and order its port incident vertices as
\[
                 (5,6,4,2,2,4,8).
\]
Its defect vector is \(N'=24(e_0+e_1)\). Both edge-colouring enumeration
and the independent matching formula (3) verify these values.
They have different extension tensors and defect vectors, but
\[
             G(N)=24\phi^2=24(1+\phi)=G(N').
\]
Their total boundary colouring counts are 96 and 120, respectively.
No equality of their 105 junction readouts is asserted.

The missing information is exactly the relation
\(\phi^{i+2}=\phi^{i+1}+\phi^i\). A Zeckendorf normal form can select
one representation only after constraints on the coefficients are
specified; these actual physical defect vectors violate the conditions
needed to infer the original seven coordinates from one golden number.

## 5. A precise golden diagram map, and two decisive hostiles

The primary source opened for this comparison is Fendley--Krushkal,
[Tutte chromatic identities from the Temperley--Lieb algebra](https://arxiv.org/pdf/0711.0016),
arXiv:0711.0016v3. Equation (1.1) is the golden identity for planar
triangulations. Section 2 relates the chromatic parameter \(Q\) to the
loop parameter \(d\) by \(Q=d^2\); Section 5.3 discusses the trace
radical at the special loop parameters. The present finite calculation
is independent of a general diagram-algebra rank theorem.

Take all 105 perfect matchings of eight cyclically ordered labels.
For two matchings let \(\ell(M,N)\) count components of their union,
including a doubled common edge. Define
\[
                     G_d(M,N)=d^{\ell(M,N)}.
\]
The 14 noncrossing matchings form a specified subset of this *same*
carrier. The embedding preserves the entire loop-count polynomial.
It does not make crossing probes redundant.

The inherited result and the new exact checker together give:

| Probe carrier / loop parameter | Rank |
|---|---:|
| All 105 matchings, \(d=3\) | 91, inherited from the 105/17 packet |
| The 14 noncrossing matchings, \(d=3\) | 14 |
| The 14 noncrossing matchings, \(d=\phi\) | 13 |
| All 105 matchings, \(d=\phi\) | 105 |

For the golden 14-by-14 Gram matrix, the checker exhibits a nonzero
exact null relation over \(\mathbb Q(\phi)\), verifies it against every
planar probe, and prints its coefficients in its fixed matching order.
Pad this vector by zero on the other 91 matchings. The crossing probe
\[
                 (02)(13)(46)(57)
\]
evaluates it to \(-3+\phi\ne0\). Even within the planar carrier, changing
the loop parameter to three makes the probe
\[
                 (01)(23)(45)(67)
\]
evaluate it to \(-42-6\phi\ne0\). These are explicit losses of the
proposed transport, not just different matrix dimensions.

The field calculation uses the representation \(a+b\phi\) and its exact
28-by-28 rational multiplication matrix. Rank 13 has both an exact null
vector and a lower bound modulo 11, under \(\phi\mapsto4\).
The full 105-by-105 golden matrix is nonsingular under that same
specialization, proving its field nonsingularity. Thus a planar golden
trace relation cannot be imported unchanged into the ordinary ternary
or crossing-permitted boundary problem.

There is nevertheless an exact Fibonacci grammar on these diagrams:
noncrossing matchings of \(2m\) ordered ports whose nesting height is
at most three correspond to Dyck walks on heights \(0,1,2,3\).
Their two-step transition matrix on even heights is
\[
                M=\begin{pmatrix}1&1\\1&2\end{pmatrix}.
\]
Hence their number is \((M^m)_{00}=F_{2m-1}\) for \(m\geq1\).
For 2 through 14 ports it is
\[
                     1,2,5,13,34,89,233.
\]
This gives 233 a precise address: fourteen ports in this restricted
nesting grammar. At eight ports its thirteen diagrams omit only the
fully nested rainbow from the fourteen planar diagrams. This packet
does not infer an all-order Gram-rank formula from that observation,
nor identify this grammar with tournament classes or Collatz orbits.

## 6. Reproduction and remaining strongest target

From the repository root:
~~~text
python 05-knowledge/results/deficiency_three_20261020/boundary_check.py
python -O 05-knowledge/results/deficiency_three_20261020/boundary_check.py
~~~
The receipt is boundary_check.out. Checks use explicit failures, not
Python assertions; the two outputs must agree byte for byte. Dependencies
are NetworkX and SymPy. No numerical optimizer is needed.

The universe includes all 2187 seven-port words; all 105 target
relabelings; all 630 path-tree labels; all 105 eight-port matchings and
their 14 planar members; and all 19 connected simple internal graphs
of orders five or seven in the NetworkX graph atlas that become cubic
seven-poles after adding their missing incident ports. Independent
internal-edge colouring and matching/path-parity computations agree on
every one of those poles. The connected 17-vertex realization is also
checked directly: its tensor is \(4J\), all fourteen probes vanish,
and its target self-closure is \(216=54\cdot4\).

The positive direction is now specific: use the fifteen target tests
while retaining a separate, checkable rooted-minor or disjoint-linkage
certificate. General full-tensor recovery is unnecessary for this
replacement target. The remaining graph-side obligation is not removed
by the new count certificate, a golden specialization, or an aggregate
zero test.
