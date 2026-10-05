# Four tournament states, recursive storage, and typed doubling

**Status:** PROVED for the stated graph and matrix identities; FINITE-EXACT for
the enumerations below. The graph codec stores words. It does not turn an
arbitrary source into a completed Collatz certificate.

## 1. Inheritance and the four-state boundary

The closest earlier mechanisms are the fixed-path half-tiling in
[HYP-3244](../hypotheses/HYP-3244-tiling-half-tiling-interlocking-recursions.md),
the endpoint padding in
[THM-291](../../01-canon/theorems/THM-291-mode-b-multilinear-recursion.md),
the order join in
[THM-1862](../../01-canon/theorems/THM-1862-order-join-reduction-principle.md),
and the exact skew block matrix in
[THM-447](../../01-canon/theorems/THM-447-skew-sylvester-doubling.md).
Their operations act on different coordinates. The historical owner request
behind skew doubling was \(n*2\), which should not silently become \(n^2\).

The relevant hostile is already visible on four vertices: two presentations
of the same isomorphism class can respond differently to the same named tile
flip. The corrected near miss is the earlier half-band Fourier claim:
MISTAKE-034 in [MISTAKES](../../01-canon/MISTAKES.md) records a degree-four
coefficient for the six-bit \(n=5\) cube. Neither a universal half-band cutoff
nor a tile-complement/converse identification is used here. The least-used
useful sidecar is the ordered strongly connected component decomposition.

Fix the Hamiltonian path \(0\to1\to2\to3\). Start with every arc forward and
let \(a,b,c\) flip the pairs \(02,13,03\), respectively. The eight fixed-path
presentations have the following fibers:

| Class | Sorted outdegrees | Fixed-path masks | Hamiltonian paths \(H\) | Automorphisms |
|---|---|---|---:|---:|
| \(T\), transitive | \(0,1,2,3\) | \(0\) | 1 | 1 |
| \(+\), cycle over a sink | \(0,2,2,2\) | \(a\) | 3 | 3 |
| \(-\), source over a cycle | \(1,1,1,3\) | \(b\) | 3 | 3 |
| \(S\), strongly connected | \(1,1,2,2\) | \(ab,c,ac,bc,abc\) | 5 | 1 |

The script independently classifies all 64 labelled four-vertex tournaments
under all 24 permutations. In this fixed-path chart the fiber size is
\(H/|\operatorname{Aut}|\): an ordering that exhibits the prescribed path is
a Hamiltonian path, and automorphisms act freely on these orderings.

Both \(ab\) and \(c\) represent \(S\), but flipping \(a\) gives \(b\), of class
\(-\), and \(ac\), still of class \(S\). Thus the four classes are not a
quotient on which all three named flips act.

In contrast, the section \(c=0\), with masks \(0,a,b,ab\), is closed under
the two flips \(a,b\) and contains one representative of each class. It
therefore carries a valid four-state XOR action. This is a statement about
that section and its two allowed moves. Moving through \(c\), or replacing
an arbitrary presentation by the section representative of its class,
requires retaining the presentation if future named flips matter.

**General predictive-memory lemma.** For a fixed Hamiltonian path on \(n\)
vertices, let \(C=\mathbb F_2^{\binom{n-1}{2}}\) be its free-arc cube. A state
map which preserves tournament isomorphism class and admits deterministic
updates for every named free-arc flip is injective on \(C\).

Proof: the zero presentation is the unique transitive presentation, since a
transitive tournament has a unique Hamiltonian ordering and the path fixes
that ordering. If two distinct masks \(u,v\) had the same state, applying the
sequence of flips \(u\) would give equal states at \(0\) and \(u+v\). Only
the former is transitive, a contradiction. This does not restrict operations
already well-defined on unlabelled isomorphism classes.

Full bit complement also fails to model converse: it sends \(0\) to \(abc\),
hence \(T\) to \(S\). True converse followed by reversal of the path labels
exchanges \(a,b\), retains \(c\), and exchanges the two diamond classes.

## 2. A faithful recursive four-symbol carrier

Choose the four representatives \(T,+,-,S\) from \(c=0\). For a word
\(w=c_1\cdots c_k\), define
\[
 \mathcal T(w)=T_{c_1}\to T_{c_2}\to\cdots\to T_{c_k},
\]
where every cross arc points from the earlier block to the later block.
The empty word has the empty graph. This is a four-symbol coding convention;
the symbols have no intrinsic Collatz meaning.

**Lossless decoding theorem.** The unlabelled tournament \(\mathcal T(w)\)
determines \(w\), without the original labels or the construction partition.

Every tournament's strongly connected components have a unique total order.
An SCC cannot cross an order-join cut. Traverse this order and cut after each
four vertices; the original cuts guarantee that no SCC crosses a required
cut. Classify each resulting four-vertex block. This recovers every symbol.
Conversely the recovered word reconstructs an isomorphic tournament. For a
general input graph the decoder must reject if an SCC crosses a required
four-vertex cut; arbitrary tournaments are not asserted to belong to this
language.

Concatenation is exact:
\[
 \mathcal T(uv)=\mathcal T(u)\to\mathcal T(v).
\]
A recursively repeated word can be stored as a shared expression, and its
graph can be expanded when required. This saves storage for structured
repetition, not for arbitrary unrelated words; the expanded graph still has
four vertices per letter.

An order-join Hamiltonian path must finish the left block before entering
the right one, so
\[
 H(T\to S)=H(T)H(S),\qquad
 H(\mathcal T(w))=3^{\#\{+,-\}}5^{\#S}.
\]
This invariant discards both order and the diamond distinction. For example,
\((+,-)\) and \((-,+)\) both have \(H=9\), while their intrinsic SCC size
words are \((3,1,1,3)\) and \((1,3,3,1)\).

For an arithmetic word one may supply an exact serialization into four
symbols, including any delimiters needed for unbounded exponents. The graph
then stores that serialization faithfully. If each arithmetic letter acts
by a matrix, chronological concatenation satisfies \(M_{uv}=M_vM_u\);
repeating the same word still gives \(M_{ww}=M_w^2\). The graph decoder
retains the order that a scalar invariant such as \(H\) or a trace discards.
The parser, actual source, divisibility guards, and first-hit condition are
separate typed data. This is storage for a certificate, not a certificate
obtained from a graph statistic.

The earlier [route-tournament codec](collatz_atom_route_memory_20261003.md)
uses odd regular SCC blocks and their ordered sizes to store inverse-route
digits. Its \(+2\) source-block rule and transitive cloning rule have
arithmetic meanings only with that codec's guards. Internal \(+2\) can fail
the ternary guard, for instance \((0,1)\to(0,2)\). The present four-symbol
codec is another explicit carrier, not a replacement of those guards.

## 3. Three distinct meanings of growth

Write \(n=|T|\), and let \(K_1\) be a singleton.

| Operation | Definition | Vertex count | Exact retained/transformed statistic |
|---|---|---:|---|
| Endpoint padding | \(K_1\to T\to K_1\) | \(n+2\) | \(H\) unchanged |
| Self order join | \(T\to T\) | \(2n\) | \(H\mapsto H^2\) |
| Self substitution | \(T[T]\) | \(n^2\) | ordered pair addresses; not generally \(H^2\) |
| Transitive cloning | \(T[TT_2]\) | \(2n\) | marked or intrinsically recoverable twin fibers required for inverse use |
| Skew doubling | \(\mathcal D(M)\), below | \(2n\) | exact squared-matrix identity |

In substitution \(T[S]\), a vertex is a pair \((i,j)\): between distinct
outer indices use \(T\), and within one outer index use \(S\). For the
directed three-cycle, \(H(C_3\to C_3)=9\), whereas
\(H(C_3[C_3])=3159\). Even the two doubling constructions disagree at the
smallest useful control: \(TT_2[TT_2]\) is transitive with \(H=1\), while
skew doubling \(TT_2\) is a source over a three-cycle with \(H=3\).

**Local count correction to the inherited padding statement.** The true
THM-291 restriction is \(H_{n+2}(t_{\rm inner},0)=H_n(t_{\rm inner})\).
The added free boundary has
\[
 \binom{n+1}{2}-\binom{n-1}{2}=2n-1=(n-1)+(n-1)+1
\]
bits, not \(2n-3\). With old vertices \(2,\ldots,n+1\), the arms are
\(x\to1\) for \(x=3,\ldots,n+1\) and \(n+2\to y\) for
\(y=2,\ldots,n\), plus the apex \(n+2\to1\). Setting all these free bits
to the compatible endpoint orientation gives the padding identity. The
original file's arm endpoints/count have been repaired with the correction
lineage retained in MISTAKES; the path-count proof itself survives.

## 4. The exact two-matrix, four-product cancellation

Let \(M\) be a tournament's skew dominance matrix, and put
\[
 A=\begin{pmatrix}1&1\\1&-1\end{pmatrix},\qquad
 B=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]
Then \(A^2=2I\), \(B^2=-I\), and \(AB+BA=0\). The inherited skew doubling is
\[
 \mathcal D(M)=
 \begin{pmatrix}M&M+I\\M-I&-M\end{pmatrix}
 =A\otimes M+B\otimes I .
\]
It remains skew with off-diagonal entries \(\pm1\), hence is a tournament
dominance matrix. Its four product terms give
\[
 \mathcal D(M)^2
 =A^2\otimes M^2+(AB+BA)\otimes M+B^2\otimes I
 =I_2\otimes(2M^2-I).
\]
These are four algebraic product terms, not the four isomorphism classes.
No intrinsic bijection between those two four-element descriptions is
asserted.

Iterate \(M_0=M,\ M_{d+1}=\mathcal D(M_d)\). Induction gives the exact law
\[
 \boxed{\quad
 M_d^2=I_{2^d}\otimes
       \bigl(2^dM^2-(2^d-1)I\bigr).
 \quad}
\]
Thus if a skew eigenvalue is \(i\sqrt z\), its squared magnitude transforms
as \(z\mapsto2z+1\); equivalently \(z+1\) doubles. This is not the map
\(\lambda\mapsto\lambda^2\) for a Collatz word's positive affine slope.

The [quadratic word-matrix result](quadratic_escape_rank_atlas_20261004.md)
does give \(\lambda\mapsto\lambda^2\) and
\(J=\lambda+\lambda^{-1}\mapsto J^2-2\) under repetition of the same word.
The exact bridge here is the ordered-word carrier in section 2: it can
retain the word being repeated while its arithmetic evaluator follows that
quadratic identity. The skew construction supplies a different, explicitly
typed cancellation law. Neither gives a common global orbit dynamics merely
from their similar-looking formulas.

## 5. Reproduction, universe, and connection contract

Run from the repository root:

    python 04-computation/experiments/tournament_recursive_four_state_20261004.py
    python -O 04-computation/experiments/tournament_recursive_four_state_20261004.py

The [script](../../04-computation/experiments/tournament_recursive_four_state_20261004.py)
uses explicit exceptions, so its checks remain active under the optimized
interpreter. The [saved output](tournament_recursive_four_state_20261004.out)
records:

- all 64 labelled four-vertex tournaments and all 24 relabelings;
- the eight fixed-path masks and 2,045 pair distinctions for sizes 2–5;
- all 341 four-symbol words of length 0–4, each under three relabelings,
  plus independent path-count controls through length 2;
- the typed growth operations on six seeds, the two small doubling
  hostiles, and the corrected padding boundary count for sizes 1–100;
- all 75 labelled tournaments of sizes 1–4 at skew depth 0–3, giving
  300 exact matrix identities;
- five malformed or out-of-language controls.

For the storage connection the source is an ordered finite word, target is
an unlabelled order-join tournament, map is the four-block encoder, and the
preserved predicate is exact word equality under decoding. Vertex labels
and unnecessary construction partitions are discarded. General source
parsing and arithmetic guards remain sidecars. For the named-flip quotient,
the minimal failed implication is that equal present class predicts equal
future class; retaining all free bits repairs it. These are exact finite
and general statements about representations, not extra Collatz coverage.
