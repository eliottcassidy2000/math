# Joint path information, the first Selmer obstruction, and the F8 symmetry boundary

**Status: PROVED IN NOTES, corrected B1–B3, after independent audit.**
The saved computation is **FINITE-EXACT**. Inherited proofs and cited arithmetic
inputs are identified separately. No general BSD, Hadwiger–Nelson, Hamiltonian
spectrum, or new realization theorem is claimed.

## Inheritance and the actual interfaces

The anchor is the current joint two-edge path reconstruction. The niche is
the first Selmer statistic of a *marked residue tournament*. The wildcard is
the exact restriction of the Paley-blocker F8 representation to its C3
complement. The live board is joint paths / slot polynomials / first descent /
Frobenius submodules / palette boundary relations / retained realization.

Three mechanisms, rather than a shared integer, connect the inherited work:

1. [Finite Selmer growth](../../bsd_selmer_growth_20261013a.md), §§1,3–4,
   proves a finite-growth test for the full prime-power Selmer corank and
   identifies a particular first-descent matrix. Its little-used sidecar is
   [the exact strong-component matrix decomposition](../../bsd_component_recursion_20261013c.md),
   §1. Current [THM-5125, residue tournament model](../../../../01-canon/theorems/THM-5125-residue-tournament-model-heegner-domination-and-cycle-gas-identities.md),
   §4(c), also identifies the obstruction with a marked cycle-packing parity
   balance. These are recovered results, including their converse invariance.
2. [HN palette/deck recovery](../../session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md),
   §§1,3–5, proves an exact bit-lift criterion for a **fixed** three-channel
   assignment; it retains all channels and the existential choice of that
   assignment. Its thirty Fano kernels have identical pair marginals but
   different cubic correlations. [The subsequent information packet](../../hn_metagraph_information_20261026/information/RESULT.md),
   §§3–4, proves that full HP-flip-response quadratic Walsh coefficients,
   with common labels, recover a tournament up to converse. It also proves
   that separately quotienting two faithful aligned records loses information.
3. [The corrected Paley-7 blocker bridge](../../session54_reframed_numerology_20261025/sixtythree-bridge/README.md)
   supplies a literal G-set map from the 63 minimum six-arc blockers to the
   nonzero binary augmentation module, then its F8 Fourier coordinates.
   Its five-type decoder is essential: the map does not preserve native
   deleted-arc XOR or arbitrary deletion-path counts.

Canonical hostiles remain operative: E627's pure Selmer dimension two does
not imply rank two; the geometric all-negative double cover can be bipartite
without having a section; and a native blocker family need not be a binary
code. MISTAKE-752 records the known Legendre-tournament priority, the failed
implication `nu>0 => 3 divides H`, and an unsupported geometric reading of an
order-three action. None is revived below.

The current [five-manuscript source-status reference](../../../reference/CORE-PAPERS-OPENAI-FIVE-MANUSCRIPTS-2026-10-13.md)
marks the supplied low-corank BSD manuscripts **PREPRINT CLAIM**, not proved
canon. The older BSD packets' full-BSD applications were conditional on an
accepted session input. This note imports neither that input nor its full-BSD
conclusion. The unchanged-source Tunnell terminals in
[the prior exact arithmetic packet](../../bsd_tunnell_terminals_20261013c.md)
are a different certificate from first descent; substituting a new curve is
not a certificate for the old one.

## B1. A concrete arithmetic kernel consumer of reconstruction up to converse

Fix a labelled tournament T on n vertices. Let A be its zero-diagonal binary
adjacency matrix, with A_ij=1 for i→j, and let

```
M(T) = A + diag(outdegree(T) mod 2),
Phi(X) = X²+X+I,                     nu(T) = dim_F2 ker Phi(M(T)).
```

The rows of M sum to zero. This is the row-zero matrix, not the full 2n by 2n
Monsky matrix. Passing to converse gives

```
M(T^op) = M(T)^t + ((n−1) mod 2) I.
```

Since `Phi(X+I)=Phi(X)` in characteristic two,
`Phi(M(T^op))=Phi(M(T))^t`. In fact the full matrices are equal, a stronger
identity recovered from the incoming
[Session321 BSD parity note](../../session321_selective_returns_and_fano_wiring_20261113/BSD_PARITY.md).
Here is the direct proof, retaining its source credit. With J all-one,

```
M^t=M+J+I,        MJ=0,        JM=(n−1)J,        J²=nJ.
```

Expanding Phi(M+J+I) gives Phi(M): the J coefficient is
`(n−1)+n+1=0` in F2. Hence Phi(M) is symmetric and
**Phi(M(T^op))=Phi(M(T)) on the same vertex labels**.
Any aligned observation reconstructing T up to converse therefore recovers
the full Phi matrix, its embedded kernel in F2^n, and nu. Under relabelling,
the matrix and kernel transform by the common vertex permutation, so nu also
descends to unlabelled converse classes. M itself still needs orientation.
For the full Monsky embedding the exact remaining dependence has a parity:
since Phi(M) is symmetric and Phi(M)1=1, every kernel vector x has coordinate
sum zero. From `M(T^op)=M(T)+J+(n mod2)I` we obtain
`M(T^op)x=M(T)x+(n mod2)x`. Thus the embedding `x↦(x,(M+I)x)` changes by
`(0,(n mod2)x)`: it agrees at even n, and may change at odd n. The underlying
named F2 kernel never changes.
The initial weaker proof falsely assigned that gauge loss to the Phi kernel;
the explicit repair lineage is in [CORRECTION.md](CORRECTION.md).

The actual joint observer is explicit. For distinct u,v,w retain the bit
`beta(v;u,w)` saying that some spanning cover in two ordered optional slots
uses both pairs {u,v},{v,w}. The inherited
[Session320 wedge theorem](../../session320_second_moments_and_native_contexts_20261113/WEDGE_RECONSTRUCTION.md)
identifies it with

```
beta(v;u,w) = A_uv A_vw + A_wv A_vu in {0,1}.
```

Necessity is local; sufficiency puts that directed three-vertex path in one
slot and a Hamilton path of the remaining induced tournament in the other.
The bit determines the relative orientations of the two edges. Integrating
along a spanning tree of the connected line graph L(K_n) recovers T up to
converse for n≥3. The same bit is the adjacent-pair coefficient of the
[labelled quadratic path-cover jet](../../session320_second_moments_and_native_contexts_20261113/F1_WEDGE_JET.md).
Applying Phi after this decoder is the promised literal path-data→matrix/kernel
map, with nu as a corollary. The current reconstruction lane audits and extends
the observer's image and error correction; it is not needed as an unproved
premise for this inherited exact decoder.

The consumer also applies to the older full-response decoder. These are
different observation models: no local two-edge bit is identified with a
degree-two Walsh coefficient averaged over the full flip cube.

For arithmetic, retain distinct prime labels p_i≡3 (mod 8), set
`N=product p_i`, and orient i→j when p_j is a quadratic **nonresidue** mod p_i.
Reciprocity makes this a tournament. For `E_N: y²=x³−N²x`, the cited Monsky
formula gives

```
full Monsky matrix = [ M+I  I ],
                     [  I   M ]
ker Phi(M) → ker Monsky,   x ↦ (x,(M+I)x).
```

Consequently nu is the pure 2-Selmer dimension, namely
`dim_F2 Sel_2(E_N/Q) − dim_F2 E_N(Q)[2]`, with the latter torsion dimension two.
The cited input is [Watkins–Donnelly–Elkies–Fisher–Granville–Rogers, §4.1](https://dms.umontreal.ca/~andrew/PDF/RankQuadTwist.pdf),
primary PDF checked 2026-11-13, printed p.8. The displayed specialization and
kernel map are the elementary inherited calculation. A binary matrix without
the prime labels, residue restrictions and this local-condition theorem is
not an authenticated Selmer computation for an arbitrary curve.

Finite Kummer exactness says `nu = rank E_N(Q) + dim Sha(E_N/Q)[2]`.
Thus nu=0 proves rank zero and no Sha[2]; nu=2 does not say which term supplies
the two dimensions. No full Selmer corank, analytic rank, entire Sha order,
regulator, or BSD leading term is reconstructed here.

## B2. Full slot Q and scores still lose this consumer, at every order n≥5

For j≥1, let F_j(T) count ordered covers of V(T) by j nonempty directed paths.
Define the full slot polynomial

```
Q_T(t) = sum_j F_j(T) binom(t,j).
```

For each nonnegative integer t, it counts allocations into t **numbered**
slots, with an ordered directed path in each nonempty slot. This description
retains empty slots and supplies the following product without interpolation
ambiguity:

**Ordinal-join identity.** For disjoint tournaments U and V, with every cross
arc directed U→V,

```
Q_(U→V)(t) = Q_U(t) Q_V(t).
```

Proof: in each numbered slot, a path of U→V has all its U vertices followed
by all its V vertices. It is uniquely the concatenation of a possibly empty
U-path and a possibly empty V-path. Splitting and concatenation are inverse
for each integer t≥0; hence the polynomial identity follows. This is an
operation identity for literal slots, not merely a factorization of H=Q(1).

The inherited matrix counterpart is also exact. If b=|V| and J is all-one,

```
M(U→V) = [ M(U)+(b mod2)I   J    ],
          [        0        M(V) ],
Phi(M(U→V)) = Phi(M(U)) direct-sum Phi(M(V)).
```

Indeed M(U)J=0, JM(V)=(b−1)J, and the off-diagonal term in Phi is
`bJ+(b−1)J+J=0`; the diagonal shift is killed by Phi(X+I)=Phi(X).
Thus nu(U→V)=nu(U)+nu(V). This matrix identity was already proved in the
component-recursion packet; its use with equal full Q is the present bridge.

Use vertices 0,…,4 and lexicographic unordered-pair bits, with bit one meaning
the increasing arc. Let T8 and T10 be masks 8 and 10. Equivalently T8 has
increasing arc 0→4 only, and T10 has increasing arcs 0→2 and 0→4 only;
all other arcs run from the larger label to the smaller.

| Invariant | T8 | T10 |
|---|---:|---:|
| H | 9 | 9 |
| Sorted scores | (1,1,2,3,3) | (1,1,2,3,3) |
| (F1,F2,F3,F4,F5) | (9,66,186,240,120) | (9,66,186,240,120) |
| Full Q | t⁵+6t³+2t | t⁵+6t³+2t |
| nu | 2 | 0 |

These are the same two classes as the minimal full-Q/native-incidence hostile
in the prior [incidence closure packet](../../incidence_closure_20261113/closure/PROOF.md).
Both are self-converse with automorphism order one. Those extra facts are
inherited metadata, not needed to prove B2. The present independent checker
rebuilds the path covers and both binary matrices rather than importing those
computed rows.

An explicit rank certificate is

```
Phi(M(T8)) = 11111     Phi(M(T10)) = 10110
             11001                   01101
             10101                   11100
             10011                   10011
             11111                   01011
```

The row strings list coordinates 0,…,4 left to right. The first matrix's
kernel is `{00000,11110,10001,01111}` in that coordinate order, of dimension
two. The second is invertible. The JSON also stores vector bitmasks
`[0,15,17,30]` and `[0]`, avoiding a string-endianness convention.

**All-order failure.** For n≥5, form `U_n=T8→TT_(n−5)` and
`V_n=T10→TT_(n−5)`, where TT_0 is empty. A transitive tournament on r vertices
has Q=t^r: each assignment of vertices to numbered slots has a unique valid
order within each slot. Its row-zero M is triangular with diagonal entries
zero or one, so Phi(M) is triangular with all diagonal entries one and nu=0.
Therefore

```
Q_(U_n)=Q_(V_n)=t^(n−5)(t⁵+6t³+2t),
nu(U_n)=2,                         nu(V_n)=0.
```

Their score multisets agree: add n−5 to the five old scores and append
`0,1,…,n−6`. Hence **even (n, full Q, scores) fails to determine nu for every
n≥5**. This is a proved-construction candidate, not an extrapolated census.
The starting order is minimal: the independent exhaustive controls contain
all 1+2+8+64 labelled tournaments of orders1,…,4, where even H determines nu.
At order5 the only H values with differing nu are9 and15; at15 the pair of
masks40/76 gives nu4/0. Minimality is FINITE-EXACT, separate from the all-n
construction.

The following actual arithmetic realizations make the map nonempty:

| Code | Ordered prime labels | N | Pure 2-Selmer dimension |
|---|---|---:|---:|
|8|3,19,43,139,251|85,512,939|2|
|10|3,19,83,499,1163|2,745,574,347|0|

All primes are distinct and3 mod8; every specified Legendre symbol is checked
both by Euler's criterion and by literal square enumeration. The full10 by10
Monsky matrices are checked independently against Phi(M). These are two
different curves with equal path data, not two descriptions of one curve.
The inherited CRT/Dirichlet realization theorem realizes every joined
tournament too; no prime-size bound or frequency claim is needed.

## B3. The F8→F4 comparison exists only after restricting symmetry

Fix `F8=F2[alpha]/(alpha³+alpha+1)`. Let W be the even subsets of F7 under
symmetric difference. The inherited Fourier isomorphism is

```
W → F8 ⊕ F8,
S ↦ (sum_(s in S) alpha^s, sum_(s in S) alpha^(−s)).
```

For `G=C7:C3={s↦2^k s+b}`, the exact action is

```
(u,v) ↦ (alpha^b u^(2^k), alpha^(−b) v^(2^k)).
```

Here t (b=1,k=0) is the C7 translation generator; sigma (b=0,k=1) is the C3
**Frobenius complement**, not multiplication by a scalar of order three in
F8*. The latter group has order seven.

On either coordinate `sigma²+sigma+I` equals the trace F8→F2. Consequently

```
ker Phi(sigma | W) = Z ⊕ Z,
Z=ker Tr_(F8/F2)={0,alpha,alpha²,alpha⁴}.
```

Z is a two-dimensional F2-space and a one-dimensional F4-module when a root
omega of X²+X+1 acts by sigma. For the chosen alpha, an explicit module map is
`1↦alpha`, `omega↦alpha²`, `omega²↦alpha⁴`. Thus the kernel is F4² as a module
for this C3. This is exactly the same irreducible-polynomial mechanism that
makes ker Phi(M(T)) an F4-space. It identifies representations after choosing
an operator and a frame; it does not identify them as Selmer groups.

**Sharp obstruction to extending the comparison.** Z⊕Z is not G-stable:
`(alpha²,0)` has trace zero but its t-image `(alpha³,0)` has trace one.
More generally W has no nonzero G-equivariant **linear quotient** of F2
dimension2 or4. To prove this, restrict to C7. Its two summands have distinct
irreducible minimal polynomials

```
X³+X+1,                    X³+X²+1,
```

neither of which has an F2-root. The Chinese-remainder projections in F2[t]
split every invariant subspace into its intersections with these summands;
each intersection is zero or its entire simple three-dimensional summand.
The only C7-invariant subspaces are therefore0, the first axis, the second
axis, and W. Their quotient dimensions are6,3,3,0. Requiring sigma invariance
cannot add any. In particular no G-equivariant linear projection can turn
the recovered four-dimensional C3 sector into a quotient of W.

This is a category-specific obstruction: it does not forbid set maps,
nonlinear encodings, or restriction to the C3 complement. Nor does it mean
F4 is a subfield of F8; the displayed trace-zero set even lacks the field
identity1. The roots-of-unity/module action is not an automorphism of the
elliptic curve or its torsion.

### The genuine HN/Fano connection, and the correlation-order guard

The three nonzero elements of Z are the exponents1,2,4 in the alpha frame.
Their seven multiplicative translates are exactly the seven additive Fano
lines in F8\{0}; each pair lies in one line. Thus the trace-zero sector has
a concrete incidence interpretation, with the chosen frame retained.

For the inherited thirty labelled simplex codes, write
`s_i=a dot v_i`, where a ranges uniformly over F2³ and v_i enumerate its
seven nonzero vectors. Then the means of `(-1)^s_i` and all distinct-pair
correlations are zero, while

```
E (-1)^(s_i+s_j+s_k) = 1 iff v_i+v_j+v_k=0, and 0 otherwise.
```

This recovered character orthogonality says pair marginals cannot distinguish
the thirty incidences, while cubic correlations do. The checker replays all
thirty codes. It does not conflict with second-order HP reconstruction:
those observations are responses of a **different function on a different
source**, retaining a different system of shared labels. The operation and
the observer must accompany any claim about the order of information.

For actual HN colourings, even a faithful abstract incidence is insufficient.
With a fixed ternary assignment q, the old exact lift uses the edges H_q
whose q-labels agree and solves `delta b=1` there. Six-colourability quantifies
over q. Pairwise projections of the union over q can lose the obstruction;
neither the F8 module nor an HP decoder supplies the unit-distance realization
or removes this existential quantifier. These are retained boundaries of the
old proved mechanisms, not an uncompleted promise of a new plane bound.

## Reproduction and scope

From the repository root:

```
python3 -B 05-knowledge/results/integer32_twoedge_20261113/bridges/check_bridges.py
python3 -B -O 05-knowledge/results/integer32_twoedge_20261113/bridges/check_bridges.py
```

[normal.json](normal.json) and [optimized.json](optimized.json) are identical:
**19,447 always-active gates**. The checker imports no prior research code.
It compares permutation HP counts with subset DP, row reduction with complete
binary-kernel enumeration, and word/cut path covers with ordered subset-cover
DP. It also enumerates literal two-slot permutation/cut covers at every
labelled tournament through5, compares their wedge supports with the local
formula, decodes along a fixed line-graph tree and verifies the resulting
full Phi matrix, embedded kernel and nu. Full matrix and kernel invariance
under converse are separately checked throughout the same universe.
It checks the literal join matrices and Q product on all121 ordered pairs
of labelled tournaments of orders1..3; the all-n hostile is additionally
probed at orders5..9. These finite probes supplement the displayed proof.

Its other universes are all64 augmentation states under all21 affine group
elements, every generated cyclic C7 submodule, and all30 labelled Fano codes.
Source SHA256 pins and checker SHA256 are in both receipts. The normal and
optimized agreement verifies always-active checks, not independent
implementations; the explicitly different algorithms above supply that part
of the evidence. There is no large census or imported n7 matrix dependency.

The concrete next consumer is the first-descent matrix and kernel computable
from a faithful joint path observer; even their dimension is not determined
by full slot Q. The next realization
obligations remain literal higher Selmer data on the unchanged curve and
literal unit-distance/full-colouring relations in HN. The F8 restriction
describes exactly which symmetry must be relinquished, rather than supplying
either realization.
