# A retained-endpoint consumer and a positive cycle-component compiler

**Session264, 11 November2026. PROVED:** retained-endpoint/position consumer,
triangle-wrap operator on the full coloured response, odd-component merging,
and a finite-interface compiler. **FINITE-EXACT:** the explicit eleven-vertex
consumer and the stated bounded searches. **OPEN:** a physical scalar selector
of one residue, a primitive-sixteen block, and arbitrary-target coverage.
Inherited inputs are assumed at their stated types. No external novelty claim,
canon promotion, or inherited-result audit is made.

## 1. A physical endpoint consumes the length response

For a tournament A with an existing selected vertex y, write

    P_A,y^+(z) = sum_(paths R starting at y) H(A−R) z^|R|,
    P_A,y^−(z) = sum_(paths R ending at y) H(A−R) z^|R|.

Paths are nonempty; H(empty)=1. Construct W(A) by adding u,v with

    u → A → v → u.

**Theorem 1 (retained endpoint and visible position).**

- The Hamiltonian paths of W(A) starting at the SAME vertex y are in bijection
  with a y-start path R of A and a Hamiltonian path S of A−R. The map is
  `(R,v,u,S)`, where S may be empty. The number of A vertices before v is |R|.
- The Hamiltonian paths of W(A) ending at the SAME y are in bijection with a
  y-end path R of A and a Hamiltonian path S of A−R. The map is `(S,v,u,R)`.
  The number of vertices after u is |R|.

Consequently the generating polynomial for these literal Hamiltonian paths,
recording the displayed prefix or suffix length, is respectively P^+ or P^−.
In particular the starting count is U_A(y) and the ending count is V_A(y).

**Proof.** In a Hamiltonian path starting at y∈A, the vertex u cannot occur
before v: its only incoming neighbour outside A is v and no A vertex beats u.
After v occurs its only available outgoing neighbour is u, so v,u are adjacent.
All vertices before v form a nonempty y-start path R in A; all vertices after
u form a path S in the complementary set, possibly empty. Every displayed
concatenation is directed, proving bijectivity. For a path ending at y the
same forced order and adjacency hold. Its suffix after u is nonempty and ends
y, while its prefix before v may be empty. ∎

This uses no signed summation and no Fourier inversion. The old selected
vertex survives as an actual endpoint. The new gate's position records the
previous length coordinate without loss.

### Corollary 1.1 (filtered endpoint positions)

Let C be a universal m-filter from Session263: S_m divides K_C, and put
q=K_C(1)=Q_C(2). For ANY pointed B:

- append C after B, then wrap, retaining y. The Hamiltonian paths starting
  at y have equally many possibilities in every residue of the number of
  vertices before v: exactly `q U_B(y)/m`;
- prepend C before B, then wrap. The Hamiltonian paths ending at y have
  equally many possibilities in every residue of the suffix length after u:
  exactly `q V_B(y)/m`.

Proof: use Session263's physical convolution P^+_(B⇒C)=P^+_B K_C,
and its converse P^−_(C⇒B)=K_C P^−_B, then Theorem1. The endpoint counts are
qU or qV, and the source vertex is not replaced. Uniformity of this conditioned
position distribution is all that is asserted.

### Complete initial example

Take the supplied F8, followed ordinally by the singleton y, and wrap once.
The eleven rows, where bit j of row i means i→j, are

    (1514,1516,1441,1508,1391,1472,1412,1296,1024,511,512).

Here y=8, u=9, v=10. Direct path DP, recording the actual position of u, gives

    H = 7623,
    Hamiltonian paths ending at y = 3696,
    suffix-length counts (lengths0,...,11):
      (0,231,462,462,462,462,462,462,462,231,0,0),
    each suffix-length residue mod8 = 462.

The mirrored starting-endpoint example uses `{y}⇒F8` before the wrap and has
the same counts. The complete example pays eleven vertices and the displayed
H; it does not preserve an offset from an arbitrary initial source.

**Boundary:** counting one specified residue still means a position-conditioned
Hamiltonian family. We have NOT built a tournament whose unconditioned H is
exactly that single channel. A scalar channel selector remains open. This
consumer nevertheless turns the old weighted nonspanning paths into ordinary
spanning paths with the old endpoint and a directly visible position statistic.

### Cross-lane composition: an endpoint dipole invisible to the full ordinary profile

The actual six-vertex modules from
`../returns/MULTI_RETURN_CALCULUS.md`, rooted at0, are

    A=(50,44,57,17,34,8), B=(42,52,57,18,33,8).

They share `Q(q)=q^6+12q^4+16q²` and root-deletion profile
`q^5+6q³+2q`, but have `(U,V)=(52,41)` and `(50,43)`. Theorem1 turns this
rooted response dipole into ordinary endpoint counts of the eight-vertex
W(A),W(B), retaining root0. Both have H=349 and root-deleted H=93, while their
starting counts are52 versus50 and ending counts41 versus43.

Equation(4) below proves that their ENTIRE ordinary profiles still agree:

    Q_W(q)=q^8+24q^6+152q^4+172q²,
    Q_(W−0)(q)=q^7+16q^5+58q³+18q.

This supplies two physically different endpoint controllers at identical
unmarked and root-deleted ordinary profile. It neither claims an arbitrary
endpoint count nor identifies that count with the domination excess.

There is also a concrete reason to retain the full coloured sidecar:

    F_A−F_B=4(p_3²−p_1 p_5).

This vanishes at any number of unit-weight slots, but not at the two unequal
weights(1,z). Accordingly

    K_A=(29,50,51,60,51,50,29),
    K_B=(29,54,51,52,51,54,29).

The wrapper's scalar-profile equality is valid, but it does NOT imply equal
length-filter behaviour. The cross-lane example cleanly separates ordinary
profiles, coloured profiles, and retained endpoint responses.

## 2. The full coloured response is the correct wrapping state

For a finite labelled alphabet X=(x_1,...,x_q), define

    F_A(X) = sum_(S_1⊔...⊔S_q=V(A))
                  product_i H(A[S_i]) x_i^|S_i|.

Slots may be empty. Let p_j(X)=sum_i x_i^j. The inherited odd-cycle formula,
applied independently to each slot, gives

    F_A(X) = sum_(odd-cycle packings J)
                 2^|J| p_1(X)^u product_(γ∈J) p_|γ|(X),       (1)

where u is the unused vertex count. Thus F_A is a symmetric function with a
nonnegative integer expansion over partitions into odd parts, with unused
vertices recorded as parts1. `K_A(z)=F_A(1,z)` and `Q_A(q)=F_A(1,...,1)`.
The latter two are different specializations. Neither can silently replace
(1) in an arbitrary future coloured construction.

For a labelled occurrence x_i, X+x_i means adding one more slot carrying the
same weight x_i; X−x_i means deleting that particular slot. These are changes
to path slots, not negative-coloured vertices.

**Theorem 2 (physical wrap operator).** For every finite tournament A,

    F_W(A)(X) = p_1(X)^2 F_A(X)
                + sum_i x_i² [F_A(X+x_i) − F_A(X−x_i)].       (2)

**Proof by path skeletons.** Assign u and v to coloured Hamiltonian path slots.
If they occupy different slots i,j, those slots each have one possible A-run,
and all remaining slots have their usual one A-run. This contributes
`sum_(i≠j) x_i x_j F_A(X)=(p_1²−p_2)F_A(X)`.

If they occupy the same slot i in order v,u, the edge is present. An A-run
may occur before v and after u, independently, giving two slots of weight x_i:
contribution `x_i² F_A(X+x_i)`.

If they occupy that slot in order u,v, the missing edge forces a nonempty
A-run between them, and there are no other A-runs in that slot. Its contribution
is `x_i²[F_A(X)−F_A(X−x_i)]`. This difference removes exactly the configurations
with the forced run empty. Summing the three disjoint cases proves (2). ∎

### Corollary 2.1 (subtraction-free odd-component merging)

For an odd partition λ=(λ_1,...,λ_s), let p_λ=product_j p_λj. The operator in
(2) acts by the completely positive rule

    L(p_λ) = p_1² p_λ
       + 2 sum_(nonempty I⊆{1,...,s}, |I| odd)
                p_(2+sum_(i∈I) λ_i) product_(j∉I) p_λj.     (3)

Repeated equal-sized components are still distinguishable in the sum, so
binomial multiplicities must be retained.

Indeed adding/deleting x_i replaces each p_λj by p_λj±x_i^λj. The difference
retains precisely the odd-sized subcollections, with coefficient2. Multiplying
by x_i² and summing i yields (3). The new part is odd because an odd number
of odd parts were merged and then2 was added.

This is an all-size positive compiler for actual repeated triangle wraps:
keep a weighted table of odd partitions and apply (3) once per physical wrap.
The table is an aggregated count identity, not an asserted canonical bijection
between individual Hamiltonian paths and individual cycle packings.

### Corollary 2.2 (one-variable transfer)

At q equal unit weights,

    Q_W(A)(q)=q² Q_A(q)+q Q_A(q+1)−q Q_A(q−1).               (4)

On a power q^d this acts as

    q^d ↦ q^(d+2)+2 sum_(j odd,1≤j≤d) binom(d,j) q^(d−j+1).

Thus the apparently subtractive difference operator is nonnegative in the
ordinary power basis too. In particular H(W(A))=H(A)+Q_A(2), as inherited.
Equation(4) is the new full-profile extension, not a rebranding of that scalar
identity. Iterated wrapping of a singleton has actual H values

    1, 3, 15, 123, 1479, 23955, 496479   (zero through six wraps).

### A genuine non-preservation boundary

Wrapping the existing F8 does not preserve its universal eight-filter property:

    K_W(F8) = (3927,7548,7545,7752,7584,7512,
               7584,7752,7545,7548,3927).

Its mod-eight sums are unequal. Thus the endpoint consumer is not being
mistaken for an ordinal Fourier multiplier. It retains the useful coordinate
in the particular endpoint-conditioned family from Theorem1, not in the new
unconditioned K.

## 3. A general uniform-interface compiler

Fix a tournament D on k new vertices. Specify for each d∈D whether `d→A`
or `A→d`, uniformly across A. These are actual module substitutions.

To compute K of the combined tournament, enumerate a permutation of D and
one cut, giving the ordered D-skeletons of the two labelled path slots. Each
skeleton has gaps before, between and after its core vertices.

- A gap between d,e permits a nonempty A-run iff `d→A→e`.
- It permits an empty run iff d→e.
- At a boundary, omit the missing condition. The sole gap of an empty
  skeleton permits either an empty or a nonempty A-run.
- If neither choice is allowed, discard the skeleton.

Let q_0,q_1 count allowed A-run gaps in the unweighted and z-weighted paths;
let m_0,m_1 count the ones required to be nonempty because the direct core edge
is missing. Let s be the number of D vertices in the z-weighted path.

**Theorem 3.** This skeleton contributes

    z^s sum_(i=0..m_0,j=0..m_1)
        (−1)^(i+j) binom(m_0,i) binom(m_1,j)
        F_A(1^(q_0−i), z^(q_1−j)).                         (5)

Summing over all `(k+1)k!` skeletons gives the exact K of the literal combined
tournament. Every actual path pair has a unique core skeleton and A-runs;
inclusion–exclusion imposes precisely the mandatory nonempty gaps. This proves
(5), including empty path slots. The computation is finite for each fixed core
and never needs to enumerate subsets of A once its full F-state is known.

Preserved: K and all physical two-slot counts of these interfaces. Lost:
individual paths and vertex addresses inside A. Necessary sidecar for endpoint
variants: the corresponding pointed coloured response, unless the forced
wrap geometry of Theorem1 already identifies it. Equality of F is sufficient
for equality of every such uniform-interface K; no converse is claimed.

## 4. Apex words give an affine constructor, not an exponentially large search

Let A_w consist of an increasing transitive n-vertex base and one apex a,
where w_i=1 means a→i and w_i=0 means i→a. All cycles meet a. Its number of
odd cycles with j base vertices, for even j≥2, is

    c_j(w)=sum_(i<h, w_i=1,w_h=0) binom(h−i−1,j−2).          (6)

The endpoints are forced, and the remaining j−2 base vertices are any
increasing interior subset. Therefore

    F_Aw(X)=p_1^(n+1)+2 sum_(j even≥2)c_j(w) p_(j+1)p_1^(n−j). (7)

In particular, for a unit-weight slots and b slots of weight z,

    F_Aw(1^a,z^b)=(a+bz)^(n+1)
      +2 sum_j c_j(w)(a+b z^(j+1))(a+bz)^(n−j).             (8)

Equations(5) and(8) make K for ANY fixed uniform core affine in the vector
(c_2,c_4,...). Equation(3) is a second exact compiler for iterated wraps. The
word is still required: arbitrary real solutions for the c_j are not physical
constructors. The autocorrelation endpoint counts in(6) are the realization
sidecar.

The known word10010 gives c_2=4,c_4=3; one wrap gives F8. These same equations
now supply every higher word and every repeated wrap without a subset-path
enumeration at each proposal.

### Bounded primitive-sixteen experiment

A primitive-sixteen factor is z^8+1. We tested the exact remainder of K modulo
this polynomial, with a Gray-code word update maintaining the10-pair counts.
No hit occurred in these explicitly bounded universes:

- one wrap: all binary words of lengths12 through22, 8,384,512 words;
- two wraps: lengths10 through22, 8,387,584 words;
- three wraps: lengths8 through22, 8,388,352 words.

This is a stopping signal for those finite grammar ranges, not an all-size
obstruction and not evidence against existence in other tournament families.
The fixed signed64 implementation is safe: rigorous absolute bounds on every
maintained remainder coordinate are respectively 6,740,234,340,236;
3,121,113,964,776,932; and601,510,129,432,076,988, all below2^63.
The same scanner has a positive mod-eight control at word length5; it
recovers the supplied10010 primitive-eight word and its other successful
words. Its input/output are `scanner_positive8_input.txt` and
`scanner_positive8.out`.

No further scan was pursued after the positive operator and endpoint consumer
were obtained. The constructive output is the all-size compiler and consumer.

## 5. Unrelated-source bridges and retained sidecars

Sources were sampled with `git ls-tree -r --name-only HEAD` and read with
`git show HEAD:<exact path>` at inherited main8aec2264f. They are inspiration,
not newly audited dependencies:

1. **Cyclotomic factorial selectors:**
   `01-canon/theorems/THM-4001-cyclotomic-factorial-all-arity-decoder.md`.
   Its coordinatewise projector versus diagonal average suggested keeping
   labelled path slots instead of one length scalar. The map is the exact
   coloured partition response F; the preserved object is an actual path
   decomposition. We do not import its formal Reynolds average as a tournament.
2. **Reciprocal/jet transfer:**
   `01-canon/theorems/THM-3192-reciprocal-coefficient-jet-transfer-and-z-adic-pluecker-return.md`.
   Its lesson of changing state before iterating becomes the lift K→F and
   positive component operator(3). No Pluecker or p-adic conclusion is asserted.
3. **2-adic renormalized charts:**
   `01-canon/theorems/THM-4082-mahler-renormalized-linear-chart-and-exact-bit-defect.md`.
   It motivated separating the coordinate being observed from the transformed
   object's global invariant. Here the endpoint-conditioned position channels
   survive even though K_W(F8) loses the old filter zeros. This is a structural
   analogy, with no analytic theorem transferred.
4. **Tournament automata and codes:**
   `01-canon/theorems/THM-4857-the-gf2-automaton-of-a-tournament-classifies-the-merged-metagraph-exactly-through-order-five-circulants-are-additive-cellular-automata-and-paley-tournaments-give-quadratic-residue-codes.md`.
   The finite code viewpoint suggested the word realization constraint and
   the10-pair autocorrelation in(6). Our weighted integer cycle spectrum is
   not the GF(2) adjacency code; losing the actual word loses realizability.

The inherited mathematical dependencies are Session263's literal convolution,
explicit F8, and wrap geometry, plus the odd-cycle count formula
`01-canon/theorems/THM-002-ocf.md`. The full coloured operator has its own proof
above. The concept board at close-out is: endpoint positions; full coloured
state; odd-component merging; finite interfaces; binary-word realization.

## 6. Reproduction and scope

From the repository root:

    python3 05-knowledge/results/session264_constructive_20261111/filters/controls.py

The190 new-operation gates cover31 word sources at zero/one/two wraps,
20 arbitrary sources of orders3–7 with the same three depths, the explicit
retained-endpoint consumer in both directions, the non-preservation
hostile, and the two cross-lane endpoint dipoles. Direct path DP and the component/interface compilers take independent
counting routes. `controls.out` contains the exact eleven-vertex rows and counts.

`coalescence.py` implements(3); `interface.py` implements(5),(8);
`explore.py` supplies(6) and a direct one-wrap cycle expansion. The bounded
scan source is `primitive_scan.cpp`; the three input/output pairs are retained.
Use `python3 interface.py r lo hi output_path` to regenerate an input for
r wraps, compile the scanner with `g++ -O3`, then pipe that input into it.

No selected scalar projector, arbitrary dyadic constructor, or Chan–Pak
coverage theorem follows. The new all-size ingredient is a physical positive
transfer on cycle-cover components, and the new service is a retained endpoint
whose gate positions consume the already-constructed length channels.
