# Node-coloured defects, needles, Radon data, and the typed metagraph

**Status:** exact synthesis of proved identities and finite audits, followed by
explicit conjecture targets. This is not a proof of LRC(14). The standard
fourteen-runner case remains open.

The recurring failure in this part of the repository is not a shortage of
statistics. It is a type error. A tournament class, a fixed-path tiling, a
complement line, a factorization witness, a path realizer, and an LRC safe
component answer different questions. Equal cardinalities between two of
these sets do not supply a map, and a map that is adequate for one operation
need not be a congruence for the next operation.

The useful replacement is a **typed, operation-indexed metagraph**. It retains
literal witnesses until a minimization theorem proves that a field is
redundant. The node-coloured defect algebra, the Hamiltonian-path Walsh stalk,
the path-holonomy law, and the finite-field Radon proposal are four shadows of
that one rule.

## 1. The object has several sorts

At size `n`, put

```text
S_n = {(a,b): 1 <= b < a <= n, a-b >= 2},
m_n = |S_n| = binom(n-1,2),
X_n = F_2^(S_n).
```

The fixed path is `n -> ... -> 1`. A point of `X_n` is a literal marked-path
tiling and hence a labelled tournament with that Hamiltonian path. The maps
already proved in THM-781, THM-796, and THM-830 form the span

```text
                         q_n
                 X_n ----------> M_n
                  |
                  | orbit by kappa
                  v
                 L_n = X_n/<kappa> ----> Sym^2(M_n).
```

Here `M_n` is the converse-merged tournament-node set and `L_n` is the set of
literal complement lines. The endpoint map records a node loop or an unordered
pair of nodes. It also carries a blue/black colour. These are three sorts, not
three names for one graph.

THM-830 makes the incidence numerical. If `Phi_n(G)` is the merged node of the
gap tournament `G`, then

```text
K_n^c(u,v) = #{G : Phi_n(G)=u, Phi_n(G^op)=v, colour(G)=c}.       (1.1)
```

For `u != v`, this is the literal multiplicity of the coloured node edge. For
`u=v`, the loop multiplicity is `K_n^c(u,u)/2`. Its row sums recover the
coloured tiling fibres. Thus node adjacency, tiling incidence, and line
multiplicity are all pushforwards of one literal kernel, but none determines
the other two without the fibre weights.

There are at least four further witness sorts:

```text
G_n       defect-groupoid arrows (x,delta),
F_n       two-arrow factorizations g = g_2 o g_1,
R_n       direct self-line realizers (T,P,sigma),
O_n       owner-labelled LRC wall, gap, and proof-obligation packets.
```

The first three are purely tournament-theoretic. The fourth is where the LRC
metric lives. No result presently gives a lossless map from `O_n` to a bare
tournament node.

## 2. The node-coloured defect algebra

Staircase reflection has `f` fixed tile coordinates and `r` transposed pairs.
Write

```text
h=f+r,                 m=h+r,
x in F_2^h,            delta in F_2^r.
```

The arrow `(x,delta)` has target `x+i(delta)`, where `i` inserts the defect in
the paired coordinates, and composition is

```text
(x+i(delta_1),delta_2) o (x,delta_1)
    = (x,delta_1+delta_2).                                  (2.1)
```

This is an elementary abelian pair groupoid. Reflection is arrow inversion.
Blue arrows have `delta=0`; black arrows have `delta != 0`. A seed colour may
retain the full defect and the merged node of the decoded endpoint. It is
already more informative than a node colour or a blue/black bit.

For an arrow colour `c`, define the ordered factor tensor

```text
p_ab(g) = #{(g_1,g_2): g=g_2 o g_1, c(g_1)=a, c(g_2)=b}.        (2.2)
```

The recursive refinement

```text
c_(s+1)(g) = (c_s(g), (p_ab^s(g))_(a,b))                       (2.3)
```

is the coarsest refinement of the seed for which all composition structure
constants are cellwise constant. Before quotienting by inversion, the Jordan
version retains `p_ab+p_ba`. Its symmetric-closure interpretation requires
the input colour to be inversion-invariant; omitting this hypothesis would be
false in general.

THM-851 gives the exact first-step census:

```text
n                 3    4     5      6       7
arrows            2    8    64   1,024  32,768
seed cells         2    4    23     255   7,926
ordered cells      2    6    64   1,024  32,768
Jordan cells       2    5    40     544  16,640
reflection orbits  2    6    40     544  16,640.              (2.4)
```

Thus ordered factor data reconstructs every arrow at `n=3` and `5<=n<=7`,
while Jordan data reconstructs every reflection orbit there. Both fail in the
small exceptional chart `n=4`. This exception is useful: even an exact algebra
can have a low-size accidental kernel, so an all-size claim needs a proof, not
an interpolated table.

The crucial information is **which defect carries which factor pair**. If the
midpoint defect labels are erased and only the unlabelled multiset of
intermediate node pairs is retained, the deck misses `1,1,36` distinctions at
`n=5,6,7`. A node-coloured algebra is therefore a labelled incidence tensor,
not merely a histogram of adjacent node colours.

For fixed output `(x,gamma)`, the full deck can be written as the function

```text
alpha |-> {C(x,alpha), C(x+i(alpha),alpha+gamma)},
alpha in F_2^r.                                             (2.5)
```

Because `C` retains the exact defects, the factor colours recover the index
`alpha`; the multiset presentation and this indexed function are equivalent.
Its Walsh transform in `alpha` is therefore lossless. The unlabelled midpoint
deck is its node-pair population shadow. A sharp compression problem is to
find the least Walsh sector of (2.5) that reconstructs reflection orbits and
remains closed under A/B/C and the THM-853 return semigroup.

## 3. The H drift is a second operation algebra

Let `H(T)` be the number of directed Hamiltonian paths and let `T^e` denote
one arc flip. Put `M=binom(n,2)` and

```text
b_H(T) = (1/M) sum_e (H(T^e)-H(T)).                          (3.1)
```

THM-848 decomposes `H=sum_r H_(2r)` into even Walsh levels supported on
vertex-disjoint unions of positive even paths. The exact drift is

```text
b_H = -(4/M) sum_r r H_(2r).                                (3.2)
```

If `U_1` counts the Hamiltonian needles with one boundary defect and `K` is
the odd-cycle-forest partition function, then

```text
sum_e Delta_e H = U_1-(n-1)H,
U_1 = 2K-(n+1)H,
K = nH + (M/2)b_H.                                         (3.3)
```

So the H-drift functional-form question has a sharp answer: one step needs
the degree-Euler image `K=(n-D)H`, equivalently the length-weighted Walsh
stalk. `H` alone is insufficient. Iterated means require the whole Krylov
stalk generated by the Walsh levels; even the full first-moment stalk need
not determine diffusion.

The energy of level `2r` is

```text
E_(2r)/mu_n^2
 = [w^r] ((1+w)/(1-w))^(n-2r) / (n)_(2r),
mu_n = n!/2^(n-1).                                         (3.4)
```

At `n=8`, (3.4) gives `E_6/mu_8^2=1/1680` and total variance ratio
`131/560`; the former formula first fails there. This is another next-size
guardrail against pattern extrapolation.

The multiset

```text
Grad_H(T) = multiset_e (H(T^e)-H(T))                        (3.5)
```

classifies every converse-merged node through `n=8` in the exact atlases:
`3,10,34,272,3528` cells at `n=4,...,8`. This is finite evidence, not an
all-size theorem. It is nevertheless a strong candidate node stalk because it
records all one-step directional responses without naming arbitrary labels.

Equations (2.2) and (3.5) are the same recursive idea at different operation
alphabets: retain the multiset of coloured outcomes of every legal local
factor or flip. Their common abstraction is the operation profile

```text
Prof_O(x) = multiset_(o in O) (o, colour(o x)).             (3.6)
```

## 4. Self-lines require path holonomy

Let `kappa_P` flip all pairs off the marked Hamiltonian path `P`. The black
quasi-fixed endpoint set is stable under the Klein four generated by full
flip `kappa` and staircase reflection `sigma`, and its orbits are free. Hence

```text
4 divides Q_black(n),             2 divides selfK(n).       (4.1)
```

The finite endpoint counts are

```text
n                    5    6    7     8
Q_black(n)           8   12   88   404
selfK(n)             4    6   44   202
Klein-four orbits    2    3   22   101
SC(n)                8   12   88   176.                    (4.2)
```

Therefore `Q_black(n)=SC(n)` at `n=5,6,7` is a small-size coincidence, not an
all-size self-line law. There can be no all-size two-to-one bijection with
self-converse classes in the proposed form.

The surviving all-size law is witness-level. If `pi` realizes
`pi.T=kappa_P T`, and `L_pi` is the linear relabelling action on pair masks,
then

```text
pi^k.T = T + [k odd]1 + sum_(j<k) L_pi^j p_0,              (4.3)
N_(L_pi) p_0 = [ord(pi) odd]1.                             (4.4)
```

Thus odd witness order requires a mod-2 decomposition of `K_n` into the
translates of one Hamiltonian path, and is possible only for
`n=1,2 (mod 4)`.

Most importantly, relabelling does not commute with a fixed path flip:

```text
sigma kappa_P = kappa_(sigma P) sigma,
sigma^2 T = Flip_(E(P) triangle E(sigma P)) T.             (4.5)
```

The square is an automorphism only when the two path edge sets agree. Exact
search finds no such direct black realizer for `n=5,...,8`. The quotient node
forgets precisely the obstruction it was asked to decide. The minimal stalk
is the marked pair `(P,sigma P)`, or at least its symmetric-difference mask,
not just the permutation cycle type.

In general the realizer-pair mass is

```text
sum_(t in Q_n^K) |Aut(T_t)|,                               (4.6)
```

not `2 selfK(n)`. Equality in the audited range comes from asymmetric
endpoints. This is the correct groupoid weighting behind any future counting
identity.

## 5. The three signed recursions are charts, not scalar laws

The full staircase has the exact `B3` Cech cover

```text
A+B+C-D-E-F+G.                                             (5.1)
```

THM-801 proves literal gluing, while THM-830 identifies every face as a gap-
tournament operation. With address variables `X,Y,Z`, the typed enumerator is

```text
T_n=h_(n-3)(X,Y,Z),
T_n=(X+Y+Z)T_(n-1)
    -(XY+XZ+YZ)T_(n-2)+XYZ T_(n-3).                        (5.2)
```

Setting `X=Y=Z=1` retains the cardinality and destroys the direction that
owned each coordinate. That direction is exactly what deletion, gap
contraction, and continuation use.

There is an exact invariant-ring form of this loss. Reflection exchanges
`X,Y`, and

```text
Q[X,Y,Z] = Q[s=X+Y,p=XY,Z]
           direct-sum (X-Y) Q[s,p,Z].                       (5.3)
```

In homogeneous degree `n-3`, the two summands have dimensions `h_n` and
`h_(n-1)`, so `m_n=h_n+h_(n-1)`. The half chart is the Reynolds-invariant
summand. What it destroys is the entire anti-invariant converse-current
module, not one unstructured sign bit.

The even half chart is the `B2` fold

```text
A+B-C.                                                     (5.4)
```

Its scalar recurrence is valid, but the fixed reflection line, shortcut seam,
and endpoint coupling are sidecars. THM-830 shows that the middle `B` term is
a tournament on cuts rather than deletion of one fixed original vertex.

The prompt-order odd word is

```text
A+B-C+D-E-F+G.                                             (5.5)
```

Its actual Venn chart has corners `A,D,B`, edges `A+B-C`, `A+D-E`, and
`B+D-F`, and center

```text
A+B+D-C-E-F+G.                                             (5.6)
```

The slots `C` and `D` have the same scalar size in the relevant count, so
they cancel after forgetting their roles. They are geometrically different:
one is an overlap and one is a corner. Any recursive carrier that cancels
them before transport is too small.

Algebraically, in odd homogeneous degree the invariant piece is covered by
the ideals `(s)` and `(Z)`, yielding the two-generator fold. In even degree,
the pure monomial `p^((n-3)/2)` is not covered and forces the extra corner.
This is the invariant-ring reason for the parity change and for retaining the
equal-size `C,D` roles separately.

The natural arithmetic companion is not another scalar sequence but the
product ledger

```text
Div(D) x B_r.                                              (5.7)
```

The divisor coordinate retains primitive period capacity `phi(d)` and CRT
labels; the Boolean coordinate retains which face, overlap, or far packet
owns the term. Mobius inversion is legal on each typed axis. Multiplying the
two only after their labels are retained prevents a cancellation in one axis
from erasing an obstruction in the other.

## 6. Three meanings of a needle

The repository now has three exact needle constructions. They should be
related by their incidence algebra, not identified literally.

1. A **Hamiltonian needle** is a directed vertex ordering in the tournament
   arc cube. Its Walsh boundary is an even path forest, and arc flips give
   the gradient stalk (3.5).
2. A **Farey needle** is a threshold interval moving through a rational gap.
   THM-841 shows the violation ladder has exact breakpoint profiles but not a
   homogeneous toothpick recursion.
3. An **affine finite-field needle** is a line of `F_7^2`. Its line sums give
   a discrete Radon transform.

For `v` in the eight directions `P^1(F_7)`, define

```text
R_v f(s) = sum_(x: v dot x=s) f(x),        s in F_7.       (6.1)
```

For any `d` distinct directions, over characteristic zero,

```text
rank(R_D)=1+6d,                 dim ker(R_D)=48-6d.        (6.2)
```

Proof: Fourier transform in `s` reveals `f-hat(a v)` for `a in F_7`. Every
direction supplies its six nonzero frequencies, while all directions share
the zero frequency. The eight projective frequency lines partition the 48
nonzero frequencies. Hence all eight directions are necessary and sufficient
to reconstruct an arbitrary scalar function on `F_7^2`.

The integer address simplex has an additional odometer fibre over this plane.
For

```text
S_r={(alpha,beta,gamma)>=0: alpha+beta+gamma=r},
```

choose canonical residues `xi=(a,b,c) in {0,...,6}^3` and put
`L_xi=(r-a-b-c)/7`. Then

```text
S_r ~= disjoint-union_(xi in H_(r mod 7)) Comp_3(L_xi),
(alpha,beta,gamma)=xi+7(i,j,k).                            (6.3)
```

An A/B/C increment advances the corresponding residue odometer; a `6->0`
wrap advances its carry coordinate. Thus the carry fibre recursively carries
the same B3 face structure.

At the `n=14` depth `r=11`, 33 residue fibres are singletons, 15 are triples,
and `(6,6,6)` is absent. Scalar residue pushforward has rank 48 and kernel
dimension

```text
30=15(3-1).                                                (6.4)
```

The full eight-direction Radon map has rank 49 on the abstract plane but only
rank 48 after composition with the 78-address pushforward. It reconstructs
residue sums, not the two independent A/B/C carry contrasts inside each
triple fibre. The three direction orbits under face triality have sizes
`3+3+2`; the chi charge lies in the two-direction orbit.

With no-carry and A/B/C carry treated as four channels, one common
reflection-invariant six-pencil deck has ranks `(33,15,15,15)` and total rank
78. Five pencils cannot reconstruct the 33-point no-carry channel. The exact
audit finds this full-rank six-set unique under true endpoint reflection; the
formal-triality six-set has total rank 77. Thus full triality and minimal
tomography are competing gauges, not interchangeable symmetries.

This supplies the precise Kakeya warning. One line in every direction is a
geometric hitting object; all parallel line sums in every direction are a
tomographic reconstruction object. A single chi-seven pencil has rank seven
and a 42-dimensional kernel. Seven directions still leave a six-dimensional
kernel. A proof that samples one attractive pencil cannot claim to have
preserved the plane.

The j=4 flood audit provides a structured, non-generic function. Its 21 floods
are indexed by the edges of `K_7`; Fano triples organize those edges, and the
negative masks form one Fano line. But only the identity in `GL(3,2)` preserves
the full exact numerical flood data. The Fano plane is therefore an address
system, not yet a numerical quotient. The one-sided endpoint needles visit at
most four masks and at most two masks on the negative line, so they do not
close a local chi-seven sweep.

The flood observables also occupy different edge modules. For edge `{a,b}`,

```text
r_(a,b)=x_a+x_b,       x=(16,16,14,14,10,14,6),            (6.6)
```

so `r` has orbit-span rank 7. The exact `m` and `V1` fields each have full
orbit-span rank 21. Point-star and Fano-triangle rows together span only 13
dimensions, leaving an eight-dimensional invisible edge field, and both `m`
and `V1` have nonzero curl there. A Fano transform may compress `r`; it cannot
compress the other two without retaining their residual edge cocycle.

The next experiment should replace the scalar `f` in (6.1) by a value in the
free module generated by

```text
(node, exact defect, factor colour, carry role, seam, owner, metric obligation). (6.5)
```

One then computes ranks only after applying declared forgetful functors. The
scalar plane-rank theorem is a ceiling after residue aggregation, not a
certificate for the integer address tensor. A structured image may need fewer
directions or carry contrasts, but that must follow from a proved restriction
on the valued function, not from the Fano picture alone. Radon acts on Parikh
vectors; ordered j=4 peel chronology remains a separate history stalk until a
commutation theorem proves otherwise.

## 7. The ladder has a source term

The violation ladder's dyadic breakpoint sets obey the exact relation

```text
D_(2p) = D_p union (1/2)D_p
         union {1/m : p < m <= 2p, m odd}.                 (7.1)
```

The last set is a birth term. It is the information that a pure rescaling
picture drops. The corresponding cumulative counts match the classical
toothpick totals at `1,2,4,8`, giving `1,3,11,43`, and fail at `16`, where the
ladder has `159` rather than `171`.

This is not a disappointing analogy. It identifies the exact recursive datum
that must be preserved: inherited breakpoints plus newly primitive odd
denominators. The same pattern occurs in (5.7), where divisor transport alone
does not create the new primitive-period packets. Recursion is generally
**copy plus source**, not exact self-similarity.

There is nevertheless an exact self-affine subchannel. If
`W_a^H(s)=sum_(j=0)^H zeta^(2^j) 1_(s<a/2^j)`, then

```text
W_a^H(s)=zeta 1_(s<a)+zeta^2 1_(s<a/2)+zeta^4 1_(s<a/4)
         +W_a^(H-3)(8s).                                  (7.2)
```

The full ladder is therefore a self-similar chi-labelled stalk plus the
non-self-similar primitive-odd source in (7.1). This is the same odometer
pattern as the residue-plane base plus carry-simplex stalk in (6.3).

## 8. What the four-dimensional LRC object is

There are two exact uses of "four-dimensional" in the current frontier, and
they must remain distinct.

First, after the `f<=3` branch, fix one of the `binom(14,9)=2002` nine-speed
cores. The four ordered far speeds form a rank-four semilinear cone chart. It
has radial and transverse noncompact directions. A point must carry an
owner-coloured threshold loop, rational gaps, marked relations, residue
addresses, and the next proof action. The bare `Z^4` coordinate does not decide
the metric predicate.

Second, an affine family `V_i(c)=c a_i+r_i` is the integer-slope fibre of

```text
Phi_(A,R)(u,t)=min_i ||a_i u+r_i t||,
X_(A,R)={(u,t,c,lambda): u=ct, Phi_(A,R)(u,t)>=lambda}.     (8.1)
```

This is a mixed two-torus, integer-slope, clearance suspension. For fixed `c`
it has only two continuous coordinates. It is not a second four-manifold and
does not turn the global LRC14 moduli space into a four-dimensional space.

The tournament metagraph can serve as a **base atlas** for combinatorial
types of these fibres. The metric and owner data form a stalk over it. The
target is not a map from all LRC packets to one node; it is a stratified
pullback in which a node change transports the stalk and every wall crossing
updates its owner and clearance.

At `n=14`, one literal tiling has 78 bits. Reflection rewrites these as

```text
fixed bits f=6, paired objects r=36, h=f+r=42,
defect bits delta in F_2^36,              h+r=78.          (8.2)
```

The chi-seven grading rewrites the same 78 bits as one neutral block of size
12 and six nonzero charge blocks of size 11. These are two coordinate systems
on the same bits, not 156 independent coordinates. Other necessary views are

```text
six nonconstant Walsh levels H_2,...,H_12,
91 directional H derivatives,
the 2^36 factor-defect deck at each object,
marked path holonomy,
owner/wall/carry/clearance data.                            (8.3)
```

The apparent size is a reason to seek transforms and coherent algebras, not a
license to average away labels.

## 9. The universal recursive definition

Let `O` be a declared alphabet of legal operations: face restrictions,
internal deletion with seam repair, inverse lift, CF return, arc flip,
reflection, path-relative complement, factorization, or an LRC observer
update. Let `P` be the declared terminal predicate or valued obligation.
Define

```text
x ==_(O,P) y
  iff P(w x)=P(w y) for every legal operation word w.      (9.1)
```

This is the operation-indexed Nerode equivalence. It is the smallest quotient
that is simultaneously a congruence for the chosen operations and sufficient
for the chosen task. Changing `O` or `P` changes the quotient. Static
injectivity, one-step lumpability, and LRC metric closure are three different
claims.

Across sizes, the literal objects and face maps form a presheaf on the
operation category. A proposed statistic is safe exactly when it is a natural
quotient of that presheaf for the declared subcategory. A failure of
naturality is not noise; its kernel pair is the next defect stalk.

This reframes several historical results uniformly:

```text
THM-801: literal B3 face descent is functorial.
THM-813: projected node/edge CF descent can fail functoriality.
THM-840: a kernel may be a congruence for insertion but not deletion.
THM-848: H requires its Walsh/Krylov operation stalk.
THM-851: node-defect colours close under factor composition through n=7.
THM-854: a fixed-path quotient drops the square holonomy.
THM-841: dyadic scaling needs a primitive-denominator source term.
```

## 10. Formula ledger for future tournament transfers

The following small formulas are proved or exact finite targets worth keeping
together.

```text
tiling dimension:       m_n=binom(n-1,2)
half-object dimension:  h_n=floor((n-1)^2/4)
reflection defect rank: r_n=floor((n-2)^2/4)
full reconstruction:    m_n=h_n+r_n
line count:             |L_n|=2^(m_n-1)                 (m_n>0)
cyclic triangles:       C3=binom(n,3)-sum_i binom(d_i,2)
mean Hamilton paths:    mu_n=n!/2^(n-1)
H drift:                b_H=2(K-nH)/binom(n,2)
score OU law:           E[Delta x|T]=-8(x-n(n-1))/(n(n-1))
defect composition:     delta_total=delta_1+delta_2
black endpoint parity:  4 | Q_black, 2 | selfK
witness norm:            N_(L_pi)p_0=[ord(pi) odd]1
Radon rank:              rank R_D=1+6|D|
B3 enumerator:           T_n=h_(n-3)(X,Y,Z)
deletion ownership:      (n-a)+(a-b-1)+(b-1)=n-2.
```

These formulas do not imply one another. They identify compatible coordinates
that can be joined and then subjected to (9.1).

## 11. Exact next experiments

1. Compute THM-851's factor closure at `n=8` without materializing the roughly
   billion naive factor incidences. Use defect translation/convolution and
   test whether the ordered closure reconstructs all arrows and whether the
   Jordan closure reconstructs all reflection orbits.
2. Join the `n=7,8` factor-defect colours with `Grad_H`. Determine whether
   gradient cells are already factor-coherent, or record the first crossed
   pair of fibres. This directly compares composition and flip operation
   algebras.
3. Extend the black self-line census to odd `n=9`, retaining witness cycle
   type, automorphism weight, and path-holonomy module. Test (4.4), not the
   refuted SC count.
4. Build the full eight-direction `F_7^2` valued Radon deck for the 21 j=4
   floods, with the A/B/C carry simplex attached before residue pushforward.
   Compute scalar, node-valued, defect-valued, and future-refinement ranks
   separately; compare each with the exact `48+30` address split.
5. Search for a natural role map `tau` on the 21 floods that commutes with
   subtree restriction and reflection. Fano incidence alone is not such a
   map because the exact numerical stabilizer is trivial.
6. Attach owner-event cocycles and clearance margins to the four-far scale
   action. Test whether the resulting packet is a congruence under peel and
   insertion, rather than only under static residue classification.
7. Finish THM-853's operation-semigroup census and minimize the joined
   `(node, line, defect, factor, gradient, holonomy)` state under the actual
   multiword language. Compare it with each marginal before adding LRC data.

## 12. Assumption challenge

Tournament Analysis need not put runners or tournament arcs at its vertices.
Depending on the question, legitimate vertices include:

```text
gap positions, face directions, factor colours, defect sectors,
Hamiltonian needles, path masks, affine-line directions,
owner events, Fourier levels, and proof obligations.
```

For each choice, the pairwise observable and gauge must be declared, followed
by the tie Hamiltonian path and the standard fingerprint. More importantly,
the analysis must state which LRC predicate its quotient preserves and what it
destroys. A transitive carrier tournament can still lose higher overlap,
absolute metric, or continuation truth.

The most promising recursive object is therefore not a larger untyped graph.
It is the smallest natural, operation-closed quotient of the literal typed
incidence system for the terminal proof obligation. The computations above
show several of its stalks exactly. They also show, with unusual consistency,
where every smaller description breaks.
