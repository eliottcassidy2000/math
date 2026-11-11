# Rooted seam fans have a three-mode whole-profile response

**Status: PROVED IN NOTE.** New structural extension of Session268's one-seam
identity and Session273's resource exposure. The scalar retained-root dyadic
law is a recovered specialization, not claimed as new. Chan–Pak remains OPEN.

## 1. A seam never has to leave its fixed-path tiling chart

Let T be a nonempty tournament, a a named vertex, and P any named directed
Hamilton path. Add b with b→a and u→b for all old u≠a. This is exactly
Session268's one-seam operation T ⋈_(a,b) {b}, or Session273's S(T,a).

Insert b immediately **before a** in P. The result is Hamiltonian: the new
outgoing adjacency is b→a; if a had predecessor u, u≠a and u→b. All other
adjacencies are untouched. This holds even if a was first or last.

Repeat at the same root k times, naming new vertices b1,...,bk in birth order.
The resulting tournament F_k(T,a) has b_i→b_j for i<j and every b_i→a;
every old u≠a beats all b_i. Its canonical named path is obtained by replacing
...a... with ...b1,...,bk,a.... This is a genuine fixed-Hamilton-path framed
construction at every size, not merely existence of some Hamilton path.

In a tiling chart whose consecutive path arcs are fixed, the birth operation
adds the corresponding row/column of free arc coordinates and moves the named
root position by one. It does not flip an old free tile. Any old chart flip
away from the compulsory path arcs commutes with the extension, since the new
arc rule depends only on the root label. Thus this is an injective cubical
map between the two labelled fixed-path charts; it is not an assertion about
injectivity after quotienting by ordinary tournament isomorphism.

## 2. Retain a colour-marked seam, not an entire exponentially growing fan

For a positive integer t, Q_T(t) counts colourings of all vertices by t named
colours, together with one directed Hamilton path on each nonempty colour
class. Empty colour classes contribute one. Equivalently this is the usual
falling-factorial path-cover polynomial sum p_d(T)(t)_d, where p_d counts
unordered d-path covers.

The root's colour is determined by the colouring. Define:

* E_T,a(t): the same objects, with the path of the root's colour ending at a.
* L_T,a(t): replace that path by an ordered pair (U,V), where V starts at a,
  V is nonempty, U may be empty, and together they cover the root colour.
* D_T,a(t): the subfamily of L with V the singleton a.

The other colours still each carry one path. These are polynomials: there are
finitely many covers and each object has finitely many used colour classes.
Moreover D_T,a(t)=t Q_(T−a)(t), including a singleton T and Q_empty(t)=1.
One chooses the root's colour freely, then any coloured path cover of T−a;
its path in that colour becomes U. This remains true when U is empty.

**Theorem (three-mode seam-fan response).** For all k≥0 and every positive
integer t,

    Q_Fk(t) = (L−D)(t+1)^k
            +(Q−E−L+2D)t^k
            +(E−D)(t−1)^k.                       (1)

The identity is polynomial in t, hence extends as a polynomial identity.
At k=0 interpret all zero powers as 1, including 0^0. In particular, the
whole path-cover response of arbitrarily large fans is generated from four
polynomials of the original rooted object, not by re-enumerating the fans.

**Proof.** In any coloured path cover of F_k, a new vertex can only be followed
by another new vertex or by a. Thus in a nonroot colour all new vertices form
one terminal increasing block. In the root colour they can form a block just
before a and another terminal block. Removing those blocks leaves one old
path in every nonroot colour and either one or two old paths in the root
colour. There are no further run types.

First suppose the block immediately before a is empty. The old object is
counted by Q. If its root-colour path does not end a, each of the k new
vertices independently chooses one of t terminal blocks, giving t^k. If
that old path ends a, its terminal block is forbidden because a loses to all
new vertices; there are (t−1)^k choices. The contribution is
(Q−E)t^k+E(t−1)^k.

If the block immediately before a is nonempty, deleting it leaves (U,V)
counted by L, with no requirement that U's last vertex point to a. If V is
not singleton a, its terminal vertex is not a and may enter the root's
terminal block. There are t+1 places for each new vertex: t terminal blocks
and the distinguished pre-a block. Requiring the latter nonempty gives
(t+1)^k−t^k. If V={a}, its terminal block is forbidden and the count is
 t^k−(t−1)^k. The total contribution is
(L−D)[(t+1)^k−t^k]+D[t^k−(t−1)^k].

Each assignment produces a unique increasing block in each position. Blocks
concatenate legally: old vertices other than a beat new vertices, and new
vertices beat a. These constructions and deletion are inverse. Adding the
two contributions proves (1). ∎

Consequently, for each formal t, Q_k satisfies the exact recurrence

    Q_(k+3)=3t Q_(k+2)−(3t²−1)Q_(k+1)+(t³−t)Q_k,

and

    Σ_(k≥0) Q_k z^k = (L−D)/(1−(t+1)z)
                      +(Q−E−L+2D)/(1−tz)
                      +(E−D)/(1−(t−1)z).

These are identities of rational formal generating functions, not numerical
fits. Some modes can vanish or coincide after specializing t. No minimal-order
claim is made. The coefficients can have either sign; they are a decomposition
of a positive combinatorial count, not separate tournament counts.

## 2a. Three experiments completely identify this continuation response

Put A=L−D, C=E−D and B=Q−A−C. The physically relevant packet for
this operation language is (Q,A,C), not the four individual polynomials.
If Q_0,Q_1,Q_2 are the whole profiles before extension and after one and two
seam births, set

    U=Q_1−t Q_0,
    V=Q_2−2t Q_1+t² Q_0.

Then A=(V+U)/2, C=(V−U)/2 and B=Q_0−V. Consequently two rooted
inputs have identical full profiles after EVERY fan length if and only if
they agree after lengths 0,1,2. The calculation follows immediately from the
three distinct formal eigenvalues t−1,t,t+1 and holds as a polynomial
identity, even at t=1 where the zero eigenvalue needs k=0 to be observed.

This is a finite complete continuation test for the particular seam-fan
language, independent of input size. It does not say these three profiles
classify rooted tournaments or arbitrary contexts. In particular a common
change in E,L,D is invisible when Q stays fixed. This residual kernel is the
next hidden interface coordinate that a different seam operation could expose.

## 3. Scalar boundary, first cases, and what renews

At t=1 write H=Q(1), e=E(1), L=L(1), d=H(T−a)=D(1).
For k≥1, (1) becomes the familiar retained-root form

    H(F_k)=H−e+d+(2^k−1)(L−d).

The full response formula explains its dyadic term as allocation between two
permitted rooted blocks. It also displays the lost t−1 mode, invisible to
scalar H once k≥1. Thus even the entire scalar fan sequence cannot in general
recover all four initial rooted response polynomials.

The rooted two-cover quantity evolves as L(F_k,a)=d+2^k(L−d), and e(F_k,a)=d
for k≥1. For the first identity assign every new chain vertex to the prefix
or suffix of a rooted two-cover; the suffix can accept a terminal chain unless
it is singleton a. The second follows because a path ending a must traverse
all new vertices in one final increasing block after a Hamilton path of T−a.
The latter also covers T={a}. For k=0, retain the original e rather than d.

A singleton has H=e=L=d=1 and all fans are transitive, H=1. A rooted directed
triangle has H=3,e=d=1,L=3, so H(F_k)=2^(k+1)+1 for k≥0:
3,5,9,17,33,65,129. These are exact initial constructions, not coverage.

The process renews a usable named seam indefinitely and preserves its fixed
path, but its excess L−d doubles rigidly. It is **interface renewal without
independent target selection**. This is a precise distinction that Session273's
resource-transport request needs: inexhaustible seam supply alone does not
supply arbitrary odd increments.

## 4. Relation to past mathematics and the next useful question

* Session268 `automata/PROOF.md`, Theorem1, supplies the inherited singleton
  seam identity H'=H−e+L. The present proof lifts it simultaneously to every
  colour count and every seam length, with a framed geometric implementation.
* Session273 `interfaces/PROOF.md`, TheoremD, supplies invisible full-profile
  gauge transport and one seam exposure. Formula(1) tells exactly what its
  entire repeated-seam experiment reads: Q,E,L,D, with three explicit modes.
  Equal unmarked Q is insufficient; equal four-polynomial signature is
  sufficient for all fan lengths, before any uniform-substitution context.
* The mechanism resembles finite-state transfer matrices and the repo's
  cellular-automaton boundary-memory work, but here the finite state is proved
  by the restriction that new vertices can re-enter the old graph only at a.
  The run language, not an analogy alone, yields the three modes.
* The fixed-path chart is a stable coordinate system while the colour-path
  observers change. Therefore tiling coordinates and path-cover coordinates
  are complementary: one certifies physical admissibility, the other the
  all-size response.

**OPEN constructive obligation.** Find a second physically framed operation
that changes the three-mode coefficient packet in an independently selectable
way while returning it to an admissible seam-fan class. A proof of closure of
such a packet family, together with exact integer landing, would advance
Chan–Pak. Neither the recurrence nor the canonical chart lift proves that
closure, and neither licenses subtraction of response coefficients as graphs.

Run `python3 check.py` here for exhaustive small controls of (1); the declared
universe is all labelled tournaments of orders1..4, every root, t=1..3 and
k=0..3. These are controls of the new identity, not an enlarged H-spectrum
census. The all-order proof is the path-block bijection above.
