# Creation, parity carries, and a Fano-to-E8 Collatz model

**Status: PROVED scoped constructions and decoder obstructions; FINITE-EXACT
integer certificates; CITED background; OPEN Collatz termination.**
Date: 2026-09-25. These results do not assert a solution of Collatz or novelty
for the classical parity conjugacy, binary multiplication, or elliptic descent.

## 1. Inheritance and the question that survives the numbers

The closest proved mechanisms are the [finite-word carry machine](creative_transducer_20260925.md),
the [fruit elliptic correspondence](catalan_elliptic_20260921_elliptic.md),
and the [Fano/Hamming Construction-A lattice](arithmetic_braids2_20260917_fano_code.md).
The canonical hostile is the arbitrarily long all-ones binary prefix, whose
plus orbit initially expands. The corrected near miss is the older elliptic
ternary tree: its group points are valid, but positivity of the three fruit
coordinates is not preserved. The least-used useful coordinates are the
elliptic square class and the residual carry function of a binary prefix.

The [signed Pell and thirty-six construction](thirtysix_bridge_20260925.md)
remains a separate exact model. The current numbers are the fruit integers,
not decimal strings of threes ending in one. Their ratios therefore require
a different explanation from that earlier decimal affine clock.

Anchor: a creation decoder whose certificates can be checked intrinsically.
Niche: the actual first nonlinear Collatz carry acting on the Fano plane and
its lattice. Wildcard: what eightfold Clifford periodicity preserves, and
what an iterated parity encoder must keep beyond a finite color palette.

| Live concept | Retained object | Decisive boundary |
|---|---|---|
| Fruit ratios | exact projective point and literal input | approximate 4.2/8.4 do not determine binary carries |
| Elliptic halving | square class, torsion phase, signed point | actual Collatz tripling raises canonical height |
| Fano plane | input labels and parity labels | plus parity is nonlinear, minus is linear at depth three |
| E8 and Bott | lattice frame, Clifford generators, divisibility guard | an isometry or stable module is not arithmetic descent |
| Three carry colors | finite controller and unbounded tape | scan termination differs from iteration termination |
| Creation certificate | reconstruction digits and well-founded rank | the rank must decrease along the actual target operation |

## 2. The input correction and its genuine arithmetic signal

The [complete number audit](creation_numbers_20260925.md) preserves the three
comma strings exactly. Removing commas gives `(a,c,d)`, of lengths
`(81,79,81)`, where `d=10b+9`. Replacing only `d` by `(d-9)/10` recovers
the previously verified fruit triple `(a,b,c)` with lengths `(81,80,79)`.
Its exact equation and ratios are

\[
\frac a{b+c}+\frac b{a+c}+\frac c{a+b}=4,\qquad
\frac ab=4.189186440639\ldots,\quad\frac bc=8.431275128733\ldots.
\]

The literal triple has fruit sum about `2.74374`, so the repair cannot be
omitted. The 2014 source and 2025 corrigendum are linked in the audit.

For any positive solution of fruit sum four with `a` largest, the exact
dominance bound is

\[
2+\sqrt3<\frac a{b+c}\le\frac{7+\sqrt{65}}4.
\]

The upper equality occurs for positive real triples with `b=c`, and never
for rational triples. This constrains the largest coordinate relative to
the other two together; it does not fix their split or impose digit lengths.

All four distinct supplied/repaired integers have independently replayed
ordinary plus routes to four. Their ordinary lengths are `2025,2116,1831,1778`
for `(a,b,c,d)`. On the minus sheet the repaired triple enters the three
basins with minima `(17,5,1)`. That striking assignment is a finite fact:
the next positive elliptic control `13G+T` has basin assignment `(5,5,17)`.
The fruit identity alone cannot imply one coordinate per minus basin.

## 3. A creation decoder that really terminates

The [elliptic construction](creation_elliptic_20260925.md), statements E1--E3,
works on the explicitly generated subgroup

\[
E:y^2=x^3+109x^2+224x,\quad
L=\langle G=(-4,28),T=(56,728)\rangle\cong\mathbb Z\oplus\mathbb Z/6.
\]

It reads both coefficient parities directly from the point, using the
rational square class

\[
\alpha(mG+kT)=(-1)^m14^k\quad\text{in }\mathbb Q^*/\mathbb Q^{*2}.
\]

If these parities are `r,t`, subtract `rG+tT`, halve by explicit rational
quadratic tests, and select the unique half with even torsion phase. Thus

\[
P=2P'+rG+tT.
\tag{C1}
\]

This is a coordinate algorithm; the free coefficient is not supplied to it.
The proof identifies `m'=floor(m/2)` and supplies the nonnegative integer rank

\[
D(m,k)=2\rho(m)+(k\bmod2),\qquad
\rho(m)=\begin{cases}m&m\ge0,\\-m-1&m<0.\end{cases}
\]

It strictly decreases outside a specified six-point terminal bank.
The large fruit point has the short creation certificate

\[
9G\longrightarrow4G\longrightarrow2G\longrightarrow G\longrightarrow O.
\]

Reading (C1) backwards creates the point from its terminal state and bits.
The construction is globally proved on `L`, independently of a full rational
generator assertion. It is a concrete model of the desired kind of proof:
intrinsic decoding, retained phase, exact reconstruction, and a rank.

The actual shortcut Collatz lift sends an odd `mG` to `(3m+1)G/2` with the
appropriate rational half. It sends `9G` to `14G`, not `4G`. No decreasing
rank for that replacement has been obtained. Raw doubling also cannot
produce a positive fruit point: every finite real double has nonnegative
abscissa, whereas fruit positivity requires a negative one. The signed
extension and parity correction are essential to the successful decoder.

Its division quotients give a second exact connection to the user's counts:

\[
|L/bL|=b\gcd(b,6),\qquad
|L/2L|=4,\quad|L/3L|=9,\quad|L/6L|=36,\quad|L/8L|=16.
\]

These are classes, not preimage counts. Three successive binary refinements
give `4*2*2=16` after selecting the kernel phase, rather than eight or 64.

## 4. The first Collatz carry is an exact Fano-line trade

Let

\[
T_\sigma(n)=\begin{cases}n/2&n\text{ even},\\(3n+\sigma)/2&n\text{ odd},\end{cases}
\quad
q_{\sigma,d}(n)=\sum_{j=0}^{d-1}(T_\sigma^j(n)\bmod2)2^j.
\tag{C2}
\]

The input of `q` is the residue modulo `2^d`; the output is the ordered
parity itinerary. These are different coordinates on the same finite set.
The compatible maps are the inverse of the classical 2-adic conjugacy
studied by [Bernstein--Lagarias (1996)](https://websites.umich.edu/~lagarias/doc/bernstein.pdf).
No new conjugacy claim is made here.

For `d=3`, write `n=x0+2x1+4x2`. In characteristic two,

\[
q_{+,3}(x_0,x_1,x_2)=(x_0,x_1,x_2+x_0+x_0x_1),
\tag{C3}
\]
\[
q_{-,3}(x_0,x_1,x_2)=(x_0,x_1+x_0,x_2+x_0+x_1).
\tag{C4}
\]

Thus plus is the transposition `(1 5)` on the eight residue labels, while
minus is linear. The exact nonlinear defect of plus is

\[
q(x)\mathbin\oplus q(y)\mathbin\oplus q(x\mathbin\oplus y)
=4(x_0y_1+x_1y_0).
\tag{C5}
\]

Use the usual Fano lines `{x,y,x XOR y}` on the seven nonzero binary
vectors. Applying `q` preserves the three lines through four,

\[
145,\quad246,\quad347,
\]

and exchanges the other four as follows:

\[
\{123,167,257,356\}\quad\longleftrightarrow\quad
\{127,136,235,567\}.
\tag{C6}
\]

This is an actual relation to the Fano plane, including its failure of
linearity. The repair is the transported addition

\[
x\star y=x\mathbin\oplus y\mathbin\oplus4(x_0y_1+x_1y_0).
\tag{C7}
\]

Then `q(x star y)=q(x) XOR q(y)`. This is still a commutative elementary
abelian group: the correction is a coboundary of `x0*x1`. It is not a
noncommutative Heisenberg group or evidence of a new multiplication law
on the ordinary integers. The input frame and the carry in (C5) are the
information that a bare seven-point incidence picture loses.

## 5. The carry changes the E8 frame through D8

The [lattice companion](creation_fano_20260925.md) makes (C3) geometric.
Start with the affine binary code `C=RM(1,3)` evaluated on the eight cube
points. Its Construction-A lattice is E8. Permuting the coordinate slots
by `(1 5)` gives a neighboring copy E8'. The two codes intersect in a
dimension-three subcode. The lattices intersect in an index-two D8 in each:

\[
E_8\supset D_8\subset E'_8.
\tag{C8}
\]

There are 112 common roots and 128 other roots in each E8 copy. Pairing
coordinates `(0,4),(1,5),(2,6),(3,7)` and applying a Hadamard change of basis
turns the swap `(1 5)` into a sign flip of one coordinate. This exchanges
the two spinor glue classes over D8. The minus map (C4), being linear,
preserves the original affine code instead.

The same companion constructs seven octonionic left-unit operators and
an eighth block operator on two copies. They satisfy the exact negative
Clifford relations and preserve `E8 direct-sum E8`. This is the algebraic
eightfold Clifford structure associated with real Bott periodicity;
its primary references and conventions are retained there. Stable
periodicity is not a claim that an arithmetic orbit returns every eight
steps.

For a chosen E8 root `rho`, encode the integer by `(n*rho,rho)`.
A guarded rational affine block then realizes each actual Collatz step
exactly. This supplies a faithful lattice model, including divisibility.
It also exposes why the norm alone fails: it is proportional to `n^2+1`.
The Clifford operators are isometries, but the combinations implementing
tripling and halving need not be.

## 6. A proof that the carry complexity does not reset after eight levels

**C9 (PROVED).** For every `d>=3`, the highest output bit of `q_{+,d}`
has Boolean algebraic degree exactly `d-1` in the ordinary input bits.
For `q_{-,d}` its degree is at most `d-2`.

Here algebraic degree means the degree of the unique multilinear
polynomial over F2. The theorem concerns these fixed coordinates, not
the complexity of every possible encoding or arithmetic circuit.

First, the map is triangular:

\[
(q_d(n))_{d-1}=x_{d-1}+f_d(x_0,\ldots,x_{d-2}).
\tag{C9a}
\]

Indeed two inputs differing by `2^(d-1)` have the same first `d-1`
parities; after those steps their difference is an odd multiple of one.
Their last parity therefore flips. Induction proves bijectivity and
compatibility at all depths, without a convergence assumption.

For the plus map, split the permutation into its even and odd fibers:

\[
q_d(2m)=2q_{d-1}(m),\qquad
q_d(2m+1)=1+2q_{d-1}(3m+2).
\]

The product of the two signs of `q_(d-1)` is positive, leaving the sign
of `m -> 3m+2` modulo `2^(d-1)`. On `2^M` points with `M>=2`, translation
by two has even sign. Multiplication by three has odd sign: pair `r`
with `r+H`, where `H=2^(M-1)`. The induced permutation of pairs contributes
its sign twice, and the number of pair flips has the parity of

\[
\#\{r:H\le3r<2H\}=\lceil2H/3\rceil-\lceil H/3\rceil,
\]

which is odd for every power of two `H>=2`. Hence `sign(q_(+,d))=-1`.

In the triangular form (C9a), that sign is `(-1)^(sum f_d)` over all
lower-bit assignments. This sum modulo two is precisely the coefficient
of the full monomial `x0*x1*...*x_(d-2)` in the algebraic normal form.
It is one. The upper bound and the required nonzero highest-degree
coefficient prove the assertion.

Finally `q_-(n)=q_+(-n)` at each finite depth. Negation modulo `2^d`
has `2^(d-1)-1` transpositions and thus odd sign for `d>=2`.
For `d>=3` the minus permutation has even sign, eliminating the full
lower-bit monomial and proving the stated bound. The boundary `d=2`
must be excluded. Finite controls through depth12 give minus degree
exactly `d-2` for `d>=3`; equality at all depths is not claimed here.

The degree sequences through depth12 are

```text
plus:  1,1,2,3,4,5,6,7,8,9,10,11
minus: 1,1,1,2,3,4,5,6,7,8, 9,10.
```

The proof was independently audited by a second agent. It does not
conflict with a small carry machine: bounded local state and increasing
degree of iterated output formulas are different notions.

## 7. A second obstruction: the itinerary encoder needs infinitely many residual states

**C10 (PROVED).** No deterministic finite-state, synchronous, least-significant-
bit-first Mealy transducer computes `q_+` correctly at every depth.

Let `q_+` be the compatible map on the 2-adic integers. The bijections above
make it injective. Input and output streams are infinite, with ordinary
nonnegative integers represented by zero-padded binary streams. After
reading the input prefix of `k` ones, write

\[
n=2^k m+2^k-1.
\]

Its first `k` parities are all one, and the exact shortcut identity gives

\[
T_+^k(n)=3^k m+3^k-1.
\]

Consequently the residual output map on the unread input tail is

\[
R_k(m)=q_+(3^k m+3^k-1).
\tag{C10a}
\]

These residual maps are pairwise distinct: evaluate at `m=0` and use
injectivity on the distinct integers `3^k-1`. A finite-state synchronous
machine has only finitely many residual maps, one per reachable state.
This contradiction proves the claim. The argument was independently
audited; no orbit is assumed to converge.

The inherited three-state multiplier is entirely compatible with this
result. It implements **one** Collatz step on a finite tape. The full
itinerary encoder composes unboundedly many steps and is a different
map. Finite control with an unbounded tape, stack, or integer register is
not excluded by C10. This supplies a precise reason to keep such a
register when extending the Fano diagram.

## 8. What object to pursue next, and what would constitute success

An appropriate candidate is a recursive arithmetic object with three
separate parts: a local parity/lattice frame; the finite binary tape with
its end marker; and an exact affine carry record for a variable-length
block. The existing three-state multiplier supplies the local rewrite
rules. The lattice frame records the Fano trade; the carry record retains
the scale and translation lost by its finite incidence quotient.

The following comparisons specify the maps rather than inferring
equivalences from similar counts:

| Source -> target | Map and preserved predicate | Lost data / required record | Cheapest decisive test |
|---|---|---|---|
| Integer -> E8 pair | `n -> (n*rho,rho)`, exact guarded affine step | finite frame alone loses magnitude | test actual `3n+1` integrality and norm on odd n |
| Binary cube -> Fano frame | `q_(+,3)`, exact parity labels | XOR addition changes by (C5) | `q(1) XOR q(2) != q(3)` |
| E8 -> neighboring E8 | coordinate swap `(1 5)`, norm | code frame changes | D8 intersection,112 shared roots |
| Elliptic point -> creation digits | square-class correction plus selected half | positive fruit locus is not invariant | positive `9G` decodes to nonpositive-fruit `4G` |
| Finite tape -> full parity itinerary | iterated carry scans | residual functions are unbounded | sections (C10a) at all-ones prefixes |

Three concrete research proposals follow from the successful constructions
and their hostile controls:

1. **Adaptive carry blocks.** Use the already proved maximal-run records in
   [the transducer note](creative_transducer_20260925.md), then search for a
   well-founded rank of the remaining tape plus affine carry. The local
   tests must cover every legal block parameter, including arbitrarily
   long ones prefixes. This differs from bounding all scan lengths by a
   fixed number. A fixed-horizon, finitely residue-indexed coercive
   polynomial potential is ruled out in the lattice companion.
2. **Transported frames.** Let the code/lattice frame change under each
   arithmetic carry rather than demanding that every step preserve one
   Fano labeling. The D8 intersection gives an explicit first transition.
   The new obligation is a rank across changing frames and integer
   registers; counting 112 shared roots is not such a rank. An atlas at
   deeper levels must also respect C9 and C10.
3. **Transfer the creation certificate, not just the object.** The elliptic
   construction demonstrates all four required ingredients: readable
   parity, exact inverse creation, a terminal bank, and strict rank.
   A Collatz transfer must preserve the actual arrow (odd `9 ->14`, for
   example) and supply a corresponding rank. Reusing the binary-creation
   rank while changing its operation fails at that first example.

The Fibonacci three-color connection is correspondingly precise. The
inherited [Zeckendorf charge note](duck_zeckendorf_20260925.md), Z7--Z12,
repairs overlapping representations by a carry coboundary. Formula (C5)
is another explicit coboundary repair, on a different object. In both
cases a coloring survives only with the required carry or transported
operation. This does not identify Fibonacci addition with Collatz time.

The combined positive result is a terminating intrinsic creation decoder
and a faithful Fano/Clifford model of actual arithmetic, together with
proofs of the information each finite picture loses. The outstanding
obligation is a universal well-founded certificate for repeated actual
Collatz blocks. It remains **OPEN**; none of the eightfold, lattice, or
elliptic identifications is substituted for that obligation.

## 9. Exact reproduction

```text
python 04-computation/experiments/creation_decoder_20260925.py
python -O 04-computation/experiments/creation_decoder_20260925.py
```

The [script](../../04-computation/experiments/creation_decoder_20260925.py)
checks every residue for both signs at depths1--12: 16,380 inversions by
an independent affine-word formula, triangular degree coefficients,
permutation signs, and signed conjugacy. It also checks every pair and
triple for the transported eight-element group, the full seven-line
trade, the affine-sign proof through12 bits, and768 residual-section
identities. It writes [the exact output](creation_decoder_20260925.out).
All checks use exceptions and remain active under optimization.

The three companion scripts separately reproduce the fourteen integer
trajectory certificates, the294-point elliptic decoder controls, and the
Fano/code/root/Clifford arithmetic. Their notes state their complete finite
universes and independent verification paths. All-depth statements in
this note rest on their proofs, not the finite cutoffs.
