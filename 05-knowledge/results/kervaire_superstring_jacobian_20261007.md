# Word compression, faithful extensions, and finite planar inverse receipts

**Status: CITED PREMISES + PROVED elementary transfer lemmas + FINITE-EXACT controls.**
The three supplied OpenAI theorems are accepted as premises, not re-audited.
The polynomial inverse-degree bound is a separate cited classical input.
No result here settles JC(2), supplies a new Keller counterexample, or turns
an abstract group solution into a polynomial map.

## 1. Recovered objects and the connection contract

The current two-variable program is
[the audited planar board](planar_jc_long_20260906_board.md), especially its
[full finite-row response](planar_jc_long_20260906_memory.md),
[completed Hamiltonian repairs](planar_jc_long_20260906_hamiltonian.md), and
[non-rational specialization theorem](planar_jc_long_20260906_nonrational.md).
It distinguishes finite-row compatibility, formal completion, and polynomial
termination. The new
[THM-4562 / slice lemma](../../01-canon/theorems/THM-4562-slice-lemma-and-the-jacobian-counterexample-is-fully-non-slice.md)
reduces dimension only under a locally nilpotent dual-field hypothesis.
The three-variable example is not a planar counterexample.

Closest proved mechanism: exact ordered composition with an inverse witness.
Canonical hostile: an operation can preserve every requested finite jet while
remaining globally nontrivial. Corrected near miss: a frozen-prefix response
rank is not the rank when earlier source coefficients may move. Least-used
sidecar: the full relative correction fibre, followed by a finite support
certificate.

The five live concepts are ordered words, actual inverse maps, relative
lifting, finite support, and group versus partially defined operation.
The connection is:

| Source | Target and preserved predicate | Missing information |
|---|---|---|
| Common superstring of literal operation words | Shared storage plus occurrence intervals; exact ordered map and inverse | It creates neither a factorization of a new map nor a legal domain for a partial operation |
| Unimodular group equation | Faithful abstract extension of an existing group | Realization in the original group of polynomial maps |
| Associated-graded detector | A specified filtered lifting obligation | Actual representative, extension class, and order |
| Bounded formal inverse jet | Exact polynomial inverse, once its degree support is certified | Universal vanishing of the obstruction coefficients for Keller inputs |

These are proof-obligation transfers, not equivalences between the research
problems.

## 2. What the supplied papers actually provide

The [shortest-common-superstring paper, Theorem 1.1](https://github.com/openai/math/blob/main/preprints/A-Polynomial-Time-2-Approximation-for-Shortest-Common-Superstring-September-24-2026/paper.pdf)
provides a deterministic polynomial-time algorithm with output length at most
twice optimum for finite **explicitly represented** strings. Symbol encoding
is part of the input size. The guarantee concerns its constructed algorithm,
not arbitrary greedy overlap merging. Its connection records retain actual
shared words and their alignments. These are exactly the metadata that a
proof-word storage interface must retain.

The [group Kervaire paper, Theorem 1.1](https://github.com/openai/math/blob/main/preprints/The-Kervaire-Theorem-for-Groups-September-24-2026/The-Kervaire-Theorem-for-Groups-September-24-2026.pdf)
says that if a word in \(G*\langle t\rangle\) has \(t\)-exponent sum
\(\pm1\), the coefficient map into the quotient by that relation is
injective. This supplies a solution in an overgroup. It does not assert
that the solution lies in \(G\). Its full nontriviality corollary also covers
other exponent sums, without thereby asserting coefficient injectivity in
each such case.

The [prime-three Kervaire paper, Theorem 1.1](https://github.com/openai/math/blob/main/preprints/The-Kervaire-Invariant-Problem-at-the-Prime-Three-September-24-2026/paper.pdf)
classifies the standard classes \(b_j\) as surviving precisely for
\(j=0,2,3\), with stems \(4\cdot3^{j+1}-2=10,106,322\).
Each surviving detection coset contains an actual order-three element.
This last assertion is additional to associated-graded survival. The paper
retains a correction and Moore-boundary construction to obtain it; it does
not assert that every representative has order three.

The integer \(322\) here is a stable homotopy stem. The older
[223/233 comparison](triplet_crt_223_233_20261004.md) concerns finite-field
orders, Collatz carries and source guards. No map preserving a predicate
between those objects has been supplied. Reordering their decimal digits
does not supply one.

## 3. Exact storage and gluing of operation words

Let every alphabet token denote an automorphism of one specified object,
with a registered inverse. Read words from left to right as operations;
thus \(F_{uv}=F_v\circ F_u\).

**Storage lemma.** Suppose the explicit words \(w_i\) occur in a tape \(T\)
at intervals \([\ell_i,r_i)\). Keep the tape, the intervals, and the token
dictionary. This data determines each \(F_{w_i}\) and its inverse exactly.
If \(P_j=F_{T[0:j)}\), then

\[
F_{w_i}=P_{r_i}\circ P_{\ell_i}^{-1}.                 \tag{1}
\]

The identity follows by cancelling the common prefix. It is valid for any
automorphism group and does not require that the token representation be
faithful as a free group. Applying the supplied superstring theorem therefore
gives factor-two token storage relative to the best single tape for the
same explicit words. Occurrence addresses and token payloads are additional
storage. Expanded polynomial degree, evaluation bit complexity, and
compression of already succinct input programs are different objectives.

For actual charts \(F_i\), the transitions
\(g_{ij}=F_j\circ F_i^{-1}\) satisfy
\(g_{jk}\circ g_{ij}=g_{ik}\). Their triangle products are identity by
literal cancellation. This is a constructive gluing certificate: the chart
maps and their inverse witnesses, rather than only pairwise similarity
scores, imply the cocycle.

The finite demonstration uses
\[
A(x,y)=(x+y^2,y),\qquad B(x,y)=(x,y+x^2),
\]
with lowercase tokens denoting inverses. The four words
ABa, BaB, aBb, BbA occupy intervals
\([0,3),[1,4),[2,5),[3,6)\) in ABaBbA.
Twelve input tokens become six tape tokens. A bounded exact subset algorithm
finds this tape; an independent permutation enumeration verifies its minimum
length. This implementation does **not** implement or claim the imported
polynomial-time approximation algorithm.

For partial maps, (1) needs actual domains and image membership as well.
In particular a Collatz common-future receipt is not an invertible group
generator. The
[clock-holonomy package](collatz_clock_holonomy_20261007.md)
retains its separate two-sided iteration clocks.

## 4. Finite jets do not supply the gluing certificate

For every \(m\ge2\), use the genuine polynomial automorphisms
\[
A_m(x,y)=(x+y^m,y),\qquad B_m(x,y)=(x,y+x^m).
\]
Their left-to-right commutator word is \(A_mB_mA_m^{-1}B_m^{-1}\).
Its first nonidentity homogeneous term is
\[
\bigl(-m x^m y^{m-1},\;m x^{m-1}y^m\bigr),           \tag{2}
\]
of degree \(2m-1\). One obtains (2) by expanding each shear through that
degree; terms involving two substitutions of the higher correction have
larger degree. The coefficient is nonzero in characteristic zero.
Nevertheless the exact word sends \((1,0)\) to \((0,1)\), for every \(m\).

Thus any prescribed finite jet horizon misses some globally nontrivial
holonomy, even when all local maps are polynomial, determinant one, and
equipped with explicit inverses. The repair is to retain the ordered
transition words or to verify their full polynomial composition. Increasing
a horizon without retaining a support bound does not yield a finite
universal gluing test.

## 5. An abstract extension may leave the desired operation category

Here is a six-element hostile to an overstrong reading of the group theorem.
In \(S_3\), let \(a=(23)\), \(c=(12)\), and
\[
w(t)=t a t^{-1} a t c.
\]
Its \(t\)-exponent sum is one. The supplied theorem gives a faithful
overgroup with a solution. However
\[
\{t a t^{-1} a t:t\in S_3\}=\{1,a\},
\]
so no solution lies in \(S_3\). The script checks all six substitutions.

For an actual polynomial-operation group, applying that theorem would
therefore leave a realization obligation: construct the new element as a
polynomial automorphism on the same plane, with its inverse and coefficient
action. Abstract coefficient injectivity alone does not discharge it.

Likewise, \(\mathbb Z/9\) with subgroup \(3\mathbb Z/9\) and
\((\mathbb Z/3)^2\) with its first factor have the same two associated-graded
factors. In the former, no lift of \(1\) in the quotient has order three;
in the latter every such lift does. This illustrates precisely why the
prime-three paper's representative assertion is additional evidence.

The elementary relative lifting test is useful in the JC row compiler.
Given linear maps \(\pi:V\to Q\), \(D:V\to W\), and a chosen representative
\(v\) with \(\pi(v)=q\), a correction preserving \(q\) and killing \(Dv\)
exists exactly when
\[
 -Dv\in D(\ker\pi).
\]
The full solution fibre, if nonempty, is a translate of
\(\ker\pi\cap\ker D\). This follows by writing the correction as
\(h\in\ker\pi\) and solving \(Dh=-Dv\). It is a repackaging of ordinary
linear algebra. The recovered JC work's complete raw response fibres
already implement this principle; compressed quotient data alone do not.

## 6. A finite polynomial-inverse receipt

Let \(F\in k[x,y]^2\) satisfy \(F(0)=0\) and
\(\det JF(0)\ne0\), over a characteristic-zero field. Let \(d=\deg F\)
and let \(G\) be the unique formal inverse with \(G(0)=0\).

**Finite receipt.** For an integer \(D\ge1\), if the jet of \(G\) through
total degree \(dD\) contains no terms of degrees \(D+1,\ldots,dD\), then
\(H=\operatorname{jet}_D G\) is the actual polynomial inverse.

Indeed \(G-H\) vanishes through degree \(dD\), so \(F(G)-F(H)\) does too.
The polynomial \(F(H)-\mathrm{id}\) has degree at most \(dD\); its whole
possible support therefore vanishes. Hence \(F\circ H=\mathrm{id}\).
Uniqueness of the formal inverse then gives \(H\circ F=\mathrm{id}\)
as well. Conversely any inverse of degree at most \(D\) passes this test.

The classical inverse-degree bound
\(\deg(F^{-1})\le(\deg F)^{n-1}\), when \(F\) is already known to be a
polynomial automorphism, is recorded in the primary
[Cheng--Wang--Yu paper](https://doi.org/10.1090/S0002-9939-1994-1195715-1).
Consequently in dimension two:
\[
F\text{ is a polynomial automorphism}
\quad\Longleftrightarrow\quad
[G]_{d+1},[G]_{d+2},\ldots,[G]_{d^2}=0.              \tag{3}
\]
The right side is a finite, exact obstruction block for a **given**
polynomial map. This standard degree-bound consequence is not a new
Jacobian theorem.

The executable normal form has integral coefficients and linear part
identity. It constructs the formal inverse through the required precision
by repeated correction \(G\leftarrow G+\mathrm{id}-F(G)\); each iteration
raises the first possible error degree. It then verifies both full
polynomial compositions. It correctly rejects too-small degree bounds.

Two scope controls distinguish this from a proof of JC(2).
The polynomial map \((x-x^2,y)\) has a unique formal inverse but an infinite
Catalan series; its Jacobian is not constant. The formal map
\[
\left(\frac{x}{1-x},\,y(1-x)^2\right)
\]
has determinant one and an inverse in formal power series, but is not a
polynomial input. Its compatible inverse jets satisfy both truncated inverse
identities at every finite order. They do not satisfy the vanishing-block
test (3): arbitrarily high pure-x coefficients remain nonzero. The missing
hypothesis is finite polynomial support. Neither example is a Keller
counterexample.

**Concrete next target.** For a fixed degree \(d\), retain the complete
polynomial input coefficients, impose the full constant-Jacobian equations,
and prove that their solution set annihilates the inverse block in (3).
For an inherited restricted source family, first show its parameters really
produce complete polynomial inputs of that degree; finite response jets
alone are insufficient. Ordered repair words can share their storage and
verification, while the obstruction block retains the still-unpaid
polynomial termination obligation. Universal vanishing for arbitrary Keller
inputs and arbitrary degree remains JC(2).

## 7. The full quadratic family: an exact symbolic discharge

This is a verification and elementary proof of the known quadratic case,
not a new Jacobian result. Let
\[
 F=(x,y)+H,\qquad
 H=(a x^2+bxy+c y^2,\;d x^2+exy+f y^2)
\]
over a characteristic-zero field. The complete Keller ideal in the six
coefficient variables is generated by
\[
 I=(2a+e,\ b+2f,\ 2ae-2bd,\ 4af-4cd,\ 2bf-2ce).       \tag{4}
\]
These are exactly the coefficients of \(\det JF-1\), not sampled
evaluations. An equivalent generating set is
\[
 g_1=2a+e,\quad g_2=b+2f,\quad
 g_3=2cd+ef,\quad g_4=ce+2f^2,\quad g_5=4df-e^2.      \tag{5}
\]
Writing the last three generators in (4) as \(i_3,i_4,i_5\), the explicit
identities are
\[
 i_3=e g_1-2d g_2+g_5,\quad
 i_4=2f g_1-2g_3,\quad i_5=2f g_2-2g_4.
\]
Solving these three identities for \(g_3,g_4,g_5\) verifies the reverse
ideal containment. The monic lexicographic Gröbner basis in the order
\(a,b,c,d,e,f\) is
\((g_1/2,g_2,g_3/2,g_4,g_5/4)\).

Let \(J\) denote the Jacobian matrix of \(H\). The unrestricted formal
inverse through degree four is
\[
 G_1=(x,y),\quad G_2=-H,\quad G_3=J\,H,\quad
 G_4=-J G_3-H(H).                                   \tag{6}
\]
The script reduces all eight coefficients of \(G_3\) and all ten
coefficients of \(G_4\) by the basis (5). All 18 remainders are zero.
For every coefficient it also verifies the exact identity between that
coefficient and its returned combination of the five basis elements.
A separate sparse-polynomial implementation verifies both truncated
compositions for the unrestricted six-parameter expressions in (6).
This is universal symbolic ideal membership, not a finite parameter census.

There is an elementary proof independent of Gröbner reduction. Equation (4)
says \(\operatorname{tr}J=\det J=0\), so the two-by-two Cayley--Hamilton
identity gives \(J^2=0\). Euler's homogeneous identity gives
\(J(x,y)^{\mathsf T}=2H\), whence \(JH=0\).
For the derivation \(D=H_1\partial_x+H_2\partial_y\), this is \(D(H)=0\).
Thus \(D^2x=D^2y=0\). Leibniz's rule shows that a monomial of degree \(s\)
is killed by \(D^{s+1}\), so \(D\) is locally nilpotent. Its exponentials
are inverse ring automorphisms, and on the coordinate generators they are
\(\exp(\pm D)(x,y)=(x,y)\pm H\). Consequently the polynomial inverse is
exactly \(\mathrm{id}-H\), which also proves the vanishing block in (6).

## 8. Fixed-degree JC as explicit finite algebraic obligations

Fix \(d\ge2\), use all coefficients \(\mathbf c\) of a normalized polynomial
map \(F=\mathrm{id}+\) terms of degrees \(2,\ldots,d\), and let
\(I_d\subset\mathbb Q[\mathbf c]\) be generated by all coefficients of
\(\det JF-1\). Normalization keeps constant term zero and linear part
identity; any complex Keller map can be put in this form by affine target
changes without changing whether it has a polynomial inverse.

Formal inversion produces each coefficient of \(G\) polynomially in
\(\mathbf c\). Let \(p_1,\ldots,p_s\) be the coefficients of its homogeneous
parts of degrees \(d+1,\ldots,d^2\). By (3), the degree-at-most-\(d\)
complex Jacobian statement is equivalent to all \(p_j\) vanishing on
the complex zero set of \(I_d\). The
[Hilbert Nullstellensatz](https://stacks.math.columbia.edu/tag/00FS)
therefore makes it equivalent to the finite list
\[
 p_j\in\sqrt{I_d},\qquad j=1,\ldots,s.                \tag{7}
\]
The radical may be computed over \(\mathbb Q\): membership is unchanged
after the faithfully flat extension to \(\mathbb C\).

Each obligation has a finite polynomial certificate via an added variable:
\[
 1\in (I_d,\ 1-zp_j)\subset\mathbb Q[\mathbf c,z].     \tag{8}
\]
For completeness, if \(p_j^N\in I_d\), the identity
\[
 1=(1-zp_j)\sum_{i=0}^{N-1}(zp_j)^i+z^N p_j^N
\]
proves (8). Conversely substitution \(z=1/p_j\) in the localized quotient
shows that (8) forces a power of \(p_j\) into \(I_d\). A finite expression
for (8), or an explicit identity for \(p_j^N\), is checked by exact
polynomial arithmetic. Gröbner algorithms decide each fixed finite
instance; their potentially very large cost is not bounded by the
superstring theorem.

The quadratic calculation above proves the stronger ordinary membership
\(p_j\in I_2\). In a higher-degree case a nonzero ordinary-ideal remainder
would not by itself refute the radical condition; powers or (8) must be
tested. No higher-degree radical computation or universal bound is claimed.
This gives the precise next compiler target: generate the full
coefficient ideal and inverse block for a declared degree or an honestly
scoped polynomial subfamily, then retain exact radical-membership receipts.
It does not replace the open assertion that these receipts exist for
every degree.

## 9. Reproduction and finite universe

[Script](../../04-computation/experiments/kervaire_superstring_jacobian_20261007.py)
and [saved output](kervaire_superstring_jacobian_20261007.out).

~~~powershell
python -B 04-computation/experiments/kervaire_superstring_jacobian_20261007.py
python -B -O 04-computation/experiments/kervaire_superstring_jacobian_20261007.py
~~~

The universe is four explicit three-token words, all 24 orderings for the
independent storage check, all 64 chart triples, nine integer test points,
shear degrees 2 through 9, the two positive inverse examples with their
declared degree tests, eight Catalan inverse bounds, formal jets of degrees
3 through 12, all six \(S_3\) substitutions, and six malformed API inputs.
The symbolic extension additionally verifies both directions of the five
ideal generators, the full quadratic Jacobian coefficient identity, all
18 inverse-block coefficient identities and their zero remainders, and
both unrestricted truncated compositions. It uses SymPy 1.14.0 over exact
rational coefficients. The 279 checks pass in both modes with identical
output. The universal lemmas have proofs above; the finite checks are
independent controls rather than extrapolation.
