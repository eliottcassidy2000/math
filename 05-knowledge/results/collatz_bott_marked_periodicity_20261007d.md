# Marked periodicity: stable matrices and the eight-bit child's exact parameter

**PROVED:** the elementary marked matrix, endpoint isometry, finite-precision
monodromy, Newton-polynomial and continuation pullback statements below.
**FINITE-EXACT:** the declared reproducible controls and the application to
the separately authenticated new child head. **CITED:** Bott and Clifford
periodicity in their own categories. **OPEN:** complete coverage of the
target phase and grounding every emitted child. No eight-step Collatz cycle
or new ROOT discovery is claimed.

## 1. Inheritance and what periodicity can preserve

The nearest exact bridge is [creation_fano_20260925](creation_fano_20260925.md),
F2–F3: explicit Clifford matrices carry a marked affine operation, and a
primitive lattice vector retains its integer divisibility guard. Its F5
already rules out specified fixed residue-polynomial descent potentials.
The [creation decoder](creation_decoder_20260925.md) distinguishes an
intrinsic elliptic halving decoder from the different actual Collatz map.

The relevant correction is **MISTAKE-484**, in
[the mistakes ledger](../../01-canon/MISTAKES.md): the theta series of Z8
and E8 are different. The surviving operation `L -> L direct_sum E8`
preserves the discriminant form and multiplies the theta series; equal rank
does not identify the two lattices. Likewise the score-lattice membership
in [THM-868, E8 bridge score lattice](../../01-canon/theorems/THM-868-e8-bridge-score-lattice.md)
is separate from its historical informal identification of several eights.
Neither that identification nor the historical `bott_periodic_backbone_s199`
count comparison is used here.

The anchor is the emitted child phase from
[eight-bit completion](collatz_eight_bit_completion_20261007c.md). The niche
is a lossless marked matrix corner; the wildcard is a finite-degree reader
of the exact arithmetic monodromy. The live board is source, guard, carry,
frame, binary precision and terminal obligation.

**Primary-source scope.** Bott's original theorem concerns stable homotopy
groups: unitary period2, and an orthogonal/symplectic four-degree exchange
giving period8. See [Bott, 1959, pp.314–315, (1.5)–(1.7)](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/bott4.pdf).
For the negative-square real Clifford convention, algebraic stabilization is
`Cl_(0,k+8) ≅ Cl_(0,k) tensor M16(R)`; graded-module quotients also have a
specified degree-eight isomorphism. See [Atiyah–Bott–Shapiro, Clifford
modules, §4 Table1 and §6 Proposition6.8](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/abs.pdf).
The dimension changes by a factor16. These are not return-time statements
for an arithmetic orbit. The proofs in Sections2–6 are independent of them.

| Source | Target/map | Preserved | Lost without sidecar | Cheapest hostile |
|---|---|---|---|---|
| Marked affine matrix M | M tensor I16 | Products and full entries | Chosen corner/frame if quotiented away | Same similarity class, different carries |
| Word matrix | Homotopy class | Membership in GL+(2,R) | Guard, carry, growth | `(1)` at3 versus `(2)` at9 |
| Exponent parameter t | Normalized endpoint j(t) | Exact binary precision, bijectively | Source and prefix if j stands alone | Same46 bits, different47-bit guard |
| Continuation word | Unique parameter phase | Native valuation cylinder | ROOT boundary/payment | Formal legality is not grounding |

## 2. Marked stabilization is lossless; a stable class is insufficient

For a positive valuation word w with carrier `(P,Q,B)`, define

\[
M_w=\begin{pmatrix}P/Q&B/Q\\0&1\end{pmatrix}.
\]

For any d>=1, `M -> M tensor I_d` is injective and multiplicative. Keeping
the displayed tensor frame, select coordinates `(1,1),(2,1)` to recover M.
More intrinsically, keep the marked corner `I_2 tensor e_1e_1^T` and its basis.
The public reader verifies every scalar block, not only the sampled corner,
before accepting the recovery. This is an elementary linear-algebra model
of retaining data through stabilization, not a claim that the number16
compresses an integer or proves a dynamical property.

By contrast every such matrix is homotopic to the identity through

\[
M(s)=\begin{pmatrix}1-s+sP/Q&sB/Q\\0&1\end{pmatrix},\quad0\le s\le1.
\]

The determinant stays positive. Any invariant constant on these homotopies
therefore assigns the same value to every word matrix in this carrier. It
cannot distinguish the actual expanding step `3 --(1)-->5` from the actual
decreasing step `9 --(2)-->7`. The intermediate homotopy matrices are not
asserted to preserve an integer native guard.

Even similarity is too coarse for source legality:

\[
M_{12}=\begin{pmatrix}9/8&5/8\\0&1\end{pmatrix},\quad
M_{21}=\begin{pmatrix}9/8&7/8\\0&1\end{pmatrix},\quad
M_{21}=C M_{12}C^{-1},\quad C=\begin{pmatrix}1&-2\\0&1\end{pmatrix}.
\]

Their native odd-source cells are respectively11 and9 modulo16, disjoint
at the same supplied source. Keeping a complete operator and a marked
equivalence can preserve this information; forgetting to its class does
not. This is not an objection to a categorical equivalence with its inverse
functor and all markings retained.

## 3. Exact compressed coordinate for the emitted child

Write

\[
F=924745897,\quad T=2^{32},\quad E(t)=F+Tt,\qquad t\in\mathbb Z_{\ge0}.
\]

The prior D8 rule emits `M_E=2^E-1`. Its `E-1` initial ones are followed by
the fixed actual head

```text
v=(2,2,4,1,3,3,1,3,3,1,1,1,2,2,3).
```

It has length15, cost32 and carrier
`P=14348907, Q=4294967296, B=4309583573`. If Y(t) is its endpoint after
the initial ones, then

\[
Y(t)=\frac{2P3^{E(t)-1}+B-P}{Q},\qquad Y(t)\equiv5\pmod8.
\]

The normalized coordinate is exactly

\[
\boxed{j(t)=\frac{Y(t)-5}{8}
 =\frac{3^{924745911+2^{32}t}-8589800907}{2^{34}}.}         \tag{1}
\]

This is a source-specific symbolic integer, not a new assumption about an
unknown orbit. It is derived from the authenticated fixed prefix. The
implementation reads it modulo2^h with numerator precision2^(h+34), using
modular exponentiation only. It never expands `2^E`, `3^E` or their orbits.

Let `u=3^T`. Elementary binary lifting gives `v2(u-1)=34`. For distinct
nonnegative t,s, subtraction of (1) yields

\[
\boxed{v_2(j(t)-j(s))=v_2(t-s).}                         \tag{2}
\]

Indeed `v2(u^(t-s)-1)=34+v2(t-s)` and the remaining multiplier is odd.
Consequently j extends uniquely to an isometric bijection of Z2 with Z2.
Surjectivity can be seen without any compactness assumption: modulo2^h,
(2) gives an injection of2^h residues into2^h residues, hence a bijection;
the unique inverses are compatible on reduction. This yields an exact
one-bit Hensel inverse for a requested endpoint residue.

## 4. Monodromy, and why the precision does not stabilize

With `C=8589800907`, (1) gives the affine parameter transition

\[
j(t+1)=u j(t)+C\frac{u-1}{2^{34}}.                      \tag{3}
\]

The additive coefficient is odd, and u is1 modulo4. More directly, (2)
conjugates this finite map to `t -> t+1`: modulo2^h it has one cycle of
length **exactly2^h**. Raising the precision one bit doubles the cycle
length. The native guard is therefore supported by a compatible family of
finite quotients, not by a fixed eight-state or eight-step periodicity.

Equation (3) is parameter advance, not one Collatz step. Its iterate can
be evaluated at arbitrary supplied step s using `u^s` and
`C(u^s-1)/2^34` with the same guarded modular division.

There is a sharper finite-degree reader. Put `u=1+2^34 w`, with w odd, and
`A=924745911`. For every h>=1,

\[
\boxed{j(t)\equiv j(0)+3^A\sum_{i=1}^{R}
 2^{34(i-1)}w^i\binom ti\pmod{2^h},\qquad R=\lceil h/34\rceil.} \tag{4}
\]

The binomial theorem proves (4); every omitted coefficient is divisible by
2^h. All binomial coefficients are integers, including zero for i>t.
This degree is minimal among integer-valued rational polynomials agreeing
for every nonnegative t modulo2^h. The R-th forward difference at zero of
j is

\[
\Delta^Rj(0)=3^A(u-1)^R/2^{34},\qquad
v_2(\Delta^Rj(0))=34(R-1)<h.
\]

A polynomial of degree<R has zero R-th difference, contradicting the
pointwise congruence. Thus the reader is linear through34 bits, quadratic
through68, cubic through102, and so on. At35 bits the exact hostile is
`j(2)-2j(1)+j(0)=2^34 mod2^35`. A linear reader cannot be extrapolated past
its precision. This is a controlled algebraic compression with an explicit
precision bill; it is not an approximation.

## 5. Pulling every continuation guard back to the same parameter

A nonempty positive valuation word w can start at a number5 modulo8 if
and only if its first exponent is at least3. For such a word, let its cost
be a and carrier `(p,q,b)`. Its native odd-source residue is

\[
y_w=(q-b)p^{-1}\pmod{2q},\qquad q=2^a.
\]

It lies in5 modulo8. Equations (1)–(2) therefore pull this cylinder back
to **one** exact class

\[
\boxed{t\equiv t_w\pmod{2^{a-2}}.}                    \tag{5}
\]

Every word beginning with an exponent>=3 has a class; words beginning1 or2
have none. Appending letters refines the parent's class, which the reader
checks by retaining the complete ordered word. The relative natural density
of (5) in the original parameter t is2^(2-a), as also holds for Haar
parameter measure. This is a parameter-frequency statement, not a source
convergence probability. In particular the next exponent is>=3 and
`P(next exponent=a)=2^(2-a)` for each a>=3 under that parameter measure.

The native iff includes formal continuation at ROOT. To export a strict
prefix, additionally impose

\[
E(t)-1\ge2(32+a).                                      \tag{6}
\]

For every prefix of the fixed head and w, its coefficient relative to
the immutable source M_E is at least
`(3/2)^(E-1) 2^(-(32+a)) >= (9/8)^(32+a)>1`.
The carries are nonnegative and initial ones grow. Thus every nonempty
prefix stays above the original source and cannot meet ROOT. The public
packet stores this sufficient cutoff separately from its exact formal
phase. It does not claim optimality, payment, or a ROOT suffix.

The cutoff also exposes a scoped obstruction: any fixed finite collection
of continuation words, with maximum cost a, is entirely expanding on the
members satisfying (6) wherever those words are native. More precision
can select which such word is legal but cannot turn its unchanged carrier
into a forward descent. Adaptive word costs and common-future reductions
remain outside that obstruction.

## 6. A real47-bit rule, and a missing-bit hostile

The independently owned
[eight-child routes](collatz_eight_child_routes_20261007d.md) authenticate a
39-letter source head of cost81 with a guarded deletion phase

\[
E\equiv924745897\pmod{2^{79}}.
\]

Its first15 letters are exactly v from Section3. The remaining24 letters
have cost49 and begin with3. Feeding that actual suffix into (5) gives

\[
\boxed{t\equiv0\pmod{2^{47}}.}                        \tag{7}
\]

This matches the independently derived exponent phase by ordinary CRT.
It is the same law expressed in the compressed endpoint coordinate; no
new paid coverage is added by counting it twice. The47-bit endpoint
reader in (4) is exactly quadratic. No astronomical source expansion is
needed to authenticate this interface.

The nearby parameters `t=0` and `t=2^46` have identical j residues through
46 bits by (2), but only the former satisfies (7). Thus that coarser
coordinate cannot decide this actual new rule. The missing bit must be
retained or acquired; matrix stabilization and a periodic quotient do not
manufacture it. More generally each precision level has such separating
parameter pairs. The rule's actual smaller-child proof obligation remains
owned by its marked receipt and supplied ROOT interface.

## 7. Reproduction and limits

Run both:

```text
python -B -X utf8 04-computation/experiments/collatz_bott_marked_periodicity_20261007d.py
python -B -O -X utf8 04-computation/experiments/collatz_bott_marked_periodicity_20261007d.py
```

The6,940 exact controls include every quotient permutation through10 bits,
independent replay of the fixed word on modest residue representatives,
unique bit lifting, long-step monodromy and Newton formulas through103 bits,
all155 words on alphabet1..5 of lengths1..3 (93 compatible), nested guard
refinements, the actual47-bit child rule, and exact-type/forged-packet hostiles.
Matrix controls test marked recovery for dimensions1,2,16 and the disjoint
12/21 guards. These finite computations test the interfaces; the all-depth
statements follow from the algebra and valuation proofs above.

No production method performs an orbit search or assumes an unsupplied
ROOT certificate. The new reusable operation is a precision-aware,
source-labelled guard pullback. Bott periodicity supplies a useful question
about what a stabilization retains, but supplies no arithmetic completion
premise in this argument.

Normal, optimized and saved stdout agree. LF-normalized SHA256:
`4c27ba6a5ed5ba69f07a43455c0a78268d57d69a2ae40b062e04742eeef6cbeb`.
