# Ordered affine-tree recurrences

Status: **PROVED** by the elementary arguments below; bounded independent
coefficient and native-integer checks are **FINITE-EXACT**. No external priority
claim. The clock is the sum of valuations, not an auxiliary tree depth.

## 1. Objects, chronological parsing and exact sources

For a finite nonempty valuation word `w=(a_1,...,a_J)`, put

\[
 D_w=\sum a_i,\qquad U_w(x)=(3^Jx+C_w)/2^{D_w}
       =\Lambda_w x+B_w.
\]

Each letter means `(3x+1)/2^a`. Composition follows the displayed word from
left to right. Thus `C_a=1`, `Lambda_w=3^J/2^D`, `B_w=C_w/2^D`, and

\[
 C_{uv}=3^{J_v}C_u+2^{D_u}C_v.                         \tag{1}
\]

This is inherited, not a new claim: **THM-4604**,
`01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md`.
That theorem also proves faithfulness of the full `(J,D,C)` word encoding.
The exact odd source cylinder is

\[
 n\equiv(2^D-C_w)3^{-J}\pmod {2^{D+1}}.                \tag{2}
\]

It has normalized odd-Haar mass `2^{-D}`. In (2) the right side is odd.
Successively applying the inverse odd branches shows that this cylinder has
exactly the stated valuations, not just an odd terminal value. Equivalently,
each inverse branch maps the odd 2-adics bijectively onto one valuation
cylinder. Delete poles and their backward iterates when discussing infinite
itineraries. Every finite word has positive integer representatives.

Define `Q_0=1`, `Q_j=1+sum_{i<=j}(3*1[2|a_i]-1)`. Let `T` consist of words
with `Q_j>0` before their last letter and `Q_J=0`. Then

\[
             \mathcal T=\mathcal L+\mathcal E\mathcal T^3,       \tag{3}
\]

where `L` is one odd valuation and `E` is one even valuation. The three
children are ordered. To parse a nonleaf word, after its even first letter
cut at the first successive queue levels `2,1,0`. Each segment is a member
of `T`. This proves unique parsing. The words are prefix-free.

**Representation warning.** This is the depth-first parse of the chronological
word. The prior breadth-first valuation-tree representation can have a
different shape. No valuation is permuted here. A conversion between the
two tree drawings must retain the chronological word and relabel its ports.
The invariants `(J,D,C)`, the source and the map do not change.

## 2. A triangular recurrence for every carry moment

Work in `Q[[u,z]]`, where `u` marks odd steps and `z` marks valuation clock.
For `r>=0`, define

\[
 F_r(u,z)=\sum_{w\in\mathcal T}u^{J_w}z^{D_w}C_w^r,
 \quad L={uz\over1-z^2},\quad E={uz^2\over1-z^2}.        \tag{4}
\]

These are **counting** series. Haar weights are introduced by replacing
`z` with `z/2`. Each coefficient is finite since `J<=D`. For any two word
languages, their moment-series concatenation is

\[
 (A*B)_r(u,z)=\sum_{i=0}^r {r\choose i}
       A_i(u,2^{r-i}z)B_{r-i}(3^iu,z).                 \tag{5}
\]

Proof: expand (1) to the `r`th power. The powers of `2^{D_u}` and
`3^{J_v}` produce exactly the two argument shifts. This product is
associative because word concatenation and affine composition are.

Write `H_r=(F*F*F)_r`. Explicitly,

\[
 H_r=\sum_{i+j+\ell=r}{r!\over i!j!\ell!}
 F_i(u,2^{j+\ell}z)F_j(3^iu,2^\ell z)F_\ell(3^{i+j}u,z). \tag{6}
\]

Equation (3) gives `F_0=L+EF_0^3`. For `r>=1`, put

\[
 P_r=F_0(3^r u,z),\quad Q_r=F_0(u,2^rz),\quad
 S_r=P_r^2+P_rQ_r+Q_r^2.
\]

Let `H_r^mix` be (6) with its three terms having index `r` omitted. All
its moment indices are then less than `r`. Separating those three terms
and the root letter in (5) gives the recurrence

\[
 F_r={L+E(u,2^rz)H_r^{mix}
       +\sum_{i=1}^r {r\choose i}E(u,2^{r-i}z)H_{r-i}(3^iu,z)
       \over 1-E(u,2^rz)S_r}.                         \tag{7}
\]

The denominator has constant term one, so division is formal and unique.
This determines the moments successively, using only lower moment orders
and argument shifts. In particular, with `P=F_0(3u,z)`, `Q=F_0(u,2z)`,

\[
 \boxed{F_1={P/3\over1-E(u,2z)(P^2+PQ+Q^2)}}.          \tag{8}
\]

Indeed `L+EP^3=P/3`. No asymptotic or probabilistic assumption is used in
(4)–(8). A scalar size recurrence misses the substitutions in (5), which
are the record of chronological affine order.

**Hostile.** Words `(2,1,1,3)` and `(2,1,3,1)` both have `(J,D)=(4,7)` and
are first-extinction words, but their integer carries are `103` and `151`.
They therefore have different source guards and different endpoints.

**What moments retain.** The complete carry distribution at each `(J,D)`
retains every word by faithfulness of `(J,D,C)`. All moments determine that
finite distribution, for example by a Vandermonde system after its finite
support is specified. Four moments alone do not do so in general. The
checked finite jets are exact summaries, not source-complete quotients.

## 3. Shared ports and ordered concatenation

Give each word a start port and an end port. Define `F_{pq,r}` by (4) with
these ports fixed. Formula (5) becomes

\[
 (A*B)_{pq,r}=\sum_h\sum_{i=0}^r {r\choose i}
 A_{ph,i}(u,2^{r-i}z)B_{hq,r-i}(3^iu,z).               \tag{9}
\]

This is an associative matrix version of the same operation. It is valid
for a specified finite port alphabet with composition exactly at matching
ports; sums count typed paths, including multiplicity when different port
histories have the same valuation word. For physical moves one must additionally impose the source guard,
pull back the second guard through the first map, and keep any occupied-pin
set or kernel lift that later legality tests use. Summing over `h` before
an attachment is chosen can merge physically incompatible histories.

Thus (9) supplies a recurrence template for shared-port grammars; it does
not prove that a proposed tournament grammar is relocation-closed or that
its ports may be reused. The incoming transporter hostile in
`complement_port_crt_20261028/PORT_LAWS.md`, P1–P2, explains why the full
endpoint relation must accompany a fixed-set count.

## 4. The hidden-coordinate interpretation

The old continuant theorem **THM-3186**,
`01-canon/theorems/THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall.md`,
propagates a hidden two-dimensional state and sums its contributions into
an exit coordinate. Its exit sum can cancel without the hidden state
vanishing. Here `C` is the off-diagonal coordinate of the affine matrix;
the eigenvalues alone retain only `(J,D)`. Equation (5) is the exact
translation of this keep-the-extension-coordinate mechanism to affine
word languages. It is not an identification of the two underlying
operators or of their proof targets.

**OPEN next step.** Construct a finite family of guarded, port-indexed
languages for which a carry-sensitive inequality survives (9) and implies
`Lambda*n+B<n` at the original positive source. Neither the size series
nor its first carry moment supplies that inequality.
