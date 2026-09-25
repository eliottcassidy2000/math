# Sewing arithmetic certificates: the nine-to-four bridge and square-reflection creation

**Status: PROVED elementary constructions, decoders, and scoped descent;
FINITE-EXACT controls; CITED deposited construction; OPEN universal Collatz
coverage.** Date: 2026-09-25. No novelty claim is made for affine Collatz
matrices, rational parametrization, or tensor contraction in general.

## 1. Source, inheritance, and live board

The user's DOI resolves to [He--Jing--Li, *The Symbols of Six-Gluon MHV
Amplitudes through Nine Loops*](https://doi.org/10.5281/zenodo.22800071),
a dataset deposited September17,2026. Its README constructs coefficients
by joining first-entry and last-entry recurrences with a sewing matrix.
The nine-loop symbol has weight18 and uses a12+6 split. Its nine-letter
alphabet comes in three triples; the bases carry separate parity and
dihedral representation labels. The symbol omits constants and other
information invisible to that representation. These are source properties,
not statements about integer Collatz trajectories.

The [source and group-action companion](zenodo_sewing_20260925.md) records
the exact definitions and an integral form of the dihedral action; source
integrity is recorded in section8 below. The archive was read as data; its Wolfram code was inspected
but not executed. No accompanying paper was identified by the bounded
source search, and no physical result beyond the deposit is assumed.

Closest proved mechanisms: the [guarded affine-word lattice model](creation_fano_20260925.md),
F3; [intrinsic elliptic creation](creation_elliptic_20260925.md), E1--E3;
and [marked primitive-triple normalization](ternary_triples_20260925.md).
The [previous decoder](creation_decoder_20260925.md), C9--C10, proves the
fixed-coordinate degree and finite-state itinerary obstructions. The
canonical hostile is an arbitrarily long initial odd run. The corrected
near miss is treating one-step binary creation as actual Collatz descent.
The least-used resource here is the arithmetic type of the two denominators
on opposite sides of a certificate seam.

Anchor: exact front/back certificate sewing. Niche: connect the actual
9→14→7 block to a primitive-triple parent operation. Wildcard: extend the
user's3,4,5 pattern to square-graph paths and the +2 recursion.

| Live object | Preserved structure | Required boundary |
|---|---|---|
| Source tensor | ordered letters and boundary basis contraction | no claimed transfer of physics constraints |
| Affine word | exact scale and ordered carry | integer endpoint must be retained |
| Primitive triangle | marked rational coordinate and parity gap | triangle size alone can lose the clock |
| Square path | edge square labels and vertex positions | positivity, loops, repeated vertices |
| Three-color action | integral dihedral lattice reduced mod2 | source parity disappears in the reduction |
| Creation decoder | unique word, finite reconstruction rank | word membership differs from root coverage |

The META-PATTERNS cards used are **Separate descent, ambient scale, and
regularity debt**, and **Separate unbounded local support from a
height-bounded modular cover**. Here they require, respectively, preserving
the odd-return clock and not turning a finite word census into all-height
coverage.

## 2. The local difficulty at nine is repaired by keeping the next halving

Use the shortcut map `T(n)=n/2` for even n and `(3n+1)/2` for odd n.
The actual block is

\[
9\longrightarrow14\longrightarrow7.
\tag{S1}
\]

The earlier9→14 example only refuted copying the binary-creation rank at
every individual shortcut step. It did not obstruct descent along a
correctly guarded longer block.

The [triple companion](zenodo_triples_20260925.md) proves that, on the odd
integer boundary,

\[
\Phi(n)=\left(n,\frac{n^2-1}{2},\frac{n^2+1}{2}\right)
\]

is a parabolic Berggren ray, of depth `(n-1)/2`. Thus (S1), at its two odd
endpoints, is exactly

\[
(9,40,41)\longrightarrow(7,24,25),
\]

one reverse step on that ray. The middle state is also representable
without discarding its parity: for every positive integer n use

\[
\Psi(n)=\frac{(2n,n^2-1,n^2+1)}{d_n},\qquad
d_n=\begin{cases}2&n\text{ odd},\\1&n\text{ even}.
\end{cases}
\]

The marked inverse is `n=A/(C-B)`, and the gap `C-B` is one or two.
The full block becomes

\[
(9,40,41)\to(28,195,197)\to(7,24,25).
\]

This gives the factor two in14 an actual coordinate. Merely identifying
14 with its odd part7 would lose the first arrow:14 goes to7, but7 goes
to11.

The block generalizes to the whole cylinder

\[
8a+1\longrightarrow12a+2\longrightarrow6a+1,
\qquad a\ge1.
\tag{S2}
\]

Its Berggren-ray depth falls from4a to3a. On the fruit elliptic subgroup,
write `P=8Q+G` with `Q=aG`; the same actual block gives

\[
\frac{3P+G}{4}=6Q+G.
\]

Both the triangle hypotenuse and elliptic canonical height decrease on
this cylinder. The exact proof and hostile controls in the companion
show why this is scoped descent, not coverage of every odd input.

There is also a specific reason25 appears in the user's chain. For odd
`y>0`, the equation

\[
2y+T(y)=\frac{y^2+1}{2}
\]

holds exactly at `y=7`: it reduces to `y(y-7)=0`. Hence
`14+11=25` joins twice the odd state7 and its successor11 into the
hypotenuse of `Phi(7)`. That triple is the Gaussian square of the
3-4-5 root, with `(4+3i)^2=7+24i`. The configuration is exact and locally
rigid; the equation is not a universal Collatz recurrence.

The junction does extend to a proved infinite family of complete root
certificates:

\[
n_r=\frac{7\cdot4^r-1}{3}=9,37,149,597,\ldots,\qquad r\ge1.
\]

Every member reaches7 after exactly `2r` shortcut steps and shares the
last transition14→7. The directly verified accelerated tail
`7→11→17→13→5→1` completes the certificate. This uses the inherited
inverse-fibre rule with its exact clock, now joined to the marked
square-hypotenuse configuration.

## 3. Two-channel sewing of actual Collatz words

This is a mathematical extension inspired by the deposit's method, not
an identification of its amplitudes with integers.

For a chosen sign `sigma=+1` or `-1`, encode a binary word w of shortcut
parities by

\[
M_w=\begin{pmatrix}A_w&B_w\\0&D_w\end{pmatrix},\qquad
A_w=3^{\#1(w)},\quad D_w=2^{|w|},
\]

so its formal endpoint is `(A_w n+B_w)/D_w`. Empty word means the identity.
Appending zero doubles D; appending one replaces `(A,B,D)` by
`(3A,3B+sigma D,2D)`. The carry B therefore retains the order of letters.

For a split `w=uv`, the matrix order is `M_w=M_v M_u`. Given integer
start n and target r, the two boundary values are

\[
L_u(n)=\frac{A_un+B_u}{D_u},\qquad
R_v(r)=\frac{D_vr-B_v}{A_v}.
\tag{S3}
\]

**S3a (PROVED).** They agree if and only if the word `uv` is an actual
legal trajectory from n to r. On positive integer inputs both signs
stay positive.

To see the first important mechanism, the left boundary lies in
`Z[1/2]` and the right in `Z[1/3]`. A rational whose reduced denominator
is both a power of two and a power of three is an integer. Thus equality
automatically gives an integer seam. The equality also gives the integer
endpoint of the complete word. The endpoint-integrality lemma, inherited
from F3 and independently rechecked here, then recovers every intermediate
parity guard. Explicitly, at the first letter the final numerator's parity
forces the declared parity of n, because all later multipliers are odd
and all later added offsets contain a factor two. Peel that step and
repeat. The converse follows from actual iteration.

**S3b (PROVED).** The residual of the sewing condition factors through
exactly two arithmetic coordinates:

\[
K_{n,r}(u,v)
=A_v(A_un+B_u)+D_u(B_v-D_vr)
=\begin{pmatrix}A_un+B_u&D_u\end{pmatrix}
 \begin{pmatrix}A_v\\B_v-D_vr\end{pmatrix}.
\tag{S4}
\]

For fixed positive n,r and lengths `h=|u|>=1,k=|v|>=1`, the matrix with
all prefixes as rows and suffixes as columns has rank exactly two over Q,
for either sign. The factorization bounds rank by two. Choose the all-zero
word and the word consisting of one followed by zeros on each side.
The determinant of the resulting two-by-two minor is

\[
-2^h(2n+\sigma)(\sigma+2^{k+1}r)\ne0.
\]

This proves equality. The two coordinates are unbounded integer
registers, not two finite states. It therefore does not conflict with
the earlier infinite-residual-state obstruction for the full itinerary
encoder.

For the user's9, there is a short explicit root certificate:

```text
9,14,7,11,17,26,13,20,10,5,8,4
prefix u = 10
suffix v = 111010010
M_u = (A,B,D) = (3,1,4)
M_v = (A,B,D) = (243,347,512)
```

The seam is exactly

\[
\frac{3\cdot9+1}{4}=7=
\frac{512\cdot4-347}{243}.
\tag{S5}
\]

Every number, division, and parity is certified by this equality and
the word-matrix definitions. The target4 uses the shortcut clock here;
its terminal1↔2 cycle does not itself revisit4. Thus the all-input
coverage question with target4 is posed for `n>2`, with1 and2 handled
separately (or with terminal bank `{1,2,4}`).

## 4. A new creation object with an intrinsic, terminating decoder

**S4 (PROVED).** The two elementary matrices freely generate their word
monoid, and word membership has a coordinate-only decoder.

The generators are

\[
M_0=\begin{pmatrix}1&0\\0&2\end{pmatrix},\qquad
M_1=\begin{pmatrix}3&\sigma\\0&2\end{pmatrix}.
\]

For an input matrix represented by integers `(A,B,D)`, require positive
odd A and positive power-of-two D. If `D>1`, the first letter is exactly
`e=B mod2`. Peel it by the rule

\[
e=0:\ (A,B,D)\mapsto(A,B/2,D/2),
\]
\[
e=1:\ (A,B,D)\mapsto
\left(A/3,\frac{B-\sigma A/3}{2},D/2\right),
\tag{S6}
\]

where the odd case requires `3|A`. At `D=1`, accept exactly `(1,0,1)`.
Negative B is allowed for the minus sign.

Proof: write `w=e v`, so `M_w=M_v M_e`. For first letter zero,
`B=2B_v`; for first letter one, `B=sigma A_v+2B_v`. Since A_v is odd,
their parities determine e uniquely. These product equations give (S6).
The integer rank `log2 D` decreases by one, so decoding either rejects
or reaches the identity in finitely many steps. Reversing the accepted
steps reconstructs the input matrix, proving sufficiency as well as
necessity. Uniqueness of every peeled letter proves freeness.

This is a creation object whose structure retains the entire ordered
arithmetic instruction. Its matrix dimension stays fixed while the
coefficient heights grow. An endpoint equation

\[
A_wn+B_w=4D_w
\tag{S7}
\]

then turns an accepted creation object into a root certificate for `n>2`.
Constructing such an object for every n is the remaining **OPEN** coverage
obligation. The decoder's termination does not prove that (S7) is solvable
for every n.

A minimal order-loss hostile is `01` versus `10`. Both have one odd and
one even step, but their matrices are `(3,2,4)` and `(3,1,4)`. At n9 the
first formal endpoint is29/4, while the second is the legal integer7.
Neither odd-step counts nor scale alone can replace B.

## 5. Primitive triples create square-sum paths and the +2 recursion

The [square-path companion](zenodo_square_20260925.md) proves an injective
construction from every primitive triple `(a,b,c)` with a odd and b even:

\[
a^2,\quad b^2,\quad a^2+2c+1,\quad
b^2-2c-1,\quad a^2+2,\quad b^2-2.
\tag{S8}
\]

These are six distinct positive vertices below `c^2`. Their successive
edge sums are

\[
c^2,\ (c+1)^2,\ c^2,\ (c-1)^2,\ c^2.
\]

For the root `(3,4,5)`, this gives

\[
9-16-20-5-11-14.
\]

Thus5,11,14 really belong to a square-edge path together with9. The
identity9+5=14 is a separate balance relation;9--5 itself is not a
square edge. The root uniquely satisfies the corresponding balance
`(c+b)+c=2(a+b)` among primitive triples. This preserves the user's
observation with its predicate stated precisely.

Let `R_s(x)=s^2-x`, a reflection across a square-sum edge whenever both
vertices are positive and distinct. Four such reflections satisfy

\[
R_{c-1}R_cR_{c+1}R_c(x)=x+2.
\tag{S9}
\]

The constant is the second difference of the square function. More
generally, replacing1 by d gives `x+2d^2`. This directly realizes the
earlier +2 recursion as a sequence of square edges.

Starting at9 with c5, the resulting simple path inside `1,...,24` is

```text
9,16,20,5,11,14,22,3,13,12,24,1,15,10.
```

It stops being a path in that graph when the next vertex is26. For other
seeds a repeated vertex can occur first; the companion proves the exact
positivity/collision boundary. In particular the construction is not a
Hamiltonian-path proof and does not stay within the class of PPT seeds:
the new starting coordinate `a^2+2` is not a square.

## 6. The deposited dihedral action really meets the Fibonacci colors

The source uses a rational two-dimensional representation of D3. The
[source companion](zenodo_sewing_20260925.md) changes to a preserved
integral lattice, where the generators become

\[
C=\begin{pmatrix}0&-1\\1&-1\end{pmatrix},\qquad
F=\begin{pmatrix}0&-1\\-1&0\end{pmatrix}.
\]

They preserve `a^2-ab+b^2`. Reducing this integral model modulo two gives
the Fibonacci color rotation `(a,b)->(b,a+b)` and a reflection swapping
the two coordinates. This faithfully realizes `D3=GL2(F2)` on the three
nonzero colors, using the exact convention already repaired in
[the Zeckendorf charge note](duck_zeckendorf_20260925.md).

The source's additional parity involution acts as scalar minus identity
on an odd block. It becomes invisible modulo two. Therefore this map
preserves the dihedral color action but loses that separate sign; it
does not identify source parity with odd/even integer parity. Keeping
the sign as an additional label restores the distinction.

## 7. What transfers, and what remains to construct

| Source -> target | Exact map/predicate | Lost information and repair |
|---|---|---|
| Deposited coefficient assembly -> arithmetic certificates | ordered left/right recurrences joined at a boundary | analogy of method; no physical identity is asserted |
| Forward/reverse arithmetic words -> integer seam | equality in `Z[1/2]` and `Z[1/3]` | retain exact numerator and denominator, not approximate ratios |
| Word -> affine matrix | free monoid encoding and terminating peeling | no ordered-word loss; dropping B destroys the guarantee |
| Actual9→14→7 -> triangle ray | marked coordinate and odd first return | retain parity gap or explicit return clock |
| PPT -> square path | formula(S8), exact square edge labels | iteration leaves PPT seeds; keep positivity and collisions |
| Source D3 -> three colors | integral lattice then reduction modulo two | separate parity sign disappears |

The next useful extension is an adaptive family of accepted word matrices
whose endpoint equation covers every `n>2`, with a well-founded rule for
joining shorter certificates. The source suggests how to keep separate
boundary data and contract only what is needed. It does not provide the
arithmetic selection theorem. The n≡1mod8 family closes a real infinite
portion of that task; the other residue families and unbounded growing
prefixes remain obligations.

## 8. Reproduction and audit

The downloaded archive `hexagon_mhv_symbol_weight_18.7z` has16,558,941
bytes and matches the Zenodo API checksum
`md5:5e1ddc5dc356c4ed1b46309896fd3cbb`. All six extracted textual release
files match their entries in `checksums.sha256`: README, loader,
manifest, and the three TSV tables. The26 WXF binaries were not
extracted or evaluated. The README's SHA256 is
`b158e8dac72f868a47b62a55e4e78bd488071d27a057904e175fd62413ab2c3f`;
the loader's is
`c185c15247aabed370f1190efb9ce34ec8df6b19adb252393b96eb8434f00df0`.
The published tensor table declares the final sewing shape3224×1858;
this is a reported dimension, not a recomputed tensor rank. These checks
establish source identity and textual integrity, not physical correctness.

```text
python 04-computation/experiments/zenodo_bridge_20260925.py
python -O 04-computation/experiments/zenodo_bridge_20260925.py
```

The [script](../../04-computation/experiments/zenodo_bridge_20260925.py)
uses exact integers/Fractions and explicit exceptions. Its
[output](zenodo_bridge_20260925.out) records16,382 word decodes and
injectivity controls (both signs, every word of length0..12);65,408
endpoint/intermediate parity checks (length0..8,n1..64);36,912
front/back sewing checks including incorrect targets;288 nonzero
rank-two minors; and five deliberately rejected matrices. Direct
stepwise parity replay is independent of the affine integrality check.
The all-depth sewing, rank, and peeling arguments were independently
audited by a second agent.

The three companions supply independent bounded controls for the
triangle/ray formulas, square paths and exact collision boundary, and
source representation matrices. Their full finite universes are stated
there. Finite experiments are not used as substitutes for the proofs
or for universal Collatz coverage.
