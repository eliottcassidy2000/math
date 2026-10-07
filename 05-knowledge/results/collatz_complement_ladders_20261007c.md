# Signed ladders in the complement of the finite head bank

**Status: PROVED** for the marked carrier decoder, signed-gap bounds, zero-gap
normalization, exact reserve guard and guarded common-future exports below.
**FINITE-EXACT** for the stated head universe and controls. **CONDITIONAL** for
ROOT completion: an authenticated smaller-child ROOT word is still required.
No claim covers every positive integer or the entire finite-bank complement.

## 1. Inheritance and the changed operation

The closest mechanism is the positive-gap unique carry decoder in
[the two-anchor head package](collatz_twoanchor_head_decoder_20261007b.md).
Its 225 prefix-minimal heads remain unchanged. The correction in
[THM-4601, two-anchor reduction of the residual first-reset-two branch](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md)
and `04-computation/experiments/twoanchor_20261007/audit_A/` already identifies
source-side ladders, including `(12,1,1)`, and distinguishes the short J=2
state. This note does not claim discovery of that head. It supplies its exact
additional valuation reserve, transfers it to a native Mersenne phase, and
proves that this guarded family is outside the unchanged 225-head bank.

The canonical hostile is an actual `(12,1,1,1)` source: the source head is
legal, but its last valuation cannot pay a source-side gap of one. The
corrected near miss is counting all zero-gap endpoint equalities as new
rules; Section 3 recovers their signed core. The least-used sidecar is the
marked affine offset: J=2 has `(k,d)=(4,10)`, while J>=3 has `(3,-26)`.
They must not share a decoder state accidentally.

The live concepts are full affine carry, signed sibling displacement,
terminal valuation reserve, immutable-source growth and exact exponent phase.

| Source | Target/map | Preserved predicate | Lost without sidecar | Required sidecar |
|---|---|---|---|---|
| Two affine heads | Unique carry decoder | Exact marked affine equality | Native integer legality | Ordered words and source cylinder |
| Two sibling ladders | Cancel their common powers of S | Same endpoint relation | Which side supplies binary division | Signed net gap |
| Native source head | Coarse reserve cell | Head plus sufficient last valuation | Actual last valuation value | Recompute c at supplied source |
| Residual integer | Exponent or odd-cofactor phase | Entire native prefix | First-hit safety and payment | Growth cut, deletion depth, original source |
| Paid common future | Supplied child proof splice | Literal first-hit ROOT certificate | Existence of that child proof | Authenticated child word |

These are partial arithmetic maps, not a tournament or an unguarded group
action. No literature theorem is needed for the new proofs.

## 2. Complete signed-gap decoder at a fixed marked state

For a positive valuation word w, write

\[
F_w(z)=(3^{|w|}z+B_w)/2^{A_w},\qquad A_w=\sum w.
\]

Allow the empty word with carrier `(1,1,0)`. Fix a marked relation

\[
X=3^kY+d,\qquad k\ge1,\quad d\in2\mathbb Z,
\]

and a supplied head u with length ell, cost A and carrier `(P,Q,B)`.
Put `S(z)=4z+1`. For any signed integer r, the identity

\[
F_v(Y)=S^r(F_u(X))                                      \tag{1}
\]

forces and is forced by

\[
L=|v|=\ell+k,\quad C=A_v=A-2r,\quad
B_v=B+Pd+2^C\frac{4^r-1}{3}.                           \tag{2}
\]

The signed expression is integral. For `r=-i<0` its last term is
`-Q(4^i-1)/3`. Equality of slopes gives L and C; equality of intercepts gives
the remaining formula. The inherited exact decoder recovers at most one v
from `(L,C,B_v)`: for L>1,

\[
B_v-3^{L-1}=2^{v_1}B_{\rm tail},
\]

with positive odd tail carry. Thus v1 is the exact 2-valuation of the positive
difference. At L=1 the carry must be 1 and the remaining positive cost is the
last letter. The implementation rejects nonpositive differences, impossible
remaining costs and malformed exact-integer fields; it re-encodes packets.

There are only finitely many signed gaps for each supplied head, with no
arbitrary negative cutoff:

* For r>0, positivity of all L letters requires `2r<=A-ell-k`.
* For r=-i<=0, the necessary carry bound is
  `B+Pd-Q(4^i-1)/3>=1`. This strictly decreases with i and eventually fails.

The compiler tests precisely these finite ranges. Completeness is for (1)
with this fixed marked state and positive valuation words, not for arbitrary
common-future diagrams. Placing sibling ladders on both sides introduces no
second independent parameter: cancel their common power of the invertible
affine map S, retaining the signed difference r.

## 3. Zero-gap normalization

Suppose (1) has r=0. Remove the maximal identical terminal suffix from u,v.
The remaining source word cannot be empty: equal total costs, after removing
equal letters, would require the remaining nonempty v of extra length k to
have cost zero. Let their now-distinct final letters be a,b.

Before those letters the maps have dyadic rational coefficients. Modulo 3,
each numerator `3F+1` is 1. Equality after division by `2^a,2^b` therefore
forces `a=b mod2`. Removing the two final letters gives the signed core

\[
F_{v^-}(Y)=S^{(b-a)/2}(F_{u^-}(X)),\qquad (b-a)/2\ne0.   \tag{3}
\]

The stripped suffix and the final source letter a reconstruct the original
packet uniquely. Thus a zero-gap match is a completed signed ladder with
possibly a shared continuation; its mere presence is not additional guard
coverage. The finite control `(10,1,4)` reduces to source head `(10)`, gap1,
then source letter1 and common suffix4.

## 4. The exact reserve guard and safe native exports

Set `e=max(0,-2r)`, `Aeff=A+e`. In the formal valuation domain that permits
the continuation U(1)=1, a positive odd X supports the head u and a next
exponent c with `c>=e+1` if and only if

\[
\boxed{X\equiv -(3B+Q)(3P)^{-1}\pmod{2^{Aeff+1}}.}       \tag{4}
\]

Indeed this is divisibility of `3PX+3B+Q` by `2^{A+e+1}`. Reducing it
modulo `2^{A+1}` forces the head endpoint `(PX+B)/Q` to be odd; final oddness
forces every intermediate valuation to be exact. A nonintegral dyadic
intermediate cannot regain integrality, and an even one makes the following
state nonintegral. Positivity is retained by each positive affine letter.
The remaining divisibility is exactly `v2(3F_u(X)+1)>=e+1`. This iff is not
yet a strict first-hit receipt: the valid packet u=(10), gap1 includes X341
in its coarse cell, but that head reaches1 and its next formal exponent2
would be ROOT padding. The script retains this hostile. The odd-terminal
construction below and the residual growth cut in Section5 separately
exclude earlier ROOT; neither production export accepts that padding.

For gap `-i`, merely knowing u is native is insufficient: division through
`S^{-i}` needs `c>=2i+1`. The actual child final exponent is `c+2r`.
Equation (4) retains this unbounded final-valuation tail instead of selecting
one finite value. The marked ternary relation between X and Y is still a
separate necessary input.

For a generic native export take `c=e+1`. Both c and c+2r are positive odd.
Choose Y in the unique native cylinder of `v,(c+2r)`, with

\[
Y\ge\max\left(3,\left\lfloor\frac{-d}{3^k-1}\right\rfloor+1\right).
\]

Then `X=3^kY+d>Y>=3`, and the source word `u,(c)` has the same odd endpoint.
The carrier identity and final oddness prove literal legality on both sides.
A positive integer reaching ROOT at the final step would require an even
last valuation, since `2^a=1 mod3`. An earlier ROOT would force every remaining
valuation to be 2. The odd final letter therefore excludes both cases.
This yields a strict paid common-future receipt, not a ROOT suffix.

`discharge(receipt, supplied_child_word)` independently replays a supplied
first-hit child ROOT word and checks its child prefix before splicing. It
never searches that child's orbit. A frozen inherited word for2365 grounds
the receipt `63829 -> 281 <- 2365`; wrong-source and empty words are rejected.

## 5. Marked residual transfer, including the missing J=2 type

Let `n=2^K t-1`, K>=5, t positive odd. Its first K-1 valuations are ones.
After J subsequent twos the formal source state is

\[
X=1+\frac{3^J(3^{K-1}t-1)}{2^{2J-1}}.                  \tag{5}
\]

For deletion depth D=3 or4 set `h_D=2^{K-D}t-1<n`. The independently checked
clearing calculation gives:

| Two-run length | Child clearing after its initial ones | State relation |
|---|---|---|
| J>=3, D=3 | `(4,1,1),2^(J-3)` | X=27Y-26 |
| J>=3, D=4 | `(2,2,1,1),2^(J-3)` | X=27Y-26 |
| J=2, D=3 | `(4)` | X=81Y+10 |
| J=2, D=4 | `(2,2)` | X=81Y+10 |

Here `2^q` inside a word means q repeated letters2, not an exponent value.
The script rejects applying the J>=3 relation to J=2. Require a nonempty head
with first letter not2; it then terminates the maximal two-run exactly.

Let x0 denote the residue in (4). For the Mersenne line t=1, (4) becomes

\[
\boxed{3^{K-1}\equiv1+2^{2J-1}(x_0-1)3^{-J}
                \pmod{2^{2J+Aeff}}.}                  \tag{6}
\]

The target is1 mod8 and lies in the subgroup generated by3. This follows
elementarily from `v2(3^b-1)=2+v2(b)` for positive even b: the even powers
exhaust the1 mod8 subgroup. Thus (6) has one exact exponent class of minimal
period `2^(2J+Aeff-2)`. No huge integer needs to be expanded to check it.

The first letter fixes `v2(x0-1)=1` for u1=1 and 2 for u1>=3. Consequently

\[
v_2(K-1)=\begin{cases}2J-2,&u_1=1,\\2J-1,&u_1\ge3,\end{cases}
\quad
\text{relative phase density}=
\begin{cases}2^{1-Aeff},&u_1=1,\\2^{2-Aeff},&u_1\ge3.\end{cases} \tag{7}
\]

These are natural densities of exponent parameters within the stated shell,
not densities of Mersenne integers. For general odd t, multiply the right side
of (6) by `3^{-(K-1)}` to obtain one native residue for t at the same binary
precision. This supplies small exact receipt controls without expanding a
Mersenne number of hundreds of thousands of bits.

Retain the sufficient immutable-source cutoff

\[
K-1\ge2Aeff+J.                                        \tag{8}
\]

Every prefix through u then has coefficient greater than1 relative to n.
For a head prefix cost Ai<=A<=Aeff, its coefficient is at least

\[
\frac{3^{K-1+J}}{2^{K-1+2J+Aeff}}
\ge(9/8)^{Aeff+J}>1.
\]

The initial ones grow, and intermediate twos have still larger coefficients.
All nonempty carries are positive, so the head endpoint Z exceeds n. If
r=-i, (8) also gives `n>S^i(1)`, hence `S^{-i}(Z)>1` whenever it is integral.
The reserve (4) guarantees precisely its positive odd integrality.

Set `c=v2(3Z+1)`. The actual source word is
`1^(K-1),2^J,u,c`; the child word is its initial `1^(K-D-1)`, the appropriate
clearing row, v and `c+2r`. Their final odd endpoints agree. Both preterminal
endpoints exceed1, so no earlier ROOT occurs. If the shared final endpoint is
ROOT, that endpoint is the first one. This proves actual paid deletion
receipts for both D=3,4 on the retained native parameters. Removing the finite
head below (8) for a fixed phase does not change (7).

## 6. A source-side family genuinely outside the 225-head bank

The inherited signed core

\[
u=(12,1,1),\quad v=(3,1,1,3,2,6),\quad r=-1
\]

has cost A14 but needs effective cost Aeff16. Its exact coarse reserve is

\[
\boxed{X\equiv75093\pmod{131072}.}                     \tag{9}
\]

No head in the unchanged225 bank is prefix-comparable with `(12,1,1)`.
Native positive-word cylinders overlap only when one word is a prefix of the
other. Hence (9) is entirely disjoint from their union. Its relative mass in
the first-letter>=3 branch is `2^(2-16)=1/16384`. This is added guarded
head/phase coverage, with the same child obligation retained. It is not
unconditional ROOT coverage. The J3 Mersenne phase is

\[
\boxed{K\equiv485217\pmod{1048576}.}                   \tag{10}
\]

In contrast an actual source with word `(12,1,1,1)` fails (9); attaching the
same negative-gap partner would require the impossible final exponent-1.

For the changed J2 relation the exact decoder gives

| Source head | Partner, gap1 | Mersenne exponent phase | Phase modulo16 |
|---|---|---|---|
| `(10)` | `(1,2,1,1,3)` |1289 mod4096 |9 |
| `(1,12)` | `(2,1,1,2,4,1)` |27861 mod32768 |5 |
| `(1,1,14)` | `(3,1,1,1,4,3,1)` |177949 mod262144 |13 |

Only the third row addresses emitted exponents13 mod16. The first two are
useful hostiles to forgetting the actual phase. They must not be offered as
grounding that whole branch. The separately owned two-twos package handles
its source-specific grounded examples; this package authenticates the
changed relation and native row.

## 7. Reproduction and exact scope

Run `python -B -X utf8 04-computation/experiments/collatz_complement_ladders_20261007c.py`
and repeat with `-O`. The saved output records3,261 exact controls.

The fixed nonpadding source universe is letters1..12, lengths1..3. Successful
signed packets have counts `r=-1:1, r=0:228, r=1:40, r=3:1`. These are
packet counts, not disjoint coverage counts. All228 zero-gap packets are
normalized by (3). An independent enumeration of every positive composition
with length<=7 and cost<=13 checks349 marked signed carrier instances.
There are160 literal residual receipt controls, symbolic phase and wrong
phase controls, the missing-reserve hostile, and exact-type/forged-mark
rejections. The separate J2 `(1,1,14)` row lies outside the stated1..12 head
census and is clearly a selected exact control.

The first source-side phase is symbolic: no new ROOT search for
`2^485217-1` is claimed. No sampled count is promoted to infinite coverage.
The all-height laws come from Sections2–5; the finite computation tests their
interfaces and the declared finite-bank comparison.
