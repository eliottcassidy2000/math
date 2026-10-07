# A complete decoder for a fixed two-anchor head template

**PROVED:** uniqueness/existence validation for a prescribed positive-word
carrier, complete finite gap enumeration for each supplied head, and the
native common-future receipt export. **FINITE-EXACT:** the stated head
search, antichains, and controls. **CONDITIONAL:** ROOT completion consumes
a supplied first-hit child proof. No universal head success or Mersenne
coverage follows from the decoder alone.

## 1. Inheritance and the exact question

The carry decoder is inherited from
[reciprocal four-channel kernel, section 5](reciprocal_four_channel_kernel_20261004.md).
The affine template comes from
[THM-4601, two-anchor reduction](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md),
and the current independently checked
[completion-anchor compiler](collatz_completion_anchor_20261007b.md).
[The collision grammar](collatz_collision_dp_20261007.md) searches a
different family of word relations. This package does not claim its
fixed-template decoder lists all common-future partners.

The closest mechanism is positive odd carry peeling. The canonical
hostile is the formal value `F_(10)(1)=1/256`: matching at the anchor 1
does not authenticate that word at integer ROOT. The corrected near miss
is to increase the short head's exponent 10 while expecting its partner
to persist. The least-used sidecar is the **terminal cost**, since the
carry alone does not encode the final valuation. The live board is the
supplied head, ordered carry, remaining cost, sibling gap, native cell,
and supplied child certificate.

For a chronological positive valuation word `w`, write

\[
F_w(x)=\frac{3^{|w|}x+B_w}{2^{\sum w}}.
\]

Fix a head `u` of length `ell`, cost `A`, and carry `B`. For a positive
integer gap `r`, request precisely

\[
|v|=\ell+3,
\qquad \sum v=A-2r,
\qquad F_v(1)=S^r(F_u(1)),\qquad S(x)=4x+1.
\tag{1}
\]

The empty head has no such partner. This is a formal affine equation at
1; the actual source is supplied separately or constructed by section 4.

## 2. The unique decoder and complete finite gap list

For a nonempty word of length `L>=2`, splitting its first letter gives

\[
B_w=3^{L-1}+2^{a_1}B_{\mathrm{tail}},
\tag{2}
\]

where the suffix carry is positive and odd. Therefore its first letter
must be `a1=v2(B_w-3^(L-1))`. Subtract the power of three, divide by
`2^a1`, decrease length and remaining cost, and repeat. A nonpositive
difference, zero valuation, or insufficient remaining positive-letter
cost rejects the data. At length one the carry must equal 1, and the
remaining cost must be positive; that cost is the final letter.

**The decoder is an iff.** Every actual positive word passes each forced
step by (2). Conversely, accepted data reconstruct a positive word with
the specified length, cost, and carry by reversing (2). Thus it is unique.
The algorithm makes at most `L-1` peeling decisions and uses integer
arithmetic throughout. This is linear in the number of letters as an
arithmetic-operation count, not a linear bit-complexity or bounded-memory
claim: the retained cost, powers, and carry have unbounded bit size.

Since `S^r(x)=4^r*x+(4^r-1)/3`, put `C=A-2r` and `L=ell+3`. Equation (1)
is equivalent to the target carry

\[
B_v=B-26\,3^\ell+2^C\frac{4^r-1}{3}.
\tag{3}
\]

In particular, at `r=1` it is `B-26*3^ell+2^(A-2)`. Apply the decoder to
`(L,C,B_v)`. There is no positive partner if it rejects. All positive
gaps are covered by the finite interval

\[
1\le r\le \left\lfloor\frac{A-\ell-3}{2}\right\rfloor,
\tag{4}
\]

because `C>=L` is necessary. `all_partners` tests exactly (4), so it is
complete for (1) at arbitrary fixed supplied `u`. There is at most one
partner per gap, not necessarily one partner across all gaps.

### Unrestricted short-head rigidity

For `u=(a)` at gap one, write the partner's increasing partial costs as
`0<s1<s2<s3<C=a-2`. Equation (3) becomes

\[
2^C=104+9\,2^{s_1}+3\,2^{s_2}+2^{s_3}.
\]

The next partial cost must always equal the valuation of the accumulated
constant: otherwise the right side has a unique term of least valuation,
strictly below `C`. The successive constants are
`104 -> 176 -> 224 -> 256`, of valuations `3,4,5,8`. Thus the unique
solution is `a=10`, `v=(3,1,1,3)`.

For `u=(1,a)`, the corresponding equation is

\[
2^C=310+27\,2^{s_1}+9\,2^{s_2}+3\,2^{s_3}+2^{s_4},
\quad C=a-1.
\]

Now `310 -> 364 -> 400 -> 448 -> 512` forces partial costs `1,2,4,6`
and total cost 9. The unique solution is `a=10`, `v=(1,1,2,2,3)`.
These proofs have no cost cutoff. Their finite controls independently
enumerate all 53,130 child compositions of those lengths through cost 24.
Increasing just the 10 cannot extend either short-head shape.

## 3. A bounded search with useful new heads

The declared source-head universe is **all 22,620 words of lengths 1–4
over the alphabet 1–12**. Every gap in (4) is considered. There are 274
pairs, of which 269 have gap one:

| Head length | Gap 1 | Gap 2 | Gap 3 |
|---|---:|---:|---:|
| 1 | 1 | 0 | 0 |
| 2 | 7 | 0 | 0 |
| 3 | 40 | 0 | 1 |
| 4 | 221 | 3 | 1 |

Representative nonpadding heads are:

| `u` | `r` | `v` |
|---|---:|---|
| `(3,12)` | 1 | `(4,1,5,2,1)` |
| `(4,12)` | 1 | `(3,6,1,3,1)` |
| `(5,8)` | 1 | `(3,1,3,3,1)` |
| `(8,4)` | 1 | `(3,1,1,4,1)` |
| `(9,2)` | 1 | `(3,1,1,3,1)` |
| `(1,2,9)` | 3 | `(1,1,1,1,1,1)` |
| `(1,3,5,5)` | 2 | `(1,1,3,1,1,1,2)` |
| `(4,1,5,5)` | 2 | `(3,2,1,1,1,1,2)` |
| `(6,1,4,3)` | 2 | `(3,1,1,1,1,1,2)` |

A leading 2 is redundant for **anchor evaluation** because `F_2(1)=1`.
If a head `(2)+u` has a partner, uniqueness shows that partner is
`(2)+v` for the same gap. The inverse implication follows immediately by
prefixing 2 to both words. This does not make the letter 2 an identity
on actual source integers or allow it to be discarded from a receipt.

### Exact finite native guard union

Discard heads starting with 2 only for the following comparison. Among
the remaining successes, retain the prefix-minimal heads. In this finite
universe all 221 gap-one heads already form an antichain; the all-gap
version has 225 heads and is likewise prefix-free. This is a finite
statement, not a theorem about every successful head.
`finite_head_bank()` exposes these 225 authenticated packets in deterministic
length/lexicographic head order. In this universe each head has only one
successful gap; the function specifies the least gap if a head had several.

For odd Haar sources, an actual head of cost `A` has mass `2^-A`. Distinct
prefix-free valuation words have disjoint native cylinders. Conditioning
on first valuation 1 gives weight `2^(1-A)`; conditioning on first
valuation at least 3 gives weight `2^(2-A)`. Thus:

| Bank | First valuation 1: count / mass | First valuation >=3: count / mass |
|---|---|---|
| Gap 1 | `33`, `1913/524288` | `188`, `139281/8388608` |
| Every allowed gap | `35`, `2233/524288` | `190`, `142353/8388608` |

The singleton heads `(1,10)` and `(10)` have corresponding masses
`1/1024` and `1/256`, so both branches increase strictly in this finite
native-head comparison. The cost multiplicities are:

| Cost | Gap1, first1 | Gap1, first>=3 | All gaps, first1 | All gaps, first>=3 |
|---:|---:|---:|---:|---:|
| 10 | 0 | 1 | 0 | 1 |
| 11 | 1 | 1 | 1 | 1 |
| 12 | 1 | 2 | 2 | 2 |
| 13 | 3 | 5 | 3 | 5 |
| 14 | 6 | 8 | 7 | 9 |
| 15 | 8 | 14 | 8 | 15 |
| 16 | 5 | 19 | 5 | 19 |
| 17 | 3 | 27 | 3 | 27 |
| 18 | 3 | 26 | 3 | 26 |
| 19 | 2 | 23 | 2 | 23 |
| 20 | 1 | 23 | 1 | 23 |
| 21 | 0 | 15 | 0 | 15 |
| 22 | 0 | 11 | 0 | 11 |
| 23 | 0 | 8 | 0 | 8 |
| 24 | 0 | 4 | 0 | 4 |
| 25 | 0 | 1 | 0 | 1 |

These exact masses refer to native head cylinders and are also their
finite-union natural relative densities among odd sources. Transferring
them to a Mersenne shell requires the separate exponent-phase map, maximal
run convention, first-hit proof, and payment height cut. No such transfer
is inferred merely from the displayed Kraft sums.

## 4. An all-height actual receipt interface

Let `(u,v,r)` pass the decoder. Append terminal valuations `c` and `c+2r`.
The resulting words have equal denominator `Q=2^(A+c)` and equal formal
value at 1. Their slopes differ by 27. Writing their numerators as
`P*x+B_left` and `27P*y+B_right`, respectively, gives

\[
B_{\rm right}=B_{\rm left}-26P,
\qquad
F_{u,c}(27y-26)=F_{v,c+2r}(y).
\tag{5}
\]

The export fixes **c=1**, so both terminal valuations are odd. If
`(p,Q,b)` is the carrier of `v+(1+2r,)`, its exact odd-source cell is

\[
y=y_0+2Qt,\quad
y_0\equiv(Q-b)p^{-1}\pmod{2Q},\quad 0<y_0<2Q,
\quad t\ge0.
\tag{6}
\]

The standard final-oddness argument authenticates every listed valuation:
backwards from the odd endpoint, positivity/integrality forces the
successive exact dyadic divisions. Equivalently, this follows by induction
through the nested word cylinders. Equation (5) then authenticates the
left word too, at `x=27y-26`.

There are no hidden ROOT loops. If any earlier state were 1, every
remaining actual valuation would be 2; an odd final valuation contradicts
this. A final ROOT hit also requires an even valuation, since
`(2^a-1)/3` is integral only for even `a`. Thus both actual paths avoid
ROOT entirely, `y0>=3`, their common endpoint exceeds 1, and `x>y`.
Every `t>=0` in (6) therefore gives an authenticated strictly smaller-child
common-future receipt. These sources satisfy `x=1 mod27`, and the endpoint
identity has not been claimed for arbitrary supplied source pairs.

`native_receipt` returns the existing exact receipt type; `discharge`
consumes a separately supplied first-hit ROOT word for that same `y`,
checks that its actual prefix agrees, and splices the common-endpoint
suffix onto the left word. It never discovers the child's fate.

For `u=(10)`, gap one, the least pair is

\[
63829\xrightarrow{(10,1)}281
\xleftarrow{(3,1,1,3,3)}2365.
\]

The stored child word
`(3,1,1,3,3,2,1,3,1,1,3,4,1,3,1,2,3,4)` is independently replayed to
first ROOT. It exports a 15-edge ROOT word for 63829. An empty suffix,
root-padded suffix, or the same certificate attached to another parameter
is rejected. This example checks the interface; it is not a claim of
new coverage relative to every inherited Collatz rule bank.

## 5. Reproduction and remaining obligation

```text
python -B -X utf8 04-computation/experiments/collatz_twoanchor_head_decoder_20261007b.py
python -B -O -X utf8 04-computation/experiments/collatz_twoanchor_head_decoder_20261007b.py
```

An independent table contains every 26,332 positive word with length at
most seven and cost at most 16. It is constructed by composition enumeration,
without the decoder. Every entry round-trips, and all 30,464 triples with
length 1–7, cost 0–16, and carry 0–255 are compared against that table.
The head universe is as stated in section 3; every available partner is
re-encoded and tested through independent affine identities. All 274 pairs
have two native exports literally replayed, totaling 548 instances.
The short-head rigidity controls, exact source and packet types, impossible
carriers, and supplied-child ROOT boundary have explicit hostiles.

Normal and optimized modes perform **195,307 explicit checks**. The
finite search proves only its declared counts; equations (1)–(6) prove the
unrestricted decoder and native-export statements. The next obligation
for a prescribed Mersenne source is its independently checked phase and
height interface, not a longer blind partner search. Failure of this
decoder rules out only the specified length/cost/sibling template.

Saved normalized-LF SHA256:
`357dd08ce163bea259f0d472b5cbe3f95a0517f7d4485ada47ab8142aaefebcf`.
