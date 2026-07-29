---
id: THM-2960
title: "Local Smith-jet Fitting barcode and negative-depth chamber atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A general
  DVR lemma identifies truncated block-Toeplitz
  nullity growth with local Smith bars and maximal-minor Fitting
  valuation.  On the finite first-gap bank (0,1,2,M), 6<=M<=24, the
  first five jets give the complete fixed rank-35 negative-depth
  chamber atlas, separate the local wall valuations contributed by
  q200, c300, and f400, and prove divisibility by a corrected
  B^Smith_M.  The global inherited factor remains q200^5*c300; no
  global f400 divisibility is asserted.  Only (12,5)
  remains a matrix-level sporadic after this separation.  A separate
  two-full-chart sidecar verifies the extra common seam E_M only for
  6<=M<=20.  No arbitrary-width,
  positive-core, or fixed-prime rank conclusion is claimed.
source: codex-gmc-local-smith-jet-barcode-2026-07-29
audit: >
  Author and independent exact replays completed in normal and optimized
  modes; all LF-normalized transcripts byte-match the stored output.
  The independent audit rederived the DVR lemma and degree/count laws,
  checked the local-versus-global factor scope, and independently
  eliminated determinant series modulo three primes at (M,r)=(21,17)
  and (12,5), obtaining the claimed orders 20 and 25 with nonzero
  leading coefficients.
depends_on:
  - THM-2949-fixed-rank-thirty-five-cofactor-newton-atlas
related:
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
  - THM-2944-width-nine-ten-two-chart-macaulay-resultant-closure
  - THM-2946-full-macaulay-maximal-minor-gcd-and-chart-free-resultant
  - THM-2957-first-gap-width-fifteen-twenty-modular-depth-ladder
  - THM-2959-first-gap-width-twenty-one-twenty-four-modular-depth-continuation
script: 04-computation/gmc_local_smith_jet_fitting_barcode_thm2960.py
output: 05-knowledge/results/gmc_local_smith_jet_fitting_barcode_thm2960.out
script_sha256: 6ff7360458a406994a4fefa787c2f9b231f1d577be3b9ca5fd48da90ed26d365
output_sha256: 89a6cb15fb0e24b696bb0fe4fb52453c40455b9a395fc62fe19a32ac99addbe1
thm2949_dependency_sha256: 9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54
hash_basis: LF-normalized bytes
---

# THM-2960 -- local Smith jets and the negative-depth chamber atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

For \(6\leq M\leq24\), let \(A_M(n)\) be the fixed
\(35\)-by-\(35\) polynomial matrix whose determinant \(P_M(n)\) is
THM-2949's fixed rank-thirty-five cofactor for the first-gap support

```text
(0,1,2,M).                                               (1)
```

The determinant orders at

```text
n=-1/2,-1,-2,...,-M                                     (2)
```

are not accidental factors and are not determined by ordinary rank.
They are lengths of local Smith bars.  The first five matrix jets
recover every bar for the fixed cofactor.  Pure coefficient
resonances form a degree-two/three/four ladder; after separating that
ladder, exactly one matrix-level exception to the four ratio chambers
remains.

## 1. The Smith-jet/Fitting barcode lemma

Let \(K\) be a field, \(R=K[[t]]\), and let

```text
A(t):R^v -> R^u
```

have generic rank \(s\).  Over the DVR \(R\), Smith reduction gives

```text
U A V = diag(t^alpha_1,...,t^alpha_s) direct_sum 0,      (3)
```

where \(U,V\) are invertible and every \(\alpha_i\geq0\).  Let \(T_k\)
be the lower block-Toeplitz matrix formed from the first \(k\)
coefficients of \(A(t)\).  It represents the induced map

```text
(K[t]/t^k)^v -> (K[t]/t^k)^u.                           (4)
```

Write

```text
d_k = dim_K ker(T_k) - (v-s)k,             d_0=0.       (5)
```

Then

```text
d_k = sum_i min(k,alpha_i),
e_k := d_k-d_(k-1) = #{i: alpha_i>=k},
#{i: alpha_i=k} = e_k-e_(k+1)
                 = 2d_k-d_(k-1)-d_(k+1).               (6)
```

In particular, as soon as \(e_{k+1}=0\), the finite sequence
\(d_1,\ldots,d_{k+1}\) recovers the whole Smith multiset.  The
valuation of the maximal determinantal divisor is

```text
v_t(gcd of all nonzero s-minors)=sum_i alpha_i.         (7)
```

For a square generically invertible matrix, `(7)` is the valuation of
the zeroth Fitting divisor of the cokernel and

```text
v_t(det A)=sum_i alpha_i=lim_k d_k.                     (8)
```

*Proof.*  Invertible local row and column operations induce
automorphisms modulo \(t^k\), so they preserve the adjusted kernel
dimension `(5)`.  Multiplication by \(t^\alpha\) on
\(K[t]/t^k\) has kernel dimension \(\min(k,\alpha)\).
This proves `(6)`.  Every nonzero maximal minor of the Smith form is a
unit times \(t^{\sum\alpha_i}\), and some such minor is exactly that
unit multiple, proving `(7)`.  Equation `(8)` is the square case.

Thus the graph \(k\mapsto d_k\) is a barcode: a bar of length
\(\alpha_i\) contributes one unit to each increment \(e_1,\ldots,
e_{\alpha_i}\).  This is the rigorous local ``holotopy'' bridge used
below.  It is characteristic-zero linear algebra here; it makes no
claim about reducing the full Macaulay matrix at one fixed prime.

## 2. The exact fixed-cofactor chamber table

Put \(t=n-c\) at a rational centre \(c\), expand

```text
A_M(c+t)=A_0+A_1t+A_2t^2+... ,
```

and form \(T_1,\ldots,T_5\).  Write a Smith multiset as
`1^a,2^b,3^c,4^d` when exponent \(j\) occurs with the displayed
multiplicity.  Zero exponents are omitted.

Apart from the resonances in Section 3, the complete exact table is:

| centre | condition | Smith multiset | \((d_1,\ldots,d_5)\) | order |
|:---|:---|:---|:---|---:|
| \(-1/2\) | all \(M\) | `1^5` | `(5,5,5,5,5)` | 5 |
| \(-1\) | all \(M\) | `2^1,4^1` | `(2,4,5,6,6)` | 6 |
| \(-2\) | all \(M\) | `1^1,2^5,3^2,4^1` | `(9,17,20,21,21)` | 21 |
| \(-r\) | \(3\le r<M,\ 3r\le M\) | `1^5,2^6,3^3` | `(14,23,26,26,26)` | 26 |
| \(-r\) | \(3r>M,\ 2r\le M\) | `1^5,2^8,3^1` | `(14,23,24,24,24)` | 24 |
| \(-r\) | \(2r>M,\ 3r\le2M\) | `1^1,2^5,3^3` | `(9,17,20,20,20)` | 20 |
| \(-r\) | \(3r>2M,\ r<M\) | `1^1,2^3,3^4` | `(8,15,19,19,19)` | 19 |
| \(-M\) | terminal centre | `1^9,2^5` | `(14,19,19,19,19)` | 19 |

The weak inequalities are load-bearing:

```text
3r=M      lies in the order-26 chamber,
2r=M      lies in the order-24 chamber,
3r=2M     lies in the order-20 chamber.                (9)
```

The high and terminal centres have the same total length nineteen but
different barcodes.  Therefore determinant order alone loses local
module structure which the jets retain.

## 3. Pure-coefficient resonances and sharp small hostiles

Use the THM-2943/2949 coefficient notation

```text
q200=[x0^2]Q,  c300=[x0^3]C,  f400=[x0^4]F.           (10)
```

On the complete finite bank \(6\le M\le24\),

```text
ord_(-r) q200 = 1  iff  M=1 mod 6 and 3r=2M+1,
ord_(-r) c300 = 1  iff  M=1 mod 4 and 4r=3M+1,
ord_(-r) f400 = 1  iff  M=1 mod 10 and 5r=4M+1.        (11)
```

The orders are zero otherwise, and none of the three coefficients
vanishes at \(-1/2\).  For the quartic line, putting \(d=M-r\)
rewrites the condition as

```text
M=5d+1,                    r=4d+1,                    d even.  (11a)
```

In the same \(d=M-r\) notation the full ladder is

```text
q: M=3d+1, r=2d+1, d even;
c: M=4d+1, r=3d+1;
f: M=5d+1, r=4d+1, d even.                            (11b)
```

The parity conditions are respectively equivalent to
\(M\equiv1\pmod6\) and \(M\equiv1\pmod {10}\).  These are exact laws
on the declared finite bank; no all-width extrapolation is made.

The smallest hostile to a bare four-chamber floor law is

```text
(M,r)=(7,5):
ord q200=1, ord c300=ord f400=0,
Smith multiset 1^2,2^4,3^2,4^2,
d=(10,18,22,24,24), order=24.                         (12)
```

The high-chamber floor is nineteen; `q200^5` supplies exactly the five
missing units.  The smallest cubic hostile is

```text
(M,r)=(9,7):
ord q200=ord f400=0, ord c300=1,
Smith multiset 1^2,2^3,3^4,
d=(9,16,20,20,20), order=20.                          (13)
```

Here `c300` supplies the one missing unit.  The quartic controls are

```text
(M,r)=(11,9):
  q200,c300 units; ord f400=1;
  1^2,2^3,3^4; order=20;

(M,r)=(21,17):
  q200,c300 units; ord f400=1;
  1^2,2^3,3^4; order=20.                              (14)
```

At \(M=21\), the new length-one bar is exactly THM-2959's surplus
factor \(n+17\).  These examples show why separating only `q200` and
`c300` falsely makes the quartic rung look matrix-sporadic.

After all three pure coefficients are separated, the unique remaining
exception is

```text
(M,r)=(12,5):
  q200,c300,f400 units;
  1^5,2^7,3^2; order=25.                              (14a)
```

Here one length-two bar becomes length three, so this is a genuine
matrix-jet resonance rather than a pure-coefficient root.

## 4. Corrected Smith factor and exact separation

For \(1\le r\le M\), define

```text
beta_M(r)=
  6,   r=1;
  21,  r=2;
  26,  3<=r and 3r<=M;
  24,  3r>M and 2r<=M;
  20,  2r>M and 3r<=2M;
  19,  3r>2M,                                          (15)

beta^mat_M(r)=beta_M(r)+1_{(M,r)=(12,5)},

epsilon^F_M(r)=ord_(-r) f400
 =1_{M=1 mod 10 and 5r=4M+1},

beta^Smith_M(r)=beta^mat_M(r)+epsilon^F_M(r).          (16)
```

The complete exact order decomposition is

```text
ord_(-r) P_M
 = beta^mat_M(r)
   +5 ord_(-r) q200
   +  ord_(-r) c300
   +  ord_(-r) f400,                                  (17)

ord_(-1/2) P_M=5.                                     (18)
```

Therefore the local Smith ranks prove directly, without interpolating
or factoring \(P_M\), that

```text
B^Smith_M(n)
 =(2n+1)^5 prod_(r=1)^M (n+r)^beta^Smith_M(r)         (19)
```

divides \(P_M\) in \(\mathbb Q[n]\).  Its degree is

```text
19M+2 floor(M/3)+4 floor(M/2)+floor(2M/3)-20
  +1_{M in {11,12,21}}.                               (20)
```

The quartic-root units remain part of the actual Smith divisor `(19)`;
the local separation in `(17)` explains their source without asserting
that the full polynomial `f400` divides \(P_M\).  The global inherited
pure factor is still `q200^5*c300`.  For \(15\le M\le20\), the
correction vanishes and `(19)` is exactly
the factor \(B_M\) used in THM-2957.  Thus THM-2957's negative-root
division has a local-module proof: its `26/24/20/19` staircase is a
Smith-length spectrum.  THM-2957's modular core gates remain logically
separate and are still needed for nonnegative integral-depth
nonvanishing.  At \(M=21\), `(19)` has degree \(448\): THM-2959's
floor-law degree \(447\), plus the quartic bar \(n+17\).

## 5. The two-full-determinant common-content sidecar

Let \(\Delta_{0,M}\) and \(\Delta_{1,M}\) be the original and
stable-mutated \(36\)-by-\(36\) Macaulay charts from THM-2943, specialized
to `(1)`, and let

```text
G_M=gcd(Delta_(0,M),Delta_(1,M)) in Q[n].              (21)
```

For a rational centre \(c\),

```text
ord_c G_M=min(ord_c Delta_(0,M),ord_c Delta_(1,M)).    (22)
```

The companion computes the first five jets of both full charts and,
only where an individual chart has not stabilized, extends through
the seventh jet.  At every centre in `(2)`, both chart orders are at
least the claimed value and at least one chart stabilizes at exactly
that value.  Hence `(22)` is certified without determinant
interpolation.

Define

```text
E_M(n)
 =(2n+1)
   prod_(2<=r<=floor(M/2)) (n+r)^3
   prod_(floor(M/2)<r<M)  (n+r)^4
   (n+M)^2.                                             (23)
```

On the finite bank \(6\le M\le20\), the exact local identities are

```text
ord_(-1/2) G_M=6=5+1,

ord_(-r) G_M
 =5 ord_(-r)q200 + ord_(-r)c300
  +beta^Smith_M(r)+ord_(-r)E_M.                        (24)
```

Equivalently, THM-2943's universal `q200^5*c300` factor and `(24)`
prove

```text
B^Smith_M E_M divides G_M/(q200^5 c300) in Q[n],       (25)
```

and the quotient in `(25)` is a local unit at every centre in `(2)`.
The extra \(E_M\) is therefore a full-chart seam, not part of the
fixed \(35\)-minor factor.

For widths nine and ten, PROVED THM-2944 globally identifies the
two-chart common gcd and proves the primitive chart cofactors coprime.
Equation `(24)` is the local Smith-bar refinement of that mechanism.
For widths outside those proved global atlases, `(24)` asserts only
the displayed finite wall-grid valuations.  It does **not** identify
the gcd of all \(36\)-minors of the full \(46\)-row matrix, and it
does not assert global coprimality of two primitive charts.

## 6. Why ordinary rank and row gcd are insufficient

At the exact hostile centre

```text
M=20,                         n=-3,                    (26)
```

each of the quadratic, cubic, and quartic form banks contains a local
unit.  Their minimum coefficient valuations are

```text
(0,0,0).
```

The evaluated fixed-cofactor matrix has corank fourteen, but

```text
(d_1,d_2,d_3)=(14,23,26),          ord_(-3)P_20=26.    (27)
```

Ordinary rank sees only the number of positive bars.  The other twelve
units of determinant order live in their lengths.  Row-gcd or
rank-at-the-centre arguments therefore discard exactly the information
which creates the staircase.

## 7. Exact replay and boundary

Run

```text
python 04-computation/gmc_local_smith_jet_fitting_barcode_thm2960.py
python -O 04-computation/gmc_local_smith_jet_fitting_barcode_thm2960.py
```

The companion hash-pins THM-2949 and reconstructs every matrix from
the canonical factorial forms.  It checks:

1. the barcode recovery control;
2. all
   \(\sum_{M=6}^{24}(M+1)=304\) fixed-cofactor centres, using
   \(1{,}520\) exact block-Toeplitz ranks;
3. every chamber profile, equality boundary, pure \(q/c/f\) resonance,
   and the sole matrix sporadic;
4. the exact degree law `(20)` and order decomposition `(17)`;
5. both full charts at the separate \(210\)-centre bank
   \(6\le M\le20\), with adaptive fifth-to-seventh jet stabilization;
6. the common-content identity `(24)` and the row-rank hostile `(27)`.

The frozen record digests are

```text
fixed:
c969fe43bd92e1a3e545341d89e09ad0a721ea312182e70c4e5786b9e75c3396

two-chart bridge:
035a164a768deb837d2f483c707eb124519709ae45d5fdef0ec08e3547e41f2a
```

This is a finite characteristic-zero local theorem: the fixed atlas
ends at width \(24\), while the two-chart bridge ends at width \(20\).
It proves no all-width chamber or resonance law, no nonvanishing of
the remaining cofactor core, no coefficientwise positivity of a
renormalized resultant, and no full-matrix rank claim modulo a fixed
or selected prime.  In
particular, the retracted all-width next-prime shortcut is unrelated
to and cannot be inferred from this atlas.

## 8. Independent hostile audit

The immutable candidate at commit `1320b187b` was independently replayed in
normal and optimized modes; both outputs LF-byte-match the stored transcript,
and the script, output, and THM-2949 dependency hashes agree with frontmatter.
The referee separately rederived the DVR/Toeplitz lemma, the degree formula,
and the `304`/`1,520`/`210` count invoices.

As an independent path through the two exceptional walls, truncated
determinant-series elimination modulo

```text
1000003, 1000033, 1000037
```

gave order `20` at `(M,r)=(21,17)` and order `25` at `(12,5)`, with nonzero
leading coefficients for every prime.  Direct pure-coefficient factorization
found `f400=(n+17)*unit` at the first wall while `q200,c300` are units, and
found all three pure coefficients to be units at the second.  This confirms
both the quartic-bar classification and the genuinely matrix-level boundary.
The audit found no defect.

**QED.**
