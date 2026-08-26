---
id: THM-4204
title: "Rule 30 de Bruijn resets, unique-fibre sequence, and dyadic inverse ancestry"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.  The inverse-boundary de Bruijn
  monoid has 41 transformations (40 from nonempty words) and sharp reset
  threshold four.  Its
  unique-predecessor output language has a 13-state minimal DFA and an exact
  order-five complementary counting law, so asymptotically almost every
  cyclic output has one predecessor.  Under uniform finite-ring input, the
  probability of a multiple-predecessor output is exactly V_N/2^N and the
  one-step Shannon entropy deficit is Theta((rho/2)^N).  On every ring retaining a two-cell zero
  gap, each isolated physical row from time two onward has one cyclic
  predecessor before a final sparse/dense fork;
  on rings of size four times a dyadic time, the next ancestor count alternates
  exactly 3,1.  The isolated-seed temporal center prefix reaches a constant
  product at length 15, making this unpointed spatial quotient blind to every
  later center bit.  No Rule 30 prize is solved.
source: rule30-dyadic-subsequences-20260826
audit: >
  PASS.  The local rule independently regenerates both de Bruijn matrices.
  Exhaustive transformation closure, Moore minimization, an exact matrix-
  polynomial recurrence certificate, two direct center implementations,
  direct cyclic enumeration through size nine, physical-row controls through
  size eleven, all sixteen background-defect cases, and exhaustive dyadic
  depth controls at (h,N)=(2,8),(4,16) agree.  Ordinary and optimized runs
  reproduce the stored output byte for byte.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
related:
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3488-rule30-inward-slack-monicity-and-parity-cartier-ramification
  - THM-3489-rule30-packed-restart-and-pointed-pascal-face
  - THM-3500-rule30-dyadic-section-cut-defect-and-cross-depth-valuation-carrier
  - THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge
  - THM-4048-rule30-periodicity-balance-and-model-firewalls
  - THM-4206-rule30-characteristic-address-contrast-deck-entropy-decomposition
script: 04-computation/rule30_debruijn_reset_thm4204.py
output: 05-knowledge/results/rule30_debruijn_reset_thm4204.out
script_sha256: 75610536cd469c8aa68160f172a52bf410138c9494e3cd23af960bbef8f965c2
output_sha256: 3da2f073c80bafa84b46198f4804d45950022e9e16626084333b71a42d743c30
hash_basis: raw LF bytes
---

# THM-4204 -- Rule 30 de Bruijn resets, unique-fibre sequence, and dyadic inverse ancestry

**PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.**

This theorem extracts a positive spatial inverse theorem and a sharp no-go
from THM-3458's least-used representation.  Rule 30 has a shortest
synchronizing output word of length four; that reset makes almost every
periodic-ring output uniquely invertible and makes every sufficiently deep
physical isolated-seed row backward-rigid.  The same reset also shows why the
unpointed de Bruijn product cannot see the named temporal center: its product
on the center word becomes constant after fifteen symbols.

## 1. Inheritance, portfolio, and types

The inheritance pass is:

1. closest proved mechanism: THM-3458's exact de Bruijn trace compiler;
2. canonical hostile: THM-3463's large finite 2-kernel does not prove an
   infinite kernel, and a spatial inverse count is not a temporal observer;
3. corrected near miss: THM-3471/3488 show that scalarizing before retaining
   a marked transverse channel can annihilate the relevant current; and
4. least-used sidecar: the inverse boundary-pair state before taking a trace.

The live concept board has five objects:

1. the temporal center prefix `c_0...c_(N-1)`;
2. the two de Bruijn boundary transformations;
3. their rank-one reset ideal;
4. physical isolated-seed rows on large cyclic rings; and
5. the sparse/dense two-background ancestry fork.

The **Anchor** is dyadic center information, the **Niche** is the inverse
de Bruijn monoid, and the **Wildcard** is its unique-fibre counting
sequence.  The native-operation warning is load-bearing throughout.  A word
of center values across time will sometimes be reinterpreted below as one
word across space; the theorem marks that change of type explicitly.

## 2. The 41-element inverse-boundary monoid

Retain THM-3458's ordered states `00,01,10,11` and matrices

```text
A_0 = [[1,0,0,0],             A_1 = [[0,1,0,0],
       [0,0,0,0],                    [0,0,1,1],
       [0,1,0,0],                    [1,0,0,0],
       [0,0,1,1]],                   [0,0,0,0]].       (1)
```

For an output word `w=w_0...w_(N-1)`, put

```text
M_w=A_(w_0)...A_(w_(N-1)).                            (2)
```

Every column of each generator has exactly one `1`.  This is left
permutivity read backward: an output bit and a target pair `(b,c)` determine
the unique source pair `(a,b)`.  Thus the columns encode transformations

```text
g_0=(0,2,3,3),             g_1=(2,0,1,1),             (3)
```

and matrix multiplication is transformation composition.  A direct closure
of (3), including the empty product identity, contains exactly

```text
41 monoid transformations, with image ranks
rank 1 / 2 / 3 / 4:        4 / 30 / 6 / 1.            (4)
```

This is a finite exact census: closure is checked inside the `4^4=256`
possible transformations.  Every nonempty product has rank at most three,
so the identity cannot recur: the strict nonempty-word semigroup has 40
elements.  The four rank-one elements are precisely the constant
transformations.

For `i in {0,1,2,3}`, let `E_i` be the matrix with every column equal to the
`i`th basis vector.  Column sums in (1) give

```text
E_i A_0=E_i A_1=E_i.                                  (5)
```

The exact shortest resets are

```text
M_0010=E_3,                  M_1010=E_1.              (6)
```

There is no reset word of length at most three, and these are the only two
at length four.  The complete rank lists at lengths one, two, and three are

```text
3,3;
2,3,2,3;
2,2,2,3,2,2,2,3,                                  (7)
```

in lexicographic word order.  Equations (6)--(7) prove that the reset
threshold four is sharp; no asymptotic fit is involved.

### 2.1 A universal unique-predecessor cylinder

THM-3458 proves that the number of periodic predecessors of a cyclic output
word `w` is `tr(M_w)`.  Since every `E_i` has trace one, any **linear** word
containing `0010` or `1010` has one periodic predecessor.

More invariantly, every cyclic word of length `N>=3` containing the cyclic
factor `010` has exactly one periodic Rule 30 predecessor.  For `N>=4`, rotate
the word to begin `b010`, use (6), (5), and cyclic invariance of trace.  At
`N=3`, direct multiplication gives `tr(M_010)=1`; its rotations have the same
trace.

Rank itself is not cyclically invariant.  The smallest hostile is

```text
rank(M_0101)=2,        rank(M_1010)=1,
tr(M_0101)=tr(M_1010)=1.                              (8)
```

Thus cyclic rotation may be used for the trace conclusion, not for an
unqualified claim about the displayed product's rank.

## 3. Exact closed form for the unique-fibre sequence

Let `U_N` be the number of labelled binary cyclic output words of length `N`
having exactly one periodic predecessor, and put

```text
V_N=2^N-U_N.                                          (9)
```

For the formal empty word set `V_0=1`.  Moore minimization of the 41-state
transformation automaton, accepting exactly transformations with one fixed
point, gives the following 13-state DFA.  State zero is initial; `u=1` means
unique predecessor.

```text
state   u   on 0   on 1
  0     0     1      2
  1     0     1      3
  2     0     4      5
  3     1     6      7
  4     1     4      8
  5     0     9      0
  6     1     6      6
  7     0     1      1
  8     0     6     10
  9     0     9     11
 10     0     4      4
 11     0     6     12
 12     1     9      9.                               (10)
```

Let `T_(ij)` count the symbols in `{0,1}` sending state `i` to state `j`,
and let `b_i=1-u_i`.  Exact integer multiplication gives the compact
annihilating certificate

```text
(T^5-3T^2-2T-2I)b=0.                                 (11)
```

Since `V_N=e_0 T^N b`, (11) proves, for every `N>=5`,

```text
boxed:
V_N=3V_(N-3)+2V_(N-4)+2V_(N-5),                     (12)

(V_0,V_1,V_2,V_3,V_4)=(1,2,2,5,10).                 (13)
```

Equivalently,

```text
boxed:
sum_(N>=0) V_N z^N
  =(1+2z+2z^2+2z^3+2z^4)/(1-3z^3-2z^4-2z^5),        (14)

U_N=2^N-V_N.                                         (15)
```

If `r_1,r_2,r_3` are the roots of `x^3-x^2-2`, partial fractions give the
literal closed form

```text
V_N=r_1^N+r_2^N+r_3^N
    +(-2 if 3 divides N, else 1).                    (15a)
```

The first positive-length values of `U_N` are

```text
0,2,3,6,20,41,84,190,399,822,1716,... .              (16)
```

Factor the denominator of (14) as

```text
(1+z+z^2)(1-z-2z^3).                                 (17)
```

If `rho=1.695620769...` is the positive root of

```text
rho^3=rho^2+2,                                       (18)
```

then `1/rho` is the unique smallest-modulus pole in the `z`-plane;
equivalently, the other characteristic roots have modulus below `rho`.
Formula (15a) also shows that the dominant coefficient is exactly one.  Thus

```text
V_N=Theta(rho^N),
U_N/2^N=1-Theta((rho/2)^N).                           (19)
```

So almost every labelled periodic output has exactly one predecessor, with
an exact exponentially small exceptional proportion.  This does **not** make
the global ring map bijective: THM-3458's `0^N` and `1^N` fibres remain exact
zero/multiple-predecessor hostiles.

### 3.1 Exact finite-ring Haar entropy bridge

Let `X_N` be uniform on the `2^N` labelled cyclic input rows and put
`Y_N=F_N(X_N)`.  For an output word `w`, write

```text
k_w=|F_N^(-1)(w)|,        n_k=|{w:k_w=k}|.             (19a)
```

The transformation product in Section 2 has rank at most three for every
nonempty word, so `k_w<=3`.  Moreover

```text
n_1=U_N=2^N-V_N,             sum_w k_w=2^N.            (19b)
```

Subtracting the `k=1` contribution gives the exact identity

```text
sum_(k=2)^3 k n_k=V_N.                                  (19c)
```

Thus a uniform input lands in a multiple-predecessor output with probability
exactly

```text
P(k_(Y_N)>=2)=V_N/2^N.                                  (19d)
```

Because `F_N` is deterministic,

```text
N-H(Y_N)=H(X_N|Y_N)
 =2^(-N) sum_(k=2)^3 k n_k log_2(k).                    (19e)
```

Equations `(19c)--(19e)` imply

```text
V_N/2^N <= N-H(Y_N) <= log_2(3) V_N/2^N
                         =Theta((rho/2)^N),             (19f)
H_infinity(Y_N)>=N-log_2(3).                            (19g)
```

This is a lossless inverse-fibre-to-Haar bridge for one complete spatial row
on a finite ring.  It is complementary to THM-4206, whose entropy identity
uses forward characteristic-address pivots for a finite spacetime cell
family under infinite Bernoulli Haar input.  Neither statement supplies the
other's sidecar, and `(19f)` is not a named-seed temporal entropy theorem.

## 4. Physical isolated rows reset and become backward-rigid

Let `F_N` be Rule 30 on the labelled cyclic ring `C_N`.  Let `W_t^(N)` be the
ordinary isolated-seed time-`t` row embedded without wrap, so assume

```text
N>=2t+3.                                              (20)
```

Write its left-edge bits as

```text
ell_r(t)=a_t(-t+r),             ell_r(t)=0 for r<0.   (21)
```

The Rule 30 local rule gives the exact inward recurrence

```text
ell_r(t+1)
 =ell_(r-2)(t)+ell_(r-1)(t)
  +(1+ell_(r-1)(t))ell_r(t).                         (22)
```

Induction in `r` now gives

```text
ell_0=1                         (t>=0),
ell_1=1                         (t>=1),
ell_2=0, ell_3=t mod 2, ell_4=1 (t>=2),
ell_5=ell_6=t mod 2             (t>=3),
ell_7=0                         (t>=4).                (23)
```

In particular the stable left-edge prefix is

```text
110010       for even t>=4,
11011110     for odd  t>=5.                           (24)
```

Put

```text
P=A_0^2
 =[[1,0,0,0],
   [0,0,0,0],
   [0,0,0,0],
   [0,1,1,1]].                                       (25)
```

In fact `A_0^q=P` for every `q>=2`.  Direct multiplication of the two parity
prefixes gives the same reset:

```text
P M_110010=P M_11011110=E_0.                         (26)
```

The zero gap in (20) has length `q>=2`.  Let `S_t` be the support word read
from its left edge, and choose the canonical cyclic representative

```text
widetilde W_t=0^q S_t.                                (26a)
```

Equations (5), (24)--(26) prove

```text
M_(widetilde W_t)=E_0,              t>=4.             (27)
```

Another labelled rotation can have a nonconstant product, as (8) warns, but
its trace is still one.

The two small heads have

```text
tr(PM_11001)=tr(PM_1101111)=1                       (28)
```

at `t=2,3`.  Since the physical row `W_(t-1)^(N)` is one predecessor,
(27)--(28) prove

```text
boxed:
F_N^(-1)(W_t^(N))={W_(t-1)^(N)}       for every t>=2. (29)
```

Thus the actual physical history is backward-rigid at every depth after the
first row.  For a dyadic `t=2^m>=8`, the reset in (26) occurs strictly before
the spatial center.  Flipping that center, or changing any later spatial
suffix in the canonical linearization (26a), leaves its full unpointed product
`E_0`.  This is stronger than equal predecessor counts and is already a
precise center-blindness obstruction.

## 5. The two-background fork and cubic dyadic ancestry

Work over `F_2` and write an eventually constant configuration as

```text
x_j=b+p_j,                    b in F_2,
P(y)=sum_j p_j y^j,           P finite Laurent.       (30)
```

Define

```text
E(P)=sum_j p_j p_(j+1)y^j.                            (31)
```

Substitution in the Rule 30 polynomial, using `b^2=b`, gives a zero
background after one step and the exact defect

```text
boxed:
P'=yP+(1+b)(1+y^(-1))P+E(P).                         (32)
```

For the sparse and dense defects this yields

```text
F(delta_0)
 =y^(-1)+1+y
 =F(1+delta_(-1)+delta_0).                            (33)
```

The second input is the all-one background with adjacent zeros at `-1,0`.
The analogous `t=1` calculation is

```text
tr(PM_111)=2,                                        (34)
```

so (33) lists the complete two-element predecessor fibre of `W_1` whenever
`N>=5`.

Let

```text
D_N=00 1^(N-2)                                      (35)
```

cyclically denote that dense parent.  The sparse seed has exactly one cyclic
predecessor for every `N>=2`.  The dense branch has

```text
g(N)=#F_N^(-1)(D_N)=tr(P A_1^(N-2)).                 (36)
```

THM-3458's characteristic polynomial gives

```text
chi_(A_1)(lambda)=lambda(lambda-1)(lambda^2+lambda+1),
A_1^(k+3)=A_1^k                         for k>=1.    (37)
```

At exponents `k=0,1,2,3`, the traces in (36) are `2,1,0,2`.
Cayley--Hamilton then continues the period from `k>=1`, hence

```text
boxed:
g(N)=2 if N=2 mod 3,
     1 if N=0 mod 3,
     0 if N=1 mod 3,                                 (38)

sum_(N>=2)g(N)z^(N-2)=(2+z)/(1-z^3).                 (39)
```

Iterating (29) down to the fork (33), and then applying (36), proves

```text
|F_N^(-k)(W_t^(N))|=1,                 1<=k<t,
|F_N^(-t)(W_t^(N))|=2,
|F_N^(-(t+1))(W_t^(N))|=1+g(N).                      (40)
```

Now put `h=2^m`, `m>=1`, and `N=4h=2^(m+2)`.  Condition (20) holds.  Since
`2^(m+2) mod 3` alternates, (38)--(40) give the genuine dyadic law

```text
boxed:
|F_(4h)^(-(h+1))(W_h^(4h))|
  =3,      m odd,
  =1,      m even.                                   (41)
```

The first hostile pair is exact:

```text
(h,N)=(2,8): 3,                  (h,N)=(4,16): 1.    (42)
```

The alternation is the primitive cubic packet in `F_4`: Frobenius swaps its
two nontrivial `C_3` characters.  This is an infinite theorem, not a fitted
finite prefix.

## 6. The temporal center prefix saturates at length fifteen

Return now to the temporal word `c_t=a_t(0)`.  Direct CA evolution and
THM-3458's independent packed recurrence agree on

```text
c_0...c_14=110111001100010.                          (43)
```

The final four symbols are the reset `0010`.  Exact multiplication gives

```text
M_(c_0...c_14)=E_1.                                  (44)
```

The transformation ranks of the nonempty prefixes through length fifteen are

```text
3,3,2,2,2,2,2,2,2,2,2,2,2,2,1,                     (45)
```

so fifteen is the first reset length for this word.  By (5), **without
knowing any later center bit**, one has

```text
boxed:
M_(c_0...c_(N-1))=E_1                 for every N>=15. (46)
```

In particular,

```text
M_(c_0...c_(2^m-1))=E_1               for every m>=4, (47)
```

and each such temporal prefix, when retyped as a labelled cyclic **spatial**
output, has exactly one predecessor.  At length fifteen that predecessor is

```text
100010000111010 -> 110111001100010.                  (48)
```

This is exact dyadic-prefix saturation, not center-bit regularity.

### 6.1 Precise quotient obstruction

Let `p` be the fifteen-bit word in (43).  For every finite suffix `z`,

```text
M_(pz)=E_1.                                          (49)
```

At length `15+L`, all `2^L` continuations of `p` therefore occupy one
monoid fibre.  For an infinite suffix, apply (49) to each finite continuation
prefix.  Thus the infinite extensions `p0^infinity` and `p1^infinity` have
limiting densities zero and one but the same product state at every length
after fifteen; so does `p` followed by the indicator of powers of two, which
is not eventually periodic.  Hence no argument factoring only through the
sequence of unpointed products in (46) can decide eventual periodicity,
density, or single-seed computation of the actual center continuation.

## 7. Relation to the inherited Cartier, carry, and complexity lanes

The exact connection ledger is:

```text
source:     temporal center prefix; physical isolated rows; all cyclic outputs
target:     four-state spatial inverse transformations and their trace DFA
map:        reinterpret a binary word spatially, then w -> M_w
preserves:  one-step cyclic predecessor count; a synchronized boundary pair;
            the exact all-output unique-fibre language
destroys:   temporal/spatial type; marked center address after reset;
            forward seed phase, Cartier slack, ordinary carry, current location
sidecar:    a pointed pre-reset boundary path and a lawful cross-depth temporal map
hostiles:   (8), 0^N/1^N fibres, (42), and the full continuation fibre (49)
```

This separates four superficially similar composition laws.

- THM-3489's packed restart acts **forward in time** on a growing 2-adic row
  and retains the lift sheet.
- THM-3471's uniform macroblock compiler has a widening predecessor invoice;
  it is not iteration of the fixed 41-state inverse quotient.
- THM-3488's Cartier atlas retains two marked slack channels before
  scalarization, whereas (46) has already erased every later symbol.
- THM-3516/4048 recover the named center only with the moving address and
  ordinary/Haar carries.  Those are absent here.

THM-3463's finite 2-kernel lower bound is also unaffected: (10) reads the
**symbols of one spatial word** from left to right; it is not a DFAO reading
the binary digits of the time index `n` to output `c_n`.

THM-4206 is a proved but nonduplicating lane.  Both theorems start from left
permutivity: here it makes de Bruijn columns deterministic, while there it
makes each extreme characteristic input a triangular Haar pivot.  Their
targets differ.  THM-4206 retains a directed contrast deck for a finite family
of forward spacetime outputs and computes conditional entropy; THM-4204
counts inverse fibres of one cyclic spatial output and proves that its
unpointed product saturates.  Neither theorem depends on the other, and no
entropy/preimage or contrast-deck/de Bruijn intertwiner is claimed.

The primitive-cubic factor in (37) does suggest one disciplined next test.
THM-3463's coarse Rule-150 symbol `w^(-1)+1+w` also vanishes at a primitive
cube root.  The de Bruijn packet in (41) is an inverse-ring boundary phase;
the Rule-150 packet is a forward collision-current character.  They share a
kernel but not yet a carrier.  A bridge would require evaluating both on the
same `N=4h` physical geometry while retaining the marked pre-reset cut state.

## 8. Exact verification and prize scope

The cheapest complete universe is the closure of two transformations on four
states.  The companion:

1. regenerates (1) from the local rule;
2. closes and ranks all 41 transformations;
3. proves the sharp reset threshold and minimizes (10);
4. verifies the integer certificate (11) and (12) through length 64;
5. compares matrix traces with direct ring enumeration for every output
   through `N=9`;
6. compares direct and packed center prefixes through time 64;
7. checks the physical edge resets, every background-defect truth-table row,
   and direct physical fibres through `N=11`; and
8. exhaustively verifies the depth counts in (42).

Reproduce with

```bash
python3 04-computation/rule30_debruijn_reset_thm4204.py
python3 -O 04-computation/rule30_debruijn_reset_thm4204.py
```

Both runs reproduce the stored output byte for byte.  All three advertised
Rule 30 center prizes remain **OPEN**.  The proved positive results concern
cyclic spatial inversion and one coupled dyadic ring family; the proved
negative result says that the unpointed de Bruijn product is too lossy after
its reset.  It supplies neither center nonperiodicity, center balance, nor a
fixed-seed complexity lower bound.
