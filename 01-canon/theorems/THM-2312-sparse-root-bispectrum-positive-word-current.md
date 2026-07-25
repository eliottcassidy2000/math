---
id: THM-2312
title: "Sparse-root bispectrum positivity and exact word current"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For a
  nonnegative function on a cyclic group of order p with at most two
  occupied sites, the sum of all bispectra whose three characters are
  nonzero is at least (p-2)(p-4)/4 times the cube of the total mass. The
  constant is sharp and is positive exactly from p=5 onward. On every
  positive THM-2305 blocker-word stratum, p=13, one ordered nonzero
  character pair has integrated real bispectrum at least
  (322959/16)rho^3. This retains the selected owner, prescribed clock, and
  complete target word in one gauge-invariant cubic phase-closure current.
  It is invariant under cyclic root translation; in particular, on
  one-sheet fibres it is pure mass cubed and does not recover terminal
  component phase. It does not force the gcd(m,91)=1 shell colour, an
  ordinary pair coefficient, or LRC(14).
source: codex-2026-07-25-sparse-root-bispectrum
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
script: 04-computation/lrc14_sparse_root_bispectrum_current_thm2312.py
output: 05-knowledge/results/lrc14_sparse_root_bispectrum_current_thm2312.out
script_sha256: 6d80ff4460d720eff24e4339218c65bb3a21d2464dec6ffe73d5ea8bbadd8a4f
output_sha256: 47e9f4e5e9804bd8f7bce83c853932f15aba3a697ad5be67a93069d7f4d901c7
hash_basis: working-tree bytes (LF)
---

# THM-2312 -- a whole-character bispectrum survives

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2302 identifies the first potentially sufficient rooted object as
third-order: a phase-sensitive coefficient incident to the same owner and
current-service word. THM-2299 shows why a single ordinary Fourier
coefficient is not enough. The missing positive operation is to sum a
complete character face before selecting one term.

For a two-sheet nonnegative root fibre, this whole-face sum has an exact
cubic formula. It cannot cancel. This produces a nonzero bispectrum on the
same exact THM-2305 word, with a quantitative lower bound.

## 1. The complete nonzero-character identity

Let `p>=2`, let `zeta` be a primitive `p`th root of unity, and let

```text
v=(v_r)_(r in Z/pZ),             v_r in R.
```

Use the Fourier convention

```text
M_k=sum_r v_r zeta^(-kr),        k in Z/pZ,          (1)
```

and put

```text
S_1=sum_r v_r,
S_2=sum_r v_r^2,
S_3=sum_r v_r^3.                                    (2)
```

Define the whole nonzero-root bispectrum

```text
B_p(v)
 =sum_(k,l nonzero; k+l nonzero)
    M_k M_l conjugate(M_(k+l)).                      (3)
```

All character indices in (3) are modulo `p`.

> **Whole-character identity.**
>
> ```text
> B_p(v)=p^2 S_3-3p S_1 S_2+2S_1^3.                (4)
> ```
>
> In particular, `B_p(v)` is real.

### Proof

First sum over every ordered pair `(k,l)`. Character orthogonality gives

```text
sum_(k,l) M_k M_l conjugate(M_(k+l))
 =sum_(r,s,t) v_r v_s v_t
    sum_k zeta^(k(t-r)) sum_l zeta^(l(t-s))
 =p^2 S_3.                                          (5)
```

The excluded hyperplanes are

```text
k=0,       l=0,       k+l=0.
```

Each complete hyperplane sum is `p S_1 S_2`. For example,

```text
sum_l M_0 M_l conjugate(M_l)
 =S_1 sum_l |M_l|^2
 =p S_1 S_2                                         (6)
```

by Parseval; the other two are identical because `v` is real. Every
pairwise intersection and the triple intersection is the origin term
`S_1^3`. Inclusion-exclusion therefore gives

```text
p^2S_3-3pS_1S_2+3S_1^3-S_1^3,
```

which is (4). QED.

No primality assumption is used in (4).

## 2. Sharp two-sheet positivity

Assume now that `v_r>=0`, that `v` is nonzero, and that at most two sites
are occupied. Write their masses as `a,b>=0`, not both zero. Equation (4)
becomes

```text
B_p(v)
 =(p-2)((p-1)(a^3+b^3)-3ab(a+b)).                   (7)
```

Equivalently, with

```text
V=a+b,               x=ab/V^2<=1/4,
```

we have

```text
B_p(v)/V^3=(p-2)((p-1)-3px).                        (8)
```

Hence, for every `p>=5`,

```text
B_p(v)>=(p-2)(p-4)V^3/4>0.                          (9)
```

The constant is sharp. Equality in (9) holds exactly when both sheets are
occupied with equal mass, because `x=1/4` exactly when `a=b>0`.

At thirteen,

```text
B_13(v)>=(99/4)(sum_r v_r)^3.                       (10)
```

The threshold and hypotheses are real:

```text
p=3, two equal positive sheets: B_3=-2a^3;
p=4, two equal positive sheets: B_4=0;
p=13, seven equal positive sheets: B_13=-42a^3;
p=13, signed two-sheet weights (a,-a): B_13=0.      (11)
```

Thus neither a smaller radix, arbitrary support, nor signed fibre weights
inherit the positive theorem.

## 3. Exact THM-2305 word current

Use the notation of THM-2305. Fix its selected source owner `j`, prescribed
clock, and one positive canonical terminal word

```text
Q=Q_(j,sigma),
sigma in {{a},{b},{a,b}}.                            (12)
```

Let

```text
v_r(y)=G_j((y+r)/13),
M_k(y)=sum_(r=0)^12 v_r(y)zeta^(-kr),               (13)
```

and put

```text
rho_Q=integral_Q P G_j(y)dy>0.                       (14)
```

THM-2305 proves that every active fibre in (13) has at most two nonzero
entries. To make the word restriction literal, replace `v_r(y)` by
`1_Q(y)v_r(y)`. Since `1_Q^3=1_Q`, define

```text
C_(k,l)(Q)
 =integral_Q M_k(y)M_l(y)conjugate(M_(k+l)(y))dy,   (15)
```

for `k,l,k+l` all nonzero modulo thirteen.

Every word in (12) is contained in at least one target danger comb, whose
measure is `1/7`. Therefore

```text
measure(Q)<=1/7.                                    (16)
```

Also

```text
sum_r v_r(y)=13 P G_j(y),
integral_Q sum_r v_r(y)dy=13 rho_Q.                 (17)
```

Hölder and (16) give

```text
integral_Q (sum_r v_r)^3
 >=(13rho_Q)^3/measure(Q)^2
 >=49(13rho_Q)^3.                                   (18)
```

Apply (10) pointwise and integrate. Equations (3), (15), and (18) yield

```text
sum_(k,l,k+l nonzero) C_(k,l)(Q)
 >=(10657647/4)rho_Q^3.                             (19)
```

There are

```text
12*11=132                                           (20)
```

ordered character pairs in (19). Since its sum is positive real, at least
one pair satisfies

```text
Re C_(k,l)(Q)>=(322959/16)rho_Q^3>0.                (21)
```

Thus one exact cubic phase-closure current is nonzero on the same source
owner, prescribed clock, and complete terminal blocker word.

For THM-2305's selected large word, its two mass floors give the explicit
consequences

```text
strict:
 Re C_(k,l)
 >=59330091441448204464161088965287
   /204762211570812550122621742116528000000;

repeated-first:
 Re C_(k,l)
 >=2493436238631076936119894860797
   /204762211570812550122621742116528000000.        (22)
```

The selected character pair may depend on the row and word. No common pair
for all profiles is asserted.

## 4. What phase the bispectrum retains

A cyclic translation of the root sheets by `h` changes

```text
M_k -> zeta^(-kh)M_k.
```

The three factors in (15) cancel this common character gauge:

```text
zeta^(-kh)zeta^(-lh)zeta^((k+l)h)=1.                (23)
```

The cyclic invariance has an exact continuous-gauge refinement. Recall the
periodic gauge

```text
N_k(y)=exp(-2*pi*i*k*y/13)M_k(y).
```

For representatives `1<=k,l<=12`, put

```text
c=(k+l) mod 13 in {1,...,12},
q=(k+l-c)/13 in {0,1}.
```

Then the balanced root degrees leave precisely the carry frequency:

```text
M_k M_l conjugate(M_(k+l))
 =exp(2*pi*i*q*y)N_k N_l conjugate(N_c).            (23a)
```

Thus `C_(k,l)(Q)` is the zero spatial coefficient of the raw cubic current,
equivalently the spatial coefficient at frequency `-q` of the
gauge-corrected cubic product. Exactly `66` allowed pairs have `q=0` and
`66` have `q=1`. The bispectrum can retain gauge-invariant cross-character
phase closure and this single radix-carry oscillation when two sheets are
present, while deliberately forgetting the absolute root address. It is
not a nonnegative quadratic energy and it is not an ordinary pair Fourier
coefficient.

The one-sheet boundary makes the loss exact. If the unique occupied root is
`r` with mass `A`, then

```text
M_k=A zeta^(-kr),
M_k M_l conjugate(M_(k+l))=A^3                     (23c)
```

for every allowed pair `(k,l)`. Hence the entire bispectrum collection is
blind even to the root address on a one-sheet fibre. In THM-2303's two
hostile multiplier-four carriers `F_0,F_1`, all active fibres are one-sheet,
the two component widths agree, and every integrated bispectrum in the two
states is identical. With its `epsilon=10^(-12)`, each is exactly

```text
measure(R)/13^3=52*epsilon/13^3
               =1/42250000000000.                  (23d)
```

Nevertheless the prescribed linear verdict is zero for `F_0` and nonzero
for `F_1`. Thus (21) is a genuine signed cubic current, but not a lawful
THM-2303 phase-tree edge.

This is forced by character charge, not an accident of the example. A
Fourier monomial

```text
product_i M_(a_i) product_j conjugate(M_(b_j))
```

has charge

```text
chi=sum_i a_i-sum_j b_j mod 13.
```

Its translates are multiplied by `zeta^(-h*chi)`. If a same-charge scalar
is real nonnegative on every translate and positive on one, then necessarily
`chi=0`; for nonzero charge the thirteen orbit values sum to zero. Uniform
scalar positivity is therefore confined to the invariant sector. Reaching
an ordinary nonzero-character current requires an external charged sidecar,
such as a marked root address or a unit-coloured shell edge, before the
charge-zero contraction.

This is the useful whole-face move. Individual bispectra may have either
sign or complex phase; their complete nonzero-character sum has the positive
formula (4), and only then is one nonzero term selected. The mechanism is
formally analogous to retaining a complete balanced face before reduction,
but it does not import any coefficient or valuation claim from THM-2022.

## 5. Connection and stopping boundary

The exact connection contract is

```text
source:
  THM-2305's positive canonical blocker-word stratum and its at-most-two
  nonnegative predecessor sheets;

target:
  a nonzero gauge-invariant cubic phase-closure current on that same owner,
  clock, and complete target word;

map:
  sum the complete nonzero-character bispectrum face, prove its cubic
  positivity, and select one ordered character pair;

preserved:
  source owner, prescribed time, exact pure/fork word, nonzero root
  characters, signed complex multiplication, and a quantitative margin;

destroyed:
  absolute root address, the identity of a uniform character pair, linear
  Fourier amplitude, terminal-component pair phase, shell multiplier,
  nonzero character charge, and relation-lattice ancestry;

needed sidecar:
  land (21) on a THM-2302 shell edge with gcd(m,91)=1, or couple it to
  THM-2306's exact first-collision owner-multiple atom without erasing the
  cubic phase.                                                        (24)
```

THM-2299 remains a valid boundary: an anchored rank-six carrier can have
zero ordinary pair coefficient even though the cubic whole-character
aggregate here is positive. The present theorem bypasses that cancellation;
it does not reconstruct the missing pair phase or prove a lawful phase-tree
edge.

No scalar profile is excluded, and LRC(14) remains open.

## 6. Exact companion

The companion works in exact cyclotomic group algebra over `Fraction`. It
checks (4) directly in `Q[zeta_p]` for `377` one- and two-sheet cases at
`p=5,7,13`, the sharp equality classification, all constants in
(10), (19)--(22), and every hostile control in (11). Every load-bearing
test raises explicitly under both normal and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_sparse_root_bispectrum_current_thm2312.py
python3 -O 04-computation/lrc14_sparse_root_bispectrum_current_thm2312.py
```

The two transcripts must match

```text
05-knowledge/results/lrc14_sparse_root_bispectrum_current_thm2312.out
```

byte-for-byte after LF normalization.
