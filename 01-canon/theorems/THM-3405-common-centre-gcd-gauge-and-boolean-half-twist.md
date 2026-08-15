---
id: THM-3405
title: "Common-centre gcd gauge and Boolean half-twist reduction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every sheet number
  q and every finite selected speed set U, a vanishing THM-3398 complete
  mode cochain has a common centre controlled by one scalar a modulo 2q.
  If d=gcd(U), v_i=u_i/d, and g=gcd(q,d), then the exact word is
  h_i=a v_i mod 2q.  Common cyclic sheet relabelling changes a by 2d, so the
  raw gauge is Z/(2g).  The intrinsic divisibility
  gcd(q,u_i)|h_i forces g|a and collapses every realizable cover to at most
  the two gauge classes a=0 or a=g modulo 2g.  Thus, up to common sheet
  relabelling, the whole zero-cochain cover search is the union of the two
  literal source-centre layers c_0=0 and c_1=g/(2qd).  This is not a
  classification of arbitrary physical common-time covers and gives no
  LRC(14) decrement by itself.
source: root-2608-crouzeix-puzzle-2026-08-15
audit: self-contained Bezout/gcd/mode-centre proof; exact mode, orbit, literal-cover, dilation, and q=8/16/23 hostile controls; independent half-grid, prescribed-family, quantifier, parity, and strict-endpoint audit
depends_on:
  - THM-3398-general-finite-mode-sheet-cover-cochain
related:
  - THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight
  - THM-3402-atomized-sheet-covers-and-constructive-cochain-locus
script: 04-computation/lrc_common_centre_gcd_boolean_twist_thm3405.py
output: 05-knowledge/results/lrc_common_centre_gcd_boolean_twist_thm3405.out
script_sha256: 02d789a3f209a3a418ce234f7d124c99976c933efbb8a01d78ab3738e737686c
output_sha256: 969e84720a05f48a6e9e03be971d9a63cc684e43f6eeae659a67c608bdd2e55b
semantic_sha256: 8cf22bc02893c04f54c4ac7d70121b26b94e71fc0a596604421a42421fb54ab6
hash_basis: LF-normalized bytes
---

# THM-3405 -- common-centre gcd gauge and Boolean half-twist reduction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Inheritance and exact connection

[THM-3398](THM-3398-general-finite-mode-sheet-cover-cochain.md) replaces an
arbitrary cyclic sheet cover by selected consecutive phase blocks carrying a
complete affine integral cochain.  Its zero-cochain slice says that the
selected interval centres coincide.  MISTAKE-384 records the corrected near
miss: equality of all centres leaves one common additive gauge and does not
force the physical source centre to be zero.

This theorem computes that residual gauge exactly.  The canonical hostile is
the `q=16` cover on speeds `(2,6,10,14)`: its cochain vanishes at centre
`1/32`, but the fixed-zero rank from
[THM-3401](THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight.md)
is five.  The least-used sidecar is the common speed gcd, which was invisible
after retaining only pairwise gap values.

| field | exact connection |
|---|---|
| source | THM-3398 selected-block modes with complete cochain identically zero |
| target | two literal cyclic sheet-cover layers |
| map | contract all centre congruences to the scalar `a=2qdc` |
| preserved | every selected block, owner, strict endpoint, common sheet relabelling, and cover predicate |
| destroyed | the absolute sheet origin inside one gauge orbit |
| required sidecar | `d=gcd(U)` and `g=gcd(q,d)` |
| cheapest decisive tests | the `q=16` half-twist and the fixed-zero `q=15` cover |

## 2. Setup and the scalar centre word

Fix `q>=2` and `r>=1`.  Choose distinct positive speeds

```text
U={u_1,...,u_r},       d=gcd(u_1,...,u_r),
v_i=u_i/d,             gcd(v_1,...,v_r)=1.            (1)
```

For each owner choose one THM-3398 mode, with sheet block `B_i`, centre
residue `h_i mod 2q`, and positive width.  A centre lift has the form

```text
x_i=(n_i+h_i/(2q))/u_i,       n_i in Z.               (2)
```

The complete mode cochain is zero exactly when all `x_i` equal one rational
number `c`.  In that case

```text
2q u_i c == h_i  (mod 2q).                             (3)
```

Choose Bezout integers `b_i` with `sum_i b_i v_i=1`.  Equation `(3)` makes

```text
a=2qdc=sum_i b_i(2q u_i c) in Z                       (4)
```

and gives the scalar centre word

```text
h_i == a v_i  (mod 2q)       for every i.             (5)
```

Conversely, any integer `a` satisfying `(5)` gives the common lift

```text
c=a/(2qd),
n_i=(a v_i-h_i)/(2q).                                  (6)
```

Thus `(5)` is necessary and sufficient.  It is also equivalent to the
pairwise zero fibres from THM-3398.  Indeed those fibres reduce to

```text
h_i v_j-h_j v_i == 0  (mod 2q gcd(v_i,v_j)).           (7)
```

Equation `(5)` implies `(7)`.  Conversely, define
`a=sum_i b_i h_i mod 2q`; multiplying `(7)` by the Bezout coefficients and
summing gives `a v_j-h_j == 0 mod 2q`.

## 3. The raw gcd gauge

Relabel every sheet by the same cyclic translation `ell -> ell-k`.  In the
THM-3398 mode coordinates this sends

```text
h_i -> h_i+2u_i k,
a   -> a+2dk                 (mod 2q).                 (8)
```

The cover predicate is unchanged.  Put

```text
g=gcd(q,d).                                             (9)
```

The subgroup generated by `2d` in `Z/(2q)` has index `2g`; hence the raw
common-centre gauge is

```text
Z/(2q) / <2d>  ~=  Z/(2g),                             (10)
```

with orbit label `a mod 2g`.  This is the additive coordinate omitted when
fixed source centre zero was identified with the full zero-cochain slice.

## 4. Mode divisibility collapses the gauge to two points

For one owner let

```text
e_i=gcd(q,u_i).
```

THM-3398 writes its mode centre as

```text
h_i=-e_i(2r_i+s_i-1)  (mod 2q),                        (11)
```

so `e_i|h_i`.  Since `g|e_i` and `e_i|2q`, equation `(5)` gives
`g|a v_i` for every `i`.  Bezout in `(1)` now gives

```text
g|a.                                                   (12)
```

There are therefore only two possible raw gauge classes:

```text
a == 0 or g  (mod 2g).                                 (13)
```

After a common sheet relabelling, their canonical centre representatives are

```text
c_0=0,
c_1=g/(2qd).                                           (14)
```

The second class has exact order two.  The conclusion is “at most two”: for
a specified owner and mode bank, either class may be empty.

## 5. Exact two-layer cover test

For `epsilon in {0,1}`, define the centred mode bank of owner `u_i` by

```text
M_i^epsilon={B_i:
  B_i is a THM-3398 selected block and
  h_i == epsilon g v_i (mod 2q)}.                      (15)
```

For a prescribed active family, the following are equivalent.

1. One selected mode for every owner in `U` covers `Z/qZ` and the complete
   cochain is zero.
2. For one `epsilon in {0,1}`, one can choose a block from every used bank
   `M_i^epsilon` so that the chosen blocks cover `Z/qZ`.
3. For one canonical centre `c_epsilon` in `(14)`,

   ```text
   D_(q,u)(c_epsilon) is nonempty for every u in U,
   union_(u in U) D_(q,u)(c_epsilon)=Z/qZ,             (16)
   D_(q,u)(t)={ell: ||u(t+ell/q)||<1/14}.
   ```

There is a second, equally useful ambient-pool form, but its gcd must be
recomputed.  Some zero-cochain subfamily of an ambient pool `U` covers if and
only if there is a nonempty `I subset U` for which `(16)` holds after replacing
`U,d,g` by `I,d_I,g_I`.  One may discard empty owners, but one may not compute
the two centres from the gcd of the full ambient pool and then discard them.

Only `(2)<->(3)` needs comment.  Put `e=gcd(q,u)`, `m=q/e`, and
`b=u/e mod m`.  At a canonical centre, `H=2quc` is divisible by `e`:
writing `q=gq_0,d=gd_0` with `gcd(q_0,d_0)=1` gives
`e=g gcd(q_0,v)`, which divides `gv`.  Moreover,

```text
u(c+ell/q) == (H/e+2b ell)/(2m)  (mod 1).             (17)
```

Reflection about zero preserves this half-grid.  Hence every nonempty danger
set in `(17)` is one symmetric consecutive phase block.  If its size is `s`
and its first phase is `r`, symmetry gives

```text
H/e == -(2r+s-1)  (mod 2m).                           (18)
```

The strict one-seventh inequality gives `s<=ceil(m/7)`, so `(18)` is exactly
the THM-3398 mode formula `(11)`.  The danger set is the unique maximal block
in the bank `(15)`.  Conversely, every block in that bank is dangerous at
its interval centre.  This proves the literal test without treating an
arbitrary physical time as a mode centre.

## 6. Dilation

For every positive integer `lambda`, THM-3398's dilation is

```text
(q,u_i,c) -> (lambda q,lambda u_i,c/lambda).           (19)
```

It sends `d->lambda d`, `g->lambda g`, and

```text
a=2qdc -> lambda a.                                   (20)
```

Thus it preserves `epsilon` in `(13)`.  Each old block pulls back under
`Z/(lambda q)->Z/q`, so both Boolean layers and the cover predicate are
stable under dilation.

## 7. Sharp controls and loss boundary

- At `q=15`, speeds `(1,2,3,4,5,7)` cover at `c_0=0`.  This is the fixed-zero
  layer used by THM-3401.
- At `q=8`, speeds `(1,3,5,7)` partition the sheets at `c_1=1/16`.
- Dilating by two gives the MISTAKE-384 hostile: at `q=16`, speeds
  `(2,6,10,14)` partition the sheets at `c_1=1/32`, while the fixed-zero
  minimum has rank five.
- At `q=23`, speeds `(1,4,5,7,9,11)` cover at `c_1=1/46`.  A common sheet
  relabelling identifies this with the independently found centre `1/2`.
- The first prescribed-family scope hostile is `q=8`,
  `U=(1,2,3,5,7)`, `c_1=1/16`.  The literal union is full, but speed two has
  an empty danger set.  Thus the ambient subfamily `(1,3,5,7)` is a valid
  zero-cochain cover while all five owners are not a valid active family.
- The gcd-recomputation hostile is `q=9`.  The active family
  `(2,10,12,14)` has `d=2` and covers at its half centre `1/36`.  Adding the
  unused speed one changes the ambient gcd to one; neither ambient canonical
  centre `0,1/18` covers.  Thus an ambient pool requires an explicit subfamily
  quantifier before applying `(14)`.

The theorem does not say that every arbitrary physical cover time belongs to
one of `(14)`.  A nonzero affine cochain retains relative centre gaps, and
even at a common physical time the chosen mode intervals need not be centred
there.  The Boolean reduction applies exactly to vanishing complete mode
cochains.  It changes no LRC wall, row, or ledger count.

## 8. Exact companion

Run

```bash
python3 04-computation/lrc_common_centre_gcd_boolean_twist_thm3405.py
python3 -O 04-computation/lrc_common_centre_gcd_boolean_twist_thm3405.py
```

The standard-library-only companion has no float or assertion-dependent
truth gate.  It checks `99,166` selected-mode tuples, `102,336` translation
orbits, `3,855` bounded prescribed-family cover instances, `21` dilations,
the four positive controls, and both quantifier/gcd scope hostiles.  Normal
and optimized runs are byte-identical.
