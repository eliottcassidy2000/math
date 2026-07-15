---
id: THM-804
title: Three-replacement oriented sheet-deck descent in the full-residue chart
status: CLAIMED / PROOF AUDIT IN PROGRESS — exact half-open capacity argument derived; theorem write-up and modular replay pending
source: codex-2026-07-15-S10 (recursive Hamming continuation)
depends_on: [THM-795, THM-800]
related: [THM-769, THM-770, HYP-6775, HYP-6800, HYP-6820]
verification: forthcoming exact modular replay
---

# THM-804 — Three-replacement oriented sheet-deck descent

This file reserves THM-804 while the proof and exact replay are audited.  The
claim is deliberately a common-scale descent, not scale-one Hamming-three
rigidity.

Let `13` not divide `c`, let `r_1,r_2,r_3` be distinct nonzero residues, and
replace `c r_i` in `c[12]` by positive `w_i` satisfying

```text
w_i = c r_i  (mod 13).
```

Put `D_i=c/gcd(c,w_i)`.  The claimed theorem is:

```text
M((c[12]\{c r_1,c r_2,c r_3}) union {w_1,w_2,w_3})=1/13
  implies D_1=D_2=D_3=1.
```

Equivalently, all three replacements share the AP scale and the packet
descends to a scale-one triple lift.  The scale-one triple-lift classification
remains a separate open problem.

## Derived proof spine awaiting final audit

At every lifted missing-owner splice, use the left core-safe germ.  A
replacement of deck order `D` can own only the phase interval

```text
(-1/13,1/13],
```

so its sheet capacity is `ceil(2D/13)/D`.  The half-open endpoint is essential.
For `D>=3` this capacity is at most `1/3`, with equality only at `D=3`; for
`D=2` it is `1/2`.

If one deck has order one, it covers its own owner family and no distinct
owner family.  The other two colours must therefore cover both of their owner
families, and the oriented two-colour descent of THM-800 forces both orders to
be one.

Assume no order is one and split by the number of order-two decks.

1. With three order-two decks, every owner needs an incoming cross edge whose
   label ratio is `+2` or `-2` in `F_13`.  A three-vertex digraph of positive
   indegree contains a 2- or 3-cycle, but neither a product of two reciprocal
   `+/-2` ratios nor a product of three `+/-2` ratios can be one modulo 13.
2. With two order-two decks, both must cross-cover the third owner, forcing
   their labels to be negatives.  They then do not cross-cover one another,
   while the third deck has capacity at most `1/3`; one owner remains short.
3. With exactly one order-two deck, it must cross-cover both other owners,
   forcing those two labels to be negatives.  For an order `D>=3`, own-owner
   and complementary-owner capacities are respectively

   ```text
   f(D)=ceil(2D/13)/D,       g(D)=floor(2D/13)/D.
   ```

   The two owner inequalities imply
   `(f+g)(D_2)+(f+g)(D_3)>=1`, whereas
   `(f+g)(D)<=3/7` for every `D>=3` (check `3<=D<=8`, then use
   `4D/13+1<=3D/7`).  This is impossible.
4. With no order-two deck, all three capacities are at most `1/3`, so equality
   forces three order-three decks and full cross incidence.  Mutual
   order-three incidence allows only ratios `5` and `8=5^(-1)`; three distinct
   labels cannot be pairwise related by those two ratios.

The pending replay will enumerate the exact residue/deck incidence for the
finite exceptional orders used above and stress the formulas over a long deck
range.  It will also record the obligation-incidence tournament and the
information lost by quotienting its half-open germs to runner vertices.
