---
id: THM-859
title: Hamming-six dilation conjugacy and the order-one ramification gate
status: PROVED STRUCTURAL — common dilation conjugates the exact component action; effective order is the exact ramification degree; every H6 common-sheet branch containing order one is a uniformly finite family of finite exact recursions
source: codex-2026-07-15-S10 continuation-congruence audit
depends_on: [THM-815, THM-823, THM-840, THM-857, THM-858]
related: [THM-844, THM-847, HYP-6820]
---

# THM-859 — dilation conjugacy and the order-one gate

For a finite positive speed set `Q`, put

```text
Safe(q)={t in R/Z:||qt||>1/13},
E(Q)=intersection_(q in Q) Safe(q),
M(Q)=max_t min_(q in Q)||qt||.                         (1)
```

For a positive integer `c`, write `pi_c(t)=ct` on `R/Z`.

## Theorem

### A. Exact dilation conjugacy

For every `Q` and every positive integer `c`,

```text
Safe(cq)=pi_c^(-1)(Safe(q)),
E(cQ)=pi_c^(-1)(E(Q)),
M(cQ)=M(Q).                                             (2)
```

More strongly, dilation conjugates every monotone insertion:

```text
pi_c^(-1)(E intersect Safe(u))
 =pi_c^(-1)(E) intersect Safe(cu).                      (3)
```

Every component of `E` has exactly `c` inverse-image components, each with
length divided by `c`.  THM-815's exact next-speed cap therefore scales by
`c`, and the complete THM-857 recursion restricted to speeds divisible by
`c` is canonically isomorphic to its scale-one tree.

Consequently, let `R subset [12]`, `|R|=6`, and consider a proper common-scale
packet

```text
A=c([12] minus R) union {w_r:r in R},
w_r=cr (mod 13c).                                      (4)
```

If every effective order

```text
D_r=c/gcd(c,w_r)                                       (5)
```

equals one, then `w_r=c(r+13h_r)` with `h_r>=1`.  The only tight packet in
this stratum is

```text
A=2c[12],                                               (6)
```

obtained from the six odd labels and `h_r=1`; every other packet in (4) is
strictly loose.

### B. Effective order is the continuation-ramification degree

The deck quotient `t~t+1/c` is preserved by insertion of the `w`-comb if and
only if

```text
c divides w,             equivalently D=c/gcd(c,w)=1.  (7)
```

If `D>1`, the translates of `Safe(w)` over one `1/c`-deck orbit consist of
exactly `D` distinct masks.  Thus `D` is literally the ramification degree of
the failed scale-quotient continuation congruence.  An order-one insertion
updates the scale-one quotient deterministically; an order-`D` insertion
requires a `Z/DZ` phase sidecar.

### C. The H6 order-one gate is uniformly finite

Consider the six-colour analogue of THM-823's common-sheet presentation.  If
at least one colour has effective order one, deleting such a colour leaves a
five-colour common-sheet presentation on the other five owners.

1. If at least two colours have order one, then either all six orders are one,
   already closed by Part A, or exactly two orders are one and the other four
   are order three on a coset of `{1,5,8,12}`, with equal units on opposite
   pairs.  The latter branch has exactly

   ```text
   3*C(8,2)*4=336                                      (8)
   ```

   labelled order/unit contexts at the order-three common sheet.

2. If exactly one colour has order one, the other five orders form a
   no-order-one THM-858 presentation.  Hence they are `{2,3,7}`-smooth, have
   no private maximal prime power, obey every relative-ramification cut, and
   satisfy

   ```text
   min D_i<=21,                 max D_i<=10,584.        (9)
   ```

   If all five are at most `21`, THM-823/858 reduce them to the `96`
   all-order-three H5 contexts; choosing the sixth order-one label gives
   exactly `96*7=672` marked H6 contexts.

Every fixed context in either case has a finite exact component recursion,
uniformly in its metric lift heights.  This makes the complete H6 order-one
gate finite-decidable.  It does **not** assert that the `336`, `672`, or the
remaining finite-strip context languages are empty.

## Proof

### 1. Dilation and the cap tree

The first identity in (2) follows from
`||(cq)t||=||q pi_c(t)||`; intersection gives the second.  Since `pi_c` is
surjective, maximizing the same minimum over its domain gives the third.
Identity (3) is then literal set algebra.

The inverse image of an open interval under the degree-`c` circle cover is a
disjoint union of `c` intervals of reciprocal length.  If THM-815's cap at a
component of length `L` is the real number `B(L,m)`, the lifted cap is
`cB(L,m)`.  For an integer candidate `u`,

```text
cu<=floor(cB)  iff  u<=B  iff  u<=floor(B).             (10)
```

Candidate residues and their numerical order are also carried by
`u->cu`, proving conjugacy of the whole restricted action tree.  Under (5),
`D_r=1` says `c|w_r`; dividing (4) by `c` gives the exact THM-857 packet.
Equation (6) and strict looseness of every other row follow from THM-857 and
(2).

### 2. The sharp operation boundary

Translation by `1/c` changes the phase of the `w`-comb by `w/c`.  The centred
safe arc `(1/13,12/13)` has trivial rotational stabilizer, and multiplication
by `w` maps the circle onto itself.  Therefore `Safe(w)` is invariant under
that deck translation exactly when `w/c=0 mod 1`, proving (7).

More generally, the number of distinct masks is the additive order of
`w/c mod 1`, namely `c/gcd(c,w)=D`.  This proves Part B.  THM-840 explains why
one cannot repair a ramified replacement by restoring, deleting, and
reinserting from endpoint data alone: exact endpoints are insertion-Markov,
not deletion- or replacement-Markov without the labelled tooth bank.

### 3. Reduction of an order-one colour

An order-one colour covers every sheet at its own owner and no sheet at any
other owner, directly from THM-823's oriented sheet formula.  Its order does
not change the common lcm.  After deleting an order-one colour `b`, every
remaining owner must therefore be covered by the other five colours, so they
form an H5 common-sheet presentation.

If another order-one colour remains, THM-823(E) says that the H5 row is either
all order one or is one order-one colour plus an order-three quartet on one of
the three multiplicative cosets of `{1,5,8,12}`.  Repeating the deletion for
either order-one colour gives exactly the two alternatives in Part C.1.  In
the mixed case choose the coset, the two order-one labels in its eight-element
complement, and one of the four opposite-pair unit words.  This gives (8).

If only one order-one colour existed, the remaining H5 row has none.  THM-858,
including its dyadic equality refinement, gives (9) and the stated carrier
conditions.  In the sub-bank through `21`, THM-858 reduces the no-order-one
row to THM-823's `96` all-order-three contexts.  Five labels are used, leaving
seven choices for the marked order-one label, which proves the count `672`.

Finally fix any one of the finitely many order/unit contexts.  If its original
common scale is `c_0`, put `C=lcm(D_i)` and `g=c_0/C`; every `D_i` divides
`c_0`, so `g` is an integer.  Divide every speed by the common factor `g`.
Every replacement then has the form

```text
w_i=(C/D_i)u_i,
u_i=D_i r_i (mod 13),       u_i=e_i (mod D_i),          (11)
```

so it lies in one labelled arithmetic progression modulo `13C`.  Literal
residual components, the remaining progressions, and the last speed are
Markov under insertion.  Before the sixth insertion the prefix has at most
eleven speeds; the settled lower-runner bound makes its strict-safe set
nonempty.  THM-815 then gives a finite cap for the next member of every
remaining progression, with at most six combs left.  Induction makes the
whole fixed-context tree finite, proving Part C.

## Tournament analysis and challenged vertices

The theorem-bearing relation is a cyclic group action, not a binary runner
tournament.  If colours are used as vertices, the raw observable `D_i-D_j`
and the deck-mask-count observable compare the same scalar orders; ties may be
resolved by increasing label, giving a transitive planning order but losing
which phase mask acts on which component.  No finite fingerprint is claimed
for the still-unenumerated strip (9).

The useful alternate vertices are deck translates or legal insertion
operations.  Their exact sidecar is the `Z/DZ` phase action on literal
components.  Quotienting to bare colours, pairwise order comparisons, or a
transitive tournament preserves an insertion priority but destroys the
continuation predicate.  This is the operation-theoretic reason the H6
frontier begins at ramification, not at common dilation.

## Scope guardrail

THM-859 closes all unramified common dilations of THM-857 and proves that the
H6 gate containing an order-one colour is uniformly finite-decidable.  It
does not evaluate the `336` or `672` ramified metric languages, enumerate the
full order strip (9), handle H6 presentations with no order-one colour, cross
the seven-comb wall, or settle global `n=12` sporadic-branch emptiness.
