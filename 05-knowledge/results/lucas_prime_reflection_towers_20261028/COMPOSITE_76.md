# The composite76 boundary: Frobenius survives, the nineteen-label mirror does not

2026-10-28. The ring decomposition and equivariance are **PROVED**; the
four-way ambiguity and eight-row table are **FINITE-EXACT**, exhaustively
checked on all76 legal period-nine words by `check.py`. Same-agent audit.

## Exact arithmetic rather than a repeated number

Write O=Z[phi], phi^2=phi+1. Since

    phi^9-1=20+34phi=2(10+17phi),
    Norm(10+17phi)=-19,

and10+17phi reduces to the unit phi modulo2, the two ideals are coprime.
The Chinese remainder theorem gives

    O/(phi^9-1) = F4 times F19.

The explicit map sends phi to(omega,5), where omega^2=omega+1 in F4.
Its additive group is C2 times C38, not C76. Thus no additive translation
has order76. A prime-degree cyclic argument already fails at its initial
transitivity step. This is an exact realization of4*19 in this problem,
not an identification with every other appearance of76 in the repository.

The audited all-odd-period reader identifies each of the76 legal words with

    (a,b)=(sum w_i omega^i, sum w_i 5^i).

Reflection R:i->-i satisfies

    a(Rw)=a(w)^2,

because omega^(-1)=omega^2 and squaring is additive in characteristic two.
So the four-state projection really is closed under reflection. In contrast,
for each fixed nonzero b in F19, the four words(a,b) have **four different**
reflected F19 coordinates. At b=0 all four reflected F19 coordinates are0.

Consequently no function F19->F19 alone computes this reflection. Given a
nonzero b, distinguishing its four possible outcomes requires a four-valued
tag, or at least two bits in a fixed-length binary code. Keeping a is enough.
This is a conditional ambiguity bound for this operation, not a universal
entropy bound and not a claim that every task needs the tag.

## An eight-chart presentation of all nonzero nineteen coordinates

The element5 has order9 modulo19. Every nonzero b is uniquely

    b=epsilon*5^i,     epsilon in {+1,-1}, i in Z/9.

Normalize the F4 tag as c=omega^(-i) a. Then a nonzero-b state is
(epsilon,c,i): eight chart types and nine phases each. Rotation of the word
adds1 to i, leaving(epsilon,c) fixed. Reflection reverses the phase, with
a chart-dependent shift given by these eight seed images:

| epsilon | c | a after reflection | b after reflection | epsilon' | shift j |
|---|---|---|---:|---|---:|
|+|0|0|16|+|7|
|+|1|1|1|+|0|
|+|omega|omega+1|5|+|1|
|+|omega+1|omega|12|-|6|
|-|0|0|18|-|0|
|-|1|1|12|-|6|
|-|omega|omega+1|7|+|6|
|-|omega+1|omega|10|-|5|

For an arbitrary phase, the resulting normalized coordinates are

    (epsilon', omega^(-j)c^2, j-i).

This follows from R rotation R=rotation^(-1), rather than requiring72
unrelated rules. The four b=0 states use the Frobenius rule a->a^2.
Two chart types are exchanged; six are fixed. They give six reflection-fixed
states, plus two on the b=0 fibre, recovering the total F6=8 fixed words.

The arithmetic translations do not preserve these multiplicative charts;
they require evaluating the retained(a,b) state and reassigning its chart.
Erasing that requirement would repeat the eleven-state extrapolation error.

## Why the prime hypothesis cannot be removed from the group theorem

There is a degree76 imprimitive control with a full76-cycle and an involution
fixing8 points, so1<8<sqrt(76) alone does not force A76. On pairs(b,j),
0<=b<19,0<=j<4, define

    T(b,j)=(b+1,j) for b<18, and T(18,j)=(0,j+1 mod4).
    R(b,j)=(b,j) for b<2, and R(b,j)=(b,j xor1) otherwise.

Both preserve the nineteen blocks of size4. T is a full cycle; R fixes
eight points and exchanges34 pairs. It has exactly the same cycle type as
the period-nine word reflection, so a relabeling identifies their involutions.
The generated group still preserves a nontrivial block system. This is a
hostile for **arbitrary** cyclic labelings at composite degree, not a claim
that the displayed T is golden addition.

## Connection to the marking calculation

The single-port marking quotient in `FULL_HOLONOMY.md` retains the first
nonzero ternary digit but loses a mod3 charge for pairs of ports. Here the
F19 projection retains one residue but loses a four-valued tag for a later
mirror. Both are exact failures of operation-compatible projection. Their
carriers and charges differ; neither is being called an isomorphism of the
other. The reusable move is to test the operation on every element of a
projection fibre and retain the smallest tag that separates its outputs.
