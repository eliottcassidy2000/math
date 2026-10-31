# Five source-preserving forms for any unpaid fourteenth-layer candidate

**PROVED conditional only on this packet's independently checked E10
certificate and the native return identities.** E14 remains **OPEN**.
This extends the one-core argument in incoming Session120's
[TEN_RETURN_GRAMMAR.md](../session120_major_frontier_interfaces_20261030/collatz/TEN_RETURN_GRAMMAR.md)
to the new fees. It does not claim a new general return-decomposition method.

Use the notation of COLLATZ.md C5-C6. Set

    eta=183105426/48828125,
    T11=19,531,250, T12=97,656,250,
    d=lambda14-lambda12=3/5^13.

Let w be a non-one primitive13 excursion. Its price is
delta13(w)=lambda14*A_w-S13(w). The two remaining failure thresholds are
delta13<1 at tangent0 and delta13<4/5 at tangent-(4,4).

## Two cuts retain exactly one non-one core

First cut w into primitive12 excursions. Write b12 for the number of
non-one constituents and a12 for the all-one constituents. Since w is
non-one, b12>=1. C6 gives each non-one12 excursion fee at slope lambda13
at least eta-2, and each all-one12 excursion fee1/5. Subtract w's one
terminal level13 charge, then raise the slope to lambda14. If b12>=2,
the resulting price is at least2(eta-2)-1>1. Thus any unpaid w has b12=1.

Call that unique non-one12 constituent c12. At slope lambda14, each
all-one12 constituent pays

    lambda14*T12-S12(1^T12)=1/5+1/25=6/25.

Cut c12 into primitive11 excursions, with b11 non-one constituents and
a11 all-one constituents. Each non-one11 constituent pays at least eta-1
at slope lambda12. Each all-one11 constituent pays, at slope lambda14,

    lambda14*T11-S11(1^T11)=1/5+1/25+1/125=31/125.

There is exactly one terminal level12 return in c12 and one terminal
level13 return in w. Therefore

    delta13(w) >= b11*(eta-1)-2
                  +(31/125)*a11+(6/25)*a12
                  +d*sum_(non-one11 cores c) A_c.        (1)

If b11>=2, its first term alone is2(eta-1)-2>1. Hence b11=1.
Let c be this single primitive11 non-one core. Then(1) becomes

    delta13(w) >= eta-3+(31/125)*a11+(6/25)*a12+d*A_c.    (2)

All discarded terms are nonnegative. The cuts use actual return points,
so every lower identity test is preserved and each clock unit is counted
once. No entry phase is replaced by a cheaper or maximizing phase.

## Exact candidate grammar and clock cutoffs

For delta13<1, inequality(2) forces

    (a11,a12) in {(0,0),(1,0),(0,1)}.

Indeed either coordinate at least2, or both at least1, already raises
the constant part above1. Each possible extra all-one block is therefore
before or after the unique core. The five necessary forms and strict
clock bounds are:

| candidate w | largest possible integer clock A_c |
|---|---:|
| c | 101,725,616 |
| 1^T11 c | 814,158 |
| c 1^T11 | 814,158 |
| 1^T12 c | 4,069,366 |
| c 1^T12 | 4,069,366 |

For delta13<4/5, both a11 and a12 must be zero. Thus w=c itself is
primitive11, with A_c<=20,345,408.

The rational strict ceilings before rounding are respectively
305176850/3,2442475/3,12208100/3,61036225/3. They come directly from
(target-eta+3-(31/125)a11-(6/25)a12)/d, not a sampled stopping-time fit.

These are necessary forms, **not** freely realizable constructions.
The core must still be a positive native11 word, the complete w must be
native13, and its actual full affine map modulo5^14 must have the required
tangent. Prefix and suffix placements remain distinct: although their
endpoint congruences here agree, intermediate source charges can differ.
The grammar therefore does not merge their clocks or path histories.

## Exact port constraints and checked scope

Put q=5^14. The two required endpoint maps modulo q are

    (m,t)=(1,0), or (1+5^13,5^13),

respectively. A block1^T_e has map (a,a-1), where a=(3/2)^T_e mod q.
Chronological composition is(m2*m1,m2*t1+t2). Solve the corresponding
left or right equation for the required core map; do not erase it.

native14_controls.py checks these modular compositions, every exact
fraction cutoff, both one-core eliminations, and the first excluded integer
clock in each form. It does not enumerate cores up to those large bounds
or assert that any candidate exists. The remaining task is to exclude the
feasible native words in this grammar or produce a cheaper compatible
certificate for all of them. All-layer content control and universal
Collatz remain open.

## Reusing the runner floor-sum engine on the integer clock

There is an exact necessary filter much smaller than a word enumeration.
Since3 has order M=4*5^13 modulo q=5^14, write

    a=log_3(2)=2883827743 mod M,
    b=log_3(required core multiplier) mod M.

The order follows from ord_5(3)=4 and v5(3^4-1)=1; the displayed logarithm
is checked by modular exponentiation. A core of length L and clock A has
multiplier3^L/2^A, so it must satisfy

    L = a*A+b mod M,               1<=L<=A-1.             (3)

For the clock ranges above A<M, hence at most one L is possible per A.
The exact number of possible pairs with1<=A<=N is

    sum_(A=1)^N [ floor((a*A+b-1)/M)
                 -floor(((a-1)*A+b)/M) ].               (4)

The difference counts multiples of M between a*A+b-A and a*A+b-1,
with precisely the required open/closed endpoints. Each sum is an ordinary
Euclidean floor sum, the same arithmetic engine used by the two-tail runner
compiler. Binary search in its cumulative count extracts pairs without
enumerating every integer clock. Here the source and target of the transfer
are explicit lattice counts, not an asserted equivalence of conjectures.

The exact counts, before testing translation, native status or positive
word realizability, are:

| form / threshold | necessary clock-length pairs | first pair(A,L) |
|---|---:|---:|
| bare c, price<1 | 1,059,639 | (177490,167570) |
| either placement of1^T11, price<1 | 72 per placement | (73697,50871) |
| either placement of1^T12, price<1 | 1,697 per placement | (128964,110752) |
| bare c, price<4/5 | 42,387 | (70996,67028) |

native_clock_floor.py checks the modular logarithms, compares(4) with4,320
literal small-range counts, and explicitly extracts all72 pairs for each
T11 placement. A pair represents potentially many valuation words. Matching
its multiplier does not match the required translation or certify a native
return. Those constraints remain attached to every row in the receipt.
