# Pairing the sevenfold neighbor rail with p and with odd sums

This follows SEVENFOLD_NEIGHBOR_RAIL.md. It certifies finite paired
examples and proves an all-order limitation of fixed transitive-block
clocks. It does not provide a universal construction of the quotient.

## Arithmetic and physical targets

On the m=1, k=3+3i slice of the two-port family, set

    q_i=(704·8^i+3)/7,         i>=0.

The upper tournament realizes 7q_i, with both 7q_i−2 and 7q_i+2 also
realized. The exact short return is

    q_(i+1)=q_i+7q_i−3,
    7q_(i+1)=8(7q_i)−21.

This is a direct numerical p+7p link: subtracting 3 moves to the next
lower target, and subtracting 21 from the eightfold upper count moves
to the next upper witness. It does not require, construct, or subtract
a count-21 tournament; 21 is a constant in the arithmetic. The upper
construction is physical at every i; the lower is a separate obligation.
Set p_j=q_(2j) below to retain the previous every-other-height notation.

The first eight consecutive q_i all have physical lower witnesses:

    q_i      witness
      101    F1(3,17)
      805    ordinal sum of count 35 and count 23
     6437    F1(35,123)
    51493    F1(1,36458)
   411941    F1(3,84871)
  3295525    F1(53,48470)
 26364197    F1(3,5713007)
210913573    F1(1,156879846).

The previous literal one-vertex lift supplies the 35 witness;
C3[TT_2,TT_3,1] has 23 paths. This verifies the first seven
transitions with actual tournaments at both ends, but not by a
single common gluing move.

If “p+7p” is read literally, the sum is 8p, which is even, so it
cannot itself be a tournament Hamilton count by Rédei oddness.
Its natural nearest odd targets are 8p−1 and 8p+1. Do not infer
these from an additive closure rule: ordinal joins multiply counts.

Here are exact F1 address certificates (see the inherited join
rule in POLYGON_COMPLETION_AND_PRIME_COUNTS.md). In every row
the upper witness 7p_j comes independently from the two-port rail.

| j | p_j | F1 addresses for p_j | 8p_j−1: addresses | 8p_j+1: addresses |
|--:|--:|---|---|---|
| 0 | 101 | (3,17) | 807: (3,206) | 809: (3,167) |
| 1 | 6437 | (35,123) | 51495: (3,11826) | 51497: (37,1118) |
| 2 | 411941 | (3,84871) | 3295527: (1,2283890) | 3295529: (3,666715) |
| 3 | 26364197 | (3,5713007) | 210913575: (1,105909313) | 210913577: (1,108439665) |

Only the p_j entries at j=0,2 are prime in this table. At j=0
the replay verifies the literal 101,707,807,809 tournaments by
subset DP; for the larger examples it uses the exact inherited
F1 separator formula, not enumeration of all Hamilton paths.

## Why fixed module clocks cannot supply the lower rail

Fix a quotient and an acyclic set of marked vertices. Substitute
nonempty transitive modules with lengths a_i j+b_i at those marks,
where a_i,b_i are integers and lengths stay positive for large j.
After rebasing j at an integer if necessary the b_i can be taken
positive. Uniform arcs between modules are prescribed by the quotient.

Compressing maximal visits gives finitely many legal marked words:
between unmarked vertices, the marks are a DAG, so no block recurs
in that segment. Each legal word with visit multiplicities r_i
contributes the product of r_i! S(a_i j+b_i,r_i). Each such term has
an inclusion–exclusion expansion in powers (t^a_i)^j with integer
coefficients (the shifted powers t^b_i are integral). Thus the count
is a finite sum of integer multiples of positive integer bases to
the power j, with an integer constant term. This includes every
fixed finite bank of such templates along fixed affine length schedules.

By contrast, q_i has the same fractional coefficient at base 8 as p_j does at base 64. In particular p_j has coefficient 704/7 at base 64 and constant 3/7;
8p_j±1 have coefficient 5632/7 and constant 24/7±1. None can equal
one of those integer-coefficient sums at infinitely many j. Indeed,
subtract the two expansions: at least one rational coefficient is
nonzero, and the largest base with nonzero coefficient eventually
dominates all smaller bases, so there can be only finitely many zeros.
Rebasing leaves the denominator 7 (coprime to 64) intact. This does
not rule out choosing lengths nonaffinely, recursively changing the
quotient, or using marked cycles or nontransitive modules.

## The change of state to look for

The short return also has a positive two-coordinate specification by
retaining the separate increment d_i=7q_i−3:

    (q_0,d_0)=(101,704),
    (q_(i+1),d_(i+1))=(q_i+d_i,8d_i).

Here 7q_i is an actual upper count, while d_i is only a requested
marked response, not necessarily a tournament count (it is even).
A physical compiler would need a marked coordinate carrying d_i and
another that adds it to the existing path count. Carrying
only a transitive-module width cannot do this on the stated schedules.
These integer recurrences are specifications for a physical
construction, not a proof that one exists.

Run node paired_sevenfold_sum_replay.js.
