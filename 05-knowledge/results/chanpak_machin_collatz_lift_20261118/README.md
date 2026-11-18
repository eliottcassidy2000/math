# Machin's unit correction and a physical Collatz lift for tournament pairs

**Proved:** a classification of unit corrections for repeated Gaussian
squaring (using Ljunggren); an exact Machin-to-Farey rounding exercise; a
physical centered pair at every width on the next dyadic rail, with an
arbitrary tournament source; and a separation between approximate and
exact slot ratios. This does not prove Chan–Pak or Collatz.

## 1. One exact residual in the irrational calculation

Machin's identity follows by a Gaussian multiplication, not estimation:

    (5+i)^4 = 4(119+120i) = 2(1+i)(239+i).
    pi/4 = 4 atan(1/5) - atan(1/239).

Likewise (3+i)^2=2(4+3i)=(1+i)(7-i), so
pi/4=2 atan(1/3)+atan(1/7). These give the specific parameter
triples (4,5,239) and (2,3,7) with opposite correction signs.
The arguments in these two cases are in the first quadrant, so there
is no branch ambiguity. A historical source with the formula printed
explicitly is [Shanks and Wrench](https://www.engert.us/erwin/miscellaneous/William%20Shanks%20707%20digits.pdf).

There is a rigidity theorem hiding in the apparently small step.
For odd a>=3, let U+iV=(a+i)^(2^j)/2^(2^(j-1)), j>=1. This is
integral (a+i is divisible by 1+i), with norm b^(2^j) for
b=(a^2+1)/2 odd. The signed correction relative to the diagonal uses
(V-U)/(U+V). The unit condition |U-V|=1 implies

    (U+V)^2 + 1 = 2 b^(2^j).

For j>=2 the right side is twice the fourth power of
b^(2^(j-2)). By
[Steiner and Tzanakis's proof of Ljunggren's theorem](https://users.math.uoc.gr/~tzanakis/Papers/LjunggrenEq.pdf),
the positive solutions of x^2+1=2y^4 are (1,1),(239,13).
Since b>=5, only j=2,a=5 is possible. For j=1 the equation
|(a^2-1)/2-a|=1 gives a=3 directly. Thus within this precise
one-seed, repeated-squaring model the two small numerators are not
generic, even though the exact angles always compose. This is not a
classification of general Machin-like formulas.

## 2. The exact role of 71, 21, 22 and 7

The two Archimedean bounds 223/71 < pi < 22/7 satisfy
22*71-223*7=1. Put

    r_k=(22k+3)/(7k+1) = 22/7 - 1/[7(7k+1)].

At k=10=(21-1)/2 the denominator is 71 and the numerator is 223.
Thus the relationship between 21 and 71 is literal under the odd
coordinate N=2k+1, but says nothing by itself about a Hamiltonian
count. Also 22=3*7+1, the unaccelerated Collatz numerator at 7.
The first r_k above pi is not r_10: it is r_16=355/113, with
r_15=333/106 below pi. Five alternating terms of each arctangent in
Machin's formula, bracketed by adding the next term, certify these
inequalities in exact rational arithmetic. This is an example of
using upper and lower error bounds to certify a discrete address,
rather than treating a good approximation as equality.

## 3. Odd rows, a +2 rung, and two exact ways to step down

The [prior two-mode recursion](../session35_dimensions_dualities_20261023/recursions/RESULT.md)
uses positive integers 2^k b with b odd: +2 on the base row and
doubling upward. On addresses those operations commute, but +2 on
row k of odd seeds adds 2^(k+1) at height k, not 2.

For the current application put unit rungs at N_s=2^k b+s, s=±1,
k>=2, b positive odd. Every positive odd N>=3 belongs to exactly one
rung: choose s with N=s mod 4, then compute k=v_2(N-s) and test its
odd part for b. The collated addresses are

    N       s     k     b        other endpoint
    7      -1     3     1         9
    21     +1     2     5        19
    63     -1     6     1        65
    71     -1     3     9        73
    239    -1     4    15       241

A fixed rail cannot contain both 7 and 21 on the same odd core,
and 7 and 63 sharing a core makes plain why full physical width
lifting would need an exceptional boundary.

Use the odd shortcuts of 3x+1 and 3x-1 together. For a rung, choose
3N_s-s on each branch and remove its exact factor 2: the endpoints
become 3*2^(k-1)b+s. Or choose 3N_s+s and remove its exact factor 4:
for k>=3 they become 3*2^(k-2)b+s. In each case the same signed
unit offset survives at both ends. These are equalities, including
the valuations, not merely agreement modulo powers of two.
Positive-only Collatz takes different choices on the two signs;
this does not prove convergence of either iteration.

## 4. A physical lift from slope 6 to slope 9

The [previous two-reservoir construction](../chanpak_exact_two_reservoir_center_20261118/README.md)
with one singleton reservoir gives 6P_2(X)±H(X).
A new pair of rooted six-vertex templates, root indexed 0, has
outgoing bit rows

    B- = (8,37,1,6,15,29)       H=17
    B+ = (12,33,2,6,15,29)      H=19.

Replace the root in each by any nonempty tournament X with uniform
cross arcs. The two templates are strong, so both hosts are strong.
The rooted quotient uses each outside vertex once. Count the quotient
words by their number j of root runs: the two profiles are
(17,9,0,...) and (19,9,0,...). A run between two occurrences of the
root requires at least two intervening outside vertices, so j<=3;
the finite enumeration is supplied. Filling runs by Q_j(X), ordered
nonempty path covers, gives exactly

    H(B_s[X]) = 9P_2(X)+s H(X),

since Q_2=P_2-2H. Thus the center goes from 6P_2(X) to 9P_2(X)
while the same absolute offset H(X) survives. At offset 1 this is
precisely the shared one-level Collatz transfer above, for every
transitive X and not just on one observed pair. With general offset
h, the corresponding arithmetic operation is (3N_s-s h)/2.
The change uses two core arcs, one at the root, so after inflation
it is not generally a two-arc edit. No claim about native metagraph
edges is made.

For transitive X of width m the exact endpoints are 9*2^m±1.
Given an integer target, choose its sign by the congruence test and
test the exact odd part and exponent; the construction then returns
both endpoints. At m=3 the new lower count is 71. In the *old*
two-reservoir template 71 cannot occur: any offset has h<=6 and,
as 3P_2(X)P_2(Y) is divisible by 12, h must be 1 with negative sign.
That would force both reservoirs transitive (count 1), so their
P_2 product is a power of two, not 24. The new core genuinely closes
a target excluded by the old physical language. There is also a
separate literal one-arc witness for 71/73 with outgoing rows
(6,8,106,17,103,11,43) and (6,8,122,17,99,11,43).

The lower ray includes 35=7*5 and 287=7*41, and is divisible by
7 precisely at m=2 mod 3, by the order of 2 modulo 7. The upper ray
never is. The next simultaneous arithmetic shift requires an
unbounded width family with odd center 27; no such lift is proved here.

## 5. An approximation that cannot stand in for exact placement

The ratio P_2/(2H) cannot be 223/71 for any tournament, even allowing
a common scaling of numerator and denominator. If there were such a
profile then H=71t with t odd and P_2=446t, so P_2 is 2 mod 4. The
odd-cycle packing formula of
[Grinberg and Stanley](https://arxiv.org/abs/2307.05569)
has only powers of the same parity as the order n, with no constant
for a nonempty tournament. Hence n is odd. For n<=5, H<71:
at n=5 there are at most ten triangles, and the number of Hamiltonian
cycles is at most H/5 by breaking each cycle at its five arcs.
The packing expansion gives H<=1+20+2H/5, hence H<=35.
Orders one and three are immediate.

For odd n>=7, the coefficient of the linear term is twice the number
of Hamiltonian cycles, at most 2H/n. All other terms at slot 2 carry
at least a factor 8. Thus P_2/(2H)>=4-6/n>=22/7, strictly above
223/71. The count 71 itself is perfectly realizable: a missing
profile is not a missing count. The same estimate gives P_2>=4H
for every non-singleton tournament (even n follows from lowest
degree 2). In particular the earlier blue norm read (P_2-H)^2+2H^2
cannot implement a *negative Pell unit correction*
(P_2-H)^2-2H^2=-1 except at the singleton, where the count is 3.
An arithmetic solution on a thin norm surface need not be a
physical tournament profile.

Conversely the ratios can approximate irrational numbers if one
allows dyadic renormalization. A four-vertex reservoir with (H,P_2)
=(5,32), repeated j times by ordinal sum, has ratio (32/5)^j.
Divide it by a suitable power of two to place it in [2,4):
the resulting ratios are dense by the irrationality of log_2 5
and the elementary irrational-rotation theorem. Near-returns such
as 2^7/5^3=128/125 give small ratio changes, but not arbitrary exact
integer corrections. On that one-seed family the reduced odd
denominator is always 5^j. The separation between normalized
approximation and an integer output is the proof issue, not just
the numerical precision.

Run `python3 verify.py` for exact rational bounds, rooted skeletons,
literal lifted counts and the 71/73 one-arc witness.
