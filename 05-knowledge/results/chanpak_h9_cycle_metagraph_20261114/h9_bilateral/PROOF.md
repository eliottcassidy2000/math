# An h=9 unit pair with arbitrary bilateral readers and a guarded three-arm rail

Let h=H(B), f=P_B(2), and let P,S denote the total and odd-minus-even
weighted starts at the vertices beaten by an apex. Let E,T denote the
same weighted ends at the vertices that beat the apex. Each nonempty
directed path has weight H(complement). Write K=H(B+apex).
These are the coordinates of the [ordinal cut transfer](
../../chanpak_selective_target_20261114/CUT_STATE.md).

There are three actual marked eight-vertex punctures (outgoing row masks,
then the apex out-mask):

| Rows | U | (h,f,P,S,E,T; K) |
|---|---:|---|
| (254,220,248,240,224,194,128,0) | 182 | (9,672,245,9,209,9;89) |
| (254,252,184,240,224,192,132,0) | 202 | (9,672,245,9,209,9;91) |
| (238,252,248,240,225,192,128,0) | 169 | (9,672,245,9,209,7;93) |

Take an arbitrary transitive left reader (ordinary binary contact
address x) and right reader (reverse-complement contact address y).
Write g(x) for the transitive-plus-apex count at address x, with g(0)=1,
and p(x)=x mod 2. The exact composition rule gives
$$
H=K+9(g(x)+g(y)-2)+209x+245y+672xy+T p(x)+9p(y).
$$
For clarity, g is unchanged by complement-reversal of its contact
word: reverse all arcs and relabel the underlying transitive order.
The left-only contribution is 9(g(x)-1)+209x+T p(x). The start total
after prepending that reader is 672x+245, its bias remains 9, and the
right-only contribution is 9(g(y)-1)+(672x+245)y+9p(y).
This proves the formula without an assumption on the tape bits.

Therefore the first two punctures give a unit pair for arbitrary
painted tapes in both directions. The first **three** give a unit
triple for every even left address and arbitrary right address.
It is not necessary to match the end bias T when the left reader
has signed start weight zero.

## Target-directed decoding in an independently unbounded triple

Choose independent heights i,j≥1 and set x=10·2^i−2 and y=2·2^j−2.
The recurrences g(2u+1)=g(u) and g(2u)=g(u)+u+p(u) give
g(x)=5+5·2^i and g(y)=1+2^j. Consequently the literal window is
$$
13440,2^{i+j}-11305,2^i-2189,2^j+{1905,1907,1909}.
$$
The coefficient 13440 is divisible by 2^7, the linear coefficients
are odd, and 11305+2189 has valuation 1. Thus for a prescribed target
N and each possible offset d, set t=v₂(N−d); any match has either
i=t<j, j=t<i, or i=j=t−1. In the first two cases ordinary linear
division computes the other power of two; test the indicated height
inequality. In the tied case directly test the inferred pair.
This finds all solutions to this rail without walking the tape words
or factoring a large integer. It does not decide all Chan–Pak targets.

At (i,j)=(1,4) the triple is 374351,374353,374355. The first count is
prime, and the second is 7·53479, where 53479 is prime; these are
checked by trial division and the three underlying 19-vertex counts
are recounted directly in C++. For all i≡1 mod 42, j≡4 mod 42, the
first two counts are coprime to 3 and the middle count has exactly
one factor 7. Indeed 2^42≡1 mod 147, so the outputs have the same
residues modulo 147 as at the displayed witness. The low count is
also coprime to 7. The guards are not a claim of infinitely many
primes.

A second literal two-arm rail uses x=2^i−1, y=7·2^j−1:
$$
4704,2^{i+j}-463,2^i-2989,2^j+{361,363}.
$$
At i=j=2 it gives 61817=7·8831 and 61819 prime. For i≡2 mod 42 and
j even the lower output has precisely one factor 7 and both are
coprime to 3: the coefficients of 2^j and 2^(i+j) vanish mod 49,
and their mod-three constraints see only the parity of j.

## Every graph on these rails has exactly nine Hamiltonian cycles

Use the binary addresses without leading zeros, so a positive
left address means the apex beats the very first vertex. A
positive reverse-complement right address means the last vertex
beats the apex. Both outside ordinal transitive blocks are
nonempty. In the puncture every Hamiltonian path must start at
its first singleton component and end at its last one. Thus every
one of the h=9 Hamiltonian paths of the puncture closes by adding
the apex, and every Hamiltonian cycle (counted up to rotation)
is obtained in this unique way by deleting the apex. There are
exactly nine cycles on all arms for all positive tape addresses.
If unnormalized tapes are used, these endpoint contact guards
must be checked rather than inferred from the numerical address.
The replay recounts these cycles directly by pinning the apex
as the first vertex in a second subset DP.

## Scope

This is a corrected h=9 family with actual bilateral flexibility,
unlike the single-triangle h=3 fiber. Correct placements in this
fixed family can be inferred from absolute targets, not merely their
small residues. A fixed finite bank of such exponentially growing
two-axis rails nevertheless has only O(log² X) distinct counts up
to X. What remains for Chan–Pak is to make the physical profile
bank change with the target, or find a different family that can
express the missing absolute counts.

Run python replay.py and g++ -O3 -std=c++17 recount.cpp -o /tmp/h9bilateral
then /tmp/h9bilateral to reproduce the checks.
