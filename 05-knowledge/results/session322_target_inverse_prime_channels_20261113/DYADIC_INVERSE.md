# An exact target-driven inverse in the two selective charts

The [selective contexts from Session 321](../session321_selective_returns_and_fano_wiring_20261113/SELECTIVE_RAYS.md)
give strong outputs of the form
```
N = A Q_V(2) − h H(V),         h=H(U),
A_13 = -3h + 6 Q_U(2) + 2 Q_U(3),
A_15 = -3h + 6 Q_U(2) + 6 Q_U(3).
```
Here U,V are nonempty, Q counts ordered optional slots, and the cut
is the appropriate named two-module quotient. This note gives an
exact inverse for a growing, explicitly constructible family, not
a proof that every target passes it.

## A unique dyadic shell

Assume Q_V(2)=2^m. Any nonempty tournament has Q_V(2)>=2H(V),
since either slot can hold a Hamilton path. Also A>h>0. Therefore
```
A 2^(m−1) < (A−h/2) 2^m <= N < A 2^m.
```
The absolute target and the source U determine **the unique possible
m**, the first integer with A2^m>N. The required reservoir count is
then the single integer D=(A2^m−N)/h. There is no free final bit
after taking a congruence. If D is not a positive odd integer, the
source fails.

## A multiplicatively coded physical reservoir

Use a bank of strong factors with the following (H,log_2 Q(2)) pairs:
```
(5,5), (19,8), (217,12), (223,12), (599,13), (761,14), (773,14).
```
The replay supplies literal outgoing rows and checks both counts
by subset DP. This is a bank of witnesses, not a census of dyadic
slot profiles. The H values are pairwise coprime (217=7·31,
the other six are prime). Form V by ordinally concatenating e_i
copies of these factors, and pad by singleton components. Both H
and Q(2) multiply on an ordinal sum, independently of the order.
Thus the inverse for this bank is exact:
1. Extract the maximal powers of each h_i from D by repeated
   division. Reject if a residual factor remains.
2. Check the cost sum sum_i e_i m_i <= m.
3. Assemble these components and exactly m−sum_i e_i m_i singletons.

Positive D and m>=1 exclude an empty reservoir, because each
component or singleton has positive dyadic cost. No primality
heuristic is needed by the inverse.

For example, given N=3821 and U a singleton, the A=15 context
chooses m=8 and D=19: the order-six factor is forced in this bank.
For N=61223 it chooses m=12 and D=217, forcing the other factor.
Given N=6907 and A=27 it chooses m=8 and D=5, forcing the
order-four factor and three singleton pads. In each case the
literal substitution is verified independently by subset DP.
The target is used first, not recovered after choosing a length.

## Choosing U as well as V

One exact compiler allows U=TT_l, l>=1, in both named contexts.
Then h=1 and
```
A_13(l)=-3+6·2^l+2·3^l,
A_15(l)=-3+6·2^l+6·3^l.
```
Since N>A for nonempty V, only O(log N) values of l can pass.
For each, compute the unique m and factor D within the fixed bank.
This enumerates all witnesses in the stated grammar, up to the
irrelevant permutation of ordinal components, and either returns a
literal strong tournament or certifies no hit in that grammar.
The same method extends to any prescribed finite source bank of
SCC profiles: precompute their (H,Q2,Q3), and enumerate bounded
multiplicities, not orders of blocks.

The supply limitation is sharp in its own scale. For fixed m,
the possible D are indexed injectively by e_i>=0 with weighted
cost <=m. Their number is
m^7/(7!·product m_i)+O(m^6), by comparing lattice points with the
simplex after enclosing unit boxes. This grows unboundedly, but
only polynomially in m, against exponentially many target integers
in the shell. The finite-bank inverse is exact when it succeeds;
it cannot by itself give Chan–Pak completeness.

For a prime target N there is an additional universal screen, not
restricted to the bank. Any common factor of H(V),Q_V(2) divides
N, as does any common factor of A,h. Because the output is strictly
larger than H(V) and h (use Q2>=2H(V) and A>h), both gcds must be 1.
In an ordinal V this also excludes a prime shared by the H of one
factor and the Q2 of another. For example a directed triangle
cannot appear as a reservoir SCC at a prime endpoint. The dyadic
bank automatically avoids the V-side obstruction.

## The small prime split gives a native modular pair, but not dense endings

The size-four and size-six factors above have H=5 and 19 respectively,
one prime on each side of the proposed first-ten cutoff. They suffice
for **every odd residue modulo 2^k**, using either singleton source
chart. For k>=2, powers of 5 give exactly the units congruent to 1
modulo 4: v_2(5^(2^j)−1)=j+2 by induction and squaring. Since 19 is
3 modulo 4, choose b in {0,1} and then solve 5^a 19^b = −N mod 2^k
with 0<=a<2^(k−2). Choose m at least max(k,5a+8b,1) by ordinal
padding. Then A2^m−5^a 19^b has the required residue and has a
literal strong realization. For k=1 parity alone suffices.

There is an exact cost caveat. With the full seven-factor bank and
cost bound m, at most O(m^7) distinct residues can occur through
exponents of that cost. Covering all 2^(k−1) odd residues at a
single m>=k therefore needs m at least on the order of 2^(k/7).
An exact absolute target fixes m by its shell, so arbitrary
additional padding cannot pay this cost. This is a physical
modular construction, not an exact arbitrary-target algorithm.

A small alignment with the independent prime graph is exact as
an integer check and nothing more. With singleton U, the two
coefficients are 15 and 27. Their first transitive reservoirs give
29 and 53; the next A=15 reservoir gives 59. These are respectively
the end of the first ten primes, the first unrestricted connector
threshold, and the first Ogg-only connector threshold in the
[prime-square-sum computation](PRIME_SELECTION_4_PLUS_6.md).
The overlap does not identify the graph operations or explain
the Monster filtration.
