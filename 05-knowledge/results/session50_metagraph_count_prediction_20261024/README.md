# Session50: recovering tournament counts from metagraph moments

24 October 2026. First checkpoint; marked growth and structural obstructions
are being integrated separately.

The [moment theorem](mass-predictions/THEOREM.md) turns authenticated positive
odd path-fibre sizes into exact arithmetic information about A000568. Its
[independent audit](audit-moments/AUDIT.md) checks the source atlas, the moment
normalization, the global inequalities, and the isolated decoder.

Write h=H(T)/|Aut(T)| for an ordinary class fibre, and K_r for the number of
ordered (r+1)-tuples of fixed-path words in the same ordinary class. Thus
K_r=sum_classes h^(r+1). On a converse-merged vertex of literal port mass mu,
h is mu when mu is odd and mu/2 otherwise. Loops count twice in that mass.

The first k moments K_0,...,K_(k-1) determine the class count modulo
2^(2k-popcount(k)). The proof uses the divisibility of
product_(j=0)^(k-1)(h-(2j+1)) at every odd h, whose constant coefficient is odd.
The exponent is sharp within the stated degree-k polynomial-certificate class,
not necessarily among all reconstruction procedures.

At order eight, an isolated decoder accepts exactly eight integers. They give
A000568(8) congruent to6880 modulo32768. A simple globally valid reciprocal
upper polynomial gives A<=31785.04749<32768. Positivity therefore forces6880.
No target count, atlas maximum fibre, finite support list, or Burnside formula
is an input to this decoder. All inequalities are checked exactly over rational
numbers, rather than by sampling candidate fibre sizes.

A separate universal decoder needs only modular moments: if
2k-popcount(k)>=m and sum h=2^m, the residue and1<=A<=2^m uniquely recover A.
The order-eight example uses twelve residues modulo2^22. These are exact
aggregate reconstruction statements, not an efficient way to obtain the moments.

The decisive open step is generating the needed collision moments directly
from order n or a substantially compressed native metagraph description.
The test moments were extracted from an existing atlas. Recovering a held-out
class count from them must not be advertised as a prediction from smaller
orders alone or as a faster replacement for known enumeration formulas.
