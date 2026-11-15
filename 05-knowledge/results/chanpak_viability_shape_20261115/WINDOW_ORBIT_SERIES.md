# A mixed exact repair language and its window generating function

Two arrows preserve the height-zero two-tape carrier and reduce the
count by exactly 2. Write a letter i for a rule
`F(x,A_i z+V_i)=F(A_i x+U_i,z)+2` on all nonnegative x,z.
The physical integer data are:

| Letter | A_i | U_i | V_i |
|---|---:|---:|---:|
| i=12, inherited unit slide in guard coordinates | 4096 | 2787 | 2875 |
| i=36, new primitive return from two opposite suffix swaps | 68719476736 | 46767159475 | 48234217403 |

The [suffix-chart proof](SUFFIX_CHART_SWITCHES.md) gives the
equalities and the 12-bit provenance. The 36-bit return ends in
low 12-bit residues U=2227 and V=4027, neither the short U nor V.
Thus the two codes are **suffix-free on both sides**. This proves
arbitrary mixing of the two repairs without ambiguity, not just
two independent families at fixed width. Both outer addresses are
still arbitrary. Other nonunit moves are deliberately not included
in this chosen rewrite relation.

## Canonical quotients and viable-window series

Whenever y ends in a V letter, move that low block to a U letter
on x and drop the count by 2. Whenever x ends in a U letter,
one can undo the corresponding move. Suffix-freeness gives at most
one edge each way. Each direction removes one block from one
finite address, so all components are finite chains. Starting from
any (x,y), undo all possible moves. The unique upper endpoint has
no trailing U letter on the left. Parse its right address into V
letters, from the low end, until no more remain. Conversely any
outer pair (x,z) with neither of these terminal code constraints,
and any finite word w of repair letters, gives one canonical
component with k=|w| moves and k+1 consecutive odd counts.
This classifies the chosen rewrite relation, not all coincidences
between counts of distinct components.

If Y_w is z with w appended on the right in reverse chronological
order, let H_w=F(x,Y_w) and E_w=H_w−2k be its endpoints.
At the stack level, append a letter V_i to a current y by
```
y' = A_i y+V_i;
p(y')=p(y)+p(V_i)+C_i y,   C_i=A_i−1−rev_i(V_i);
H' = H+[(Qx+S)(A_i−1)+45 C_i] y+(Qx+S)V_i+45 p(V_i).
```
The adjacent correction in the binary split vanishes because both
V words have top bit 1. This gives an exact unbounded recurrence
for endpoints as well as for the number of repair steps.

Let P_k(Z)=sum Z^E_w over all canonical (x,z,w) with |w|=k.
Then the actual window representation series is
```
R(Z) = sum_(k>=0) P_k(Z) (1+Z²+...+Z^(2k)).
```
Its coefficient counts windows containing the target, once per
window. Equivalently it counts physical ordered presentations,
since each chosen chain reaches each of its counts once. It is
an exact regrouping of the full two-address generating function,
not a quotient of all coincident path counts. It also lets one
ask for hits repaired by at least K edges, by keeping k>=K.
Positivity of the relevant coefficients is still a separate
requirement. All coefficients are finite by the positive
address bounds.

## A closed limiting generating function for length

Write rho=1/A_12+1/A_36. Among all physical ordered presentations
with F<=X the limiting probability of being in a chosen chain
of exactly k edges is
```
pi_k = (k+1)(1−rho)^2 rho^k;
sum_(k>=0) pi_k t^k = ((1−rho)/(1−rho*t))².
```
In fact the joint probability generating function of the numbers
of 12-bit and 36-bit moves in the chain is
```
(1−rho)^2 / (1−t_12/A_12−t_36/A_36)^2.
```
For a fixed word w let A_w be the product of its radices. Then
Y_w=A_w z+constant and H_w=Q A_w xz+O_w(x+z+1).
The number of eligible canonical centres through X is
`(1−rho)^2 X log X/(QA_w)+O_w(X)`.
Indeed positive bounded linear terms can be sandwiched using
fixed translations of the coordinate axes. The excluded
suffix constraints have complements of density 1−rho on each
coordinate, periodic modulo A_36. A hyperbola split gives a
leading coefficient of (1−rho)^2 after these two independent
coordinate restrictions: count permitted large coordinates
with their interval density, then sum the permitted inverse
small coordinates, whose harmonic sum has the same density.
The two axes and square intersection cost O_w(X).
Replacing upper by lower endpoint changes the cutoff by at
most 2k, independent of X at fixed w. Thus multiplying by
k+1 counts the chain's representations to leading order.
Normalize by the [full leading term](LOCAL_DENSITY.md), and
sum over words; the weights 1/A_w sum to rho^k.

For completeness, the same law has a uniform large-length
tail bound. The leading appended block gives Y_w>=A_w(z+c0)
for a fixed c0>0 when w is nonempty. Positivity yields
H_w>=c A_w(x+1)(z+1). A window meeting F<=X has H_w<=X+2k,
hence k=O(log X), since A_w>=4096^k. Its outer pairs are
under a divisor bound O(X/A_w). Summing the weighted tail
over all lengths at least K gives
```
O( X log X * sum_(k>=K)(k+1)rho^k ).
```
Thus arbitrary switching provides a real infinite repair
language, but long exact windows are exponentially rarer in
their depth. Above K=(1+epsilon)log_(1/rho) log X their
union still has density zero among targets through X.
An all-target theorem needs a much better coarse selector
than the preceding polynomial-error brackets, or another
independent source of correction.
