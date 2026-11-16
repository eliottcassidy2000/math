# Whole-tape determinants with a prescribed gap as large as X^(2/5)

**Theorem.** For every sufficiently large X and every even integer d
with |d|≤X^(2/5), there are Ω(X^η) different nonzero whole-tape
determinants of literal strong two-apex pairs with both Hamilton
counts ≤X and signed difference d by reversing the mutual apex arc,
where η=(log_3(2)−2/5)/4≈0.05773. The implied constants are uniform
in d. This is a determinant count, not an assertion that Ω(X^η)
different Hamilton counts occur.

Use the [fixed patch](DETERMINANTS.md) p=0001120. Its own signed
gap Δ=H(A→B)−H(B→A) is 8. Its marked forward address difference
L_b−L_a=10 and reverse-complement difference R_b−R_a=8.
On the left place a common-row word u of forward address T; on
the right place a common-row word v of reverse-complement address S.
The local endpoint corrections vanish (the patch endpoints are
common), so

    Δ(u p v)=8+10S−8T.

Choose S0∈{0,1,2,3} so that 2S0≡−|d| mod8, and
T0=(|d|+8+10S0)/8. For an integer t let

    T=5t+T0,    S=4t+S0.

Then Δ=−|d| exactly; swap the two apex names if d>0.
This sign convention does not affect the determinant count.

## Distinct determinant choices without paying twice for the gap

Choose L≥1 and take
t=3·16^L+Σ_(i=0)^(L−1) ε_i16^i, with ε_i∈{0,1,2,3}.
The width m_R=4(L+1) of S is fixed, and S starts in 11:
12·16^L≤S<13·16^L. Also T≥15·16^L, so T has
width m_L≥m_R and u starts in letter 3. The last site of v is
letter 0. Both orientations are strong: both apices reach the
first spine site and can be reached from the last.

At a lower hexadecimal digit, let a be the digit of T0 plus
the incoming carry, so 0≤a≤16. Pick two ε from {0,1,2,3} with
the same outgoing carry from a+5ε. In each case the following
weight differs by at least 54 for some such pair:

    -10 rev_4((a+5ε) mod16) + 8 rev_4(4ε).

The finite integer check for all 17 possible a is given in
check_asymmetric.py (its smallest best separation is 54).
The S0 in the lowest S digit adds just a constant to all of its
reversal weights, since S0<4. Choosing the two options at each
digit keeps every subsequent input carry fixed. Thus all 2^L
choices have the same high part of T beyond these digits and
the same widths. In particular
m_L≤max(m_R,ceil(log_2(T0+1)))+1.

With r_m the reverse-complement address, the **whole-tape**
signed determinant D=L_b R_a−L_a R_b satisfies

    D / 2^(m_L+m_R)
      = 2^7(10S−8T) + 10 R_a(p) − 8 L_a(p)
        + 10 r_(m_L)(T)/2^m_L − 8 r_(m_R)(S)/2^m_R.

The first three terms are constant, and all high digits outside
the variable blocks are fixed. The last terms have base-16
weights with the indicated differences at each variable nibble.
A leading difference of at least 54 dominates all later ones:
the entire raw digit range is contained in [−150,24], and
174/(16−1)<54. Distinct ε choices hence give distinct D, with
at most one zero. This is a statement about the whole two-row
address determinant, not just a large canceling cross term at a cut.

## Budget

A spine of n=m_L+m_R+7 sites has H≤2·3^n: select one
of two apex orders, then one of three increasing runs for each
spine vertex. If B=floor(log_3(X/2))−7 and
Q=ceil(log_2(T0+1)), choose m_R a multiple of four at most
min((B−2)/2,B−Q−2), losing at most four in rounding.
Then n≤log_3(X/2), and m_R≥8 for large X uniformly under
|d|≤X^(2/5). The number of nonzero determinants is at least
2^(m_R/4−1)−1, or

    Ω(min(X^((log_3 2)/8), X^((log_3 2)/4)/(2+|d|)^(1/4))).

This gives the displayed η. The bound does not assert H≈X:
when a large requested gap dictates a long left address, its
first two distinct winning contacts need not be near the boundary.
Therefore it is not a claim that the new supply alone fills the
bounded-two-port capacity gap, or that an earlier approximate
baseline can be reused unchanged. The missing step is a joint
decoder for H and the determinant at the prescribed target.
