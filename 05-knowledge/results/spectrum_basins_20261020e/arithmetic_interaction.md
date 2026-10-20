# A local interaction restores the Hamiltonian weight

**PROVED:** the full Hamiltonian weight and feedback inventory of golden
SCC words are local statistics of their parity-path occupancy strings.
Their transfer requires four states, with correct empty and short-word
boundaries. **FINITE-EXACT:** all7920 words through order18, independent
transfer generation, literal Hamiltonian counts through order8, and exact
symbolic determinant verification.

This combines our [matching dictionary](arithmetic.md) with the incoming
[Hamiltonian-weight obstruction](../rotation_owners_spectrum_20261020e/arithmetic.md).
That obstruction is correct: independent parity-row weights assign9 to
the R4 block when its actual Hamiltonian weight is5. The new coordinate
below makes the needed interaction explicitly local. Neither earlier
file or receipt is changed by this addition.

## 1. The two local statistics

For a matching on vertices0,...,n-1 with edges `{i,i+2}`, let chi_i be1
exactly when that edge is selected, for `0<=i<n-2`. Distinct selected
edges share a vertex exactly when their start positions differ by2.
Consequently the exact matching constraint is

```text
chi_i chi_(i+2)=0.                                      (1)
```

Define

```text
E=sum_i chi_i,           C=sum_i chi_i chi_(i+1).        (2)
```

Every R4 block in the greedy dictionary selects two edges with consecutive
start positions, so it contributes one to C. Conversely, if chi_i and
chi_(i+1) are both1, then (1) forces chi_(i-2)=chi_(i-1)=0 whenever those
indices exist. No earlier selected edge can consume the beginning of this
four-vertex block. The greedy decoder therefore selects precisely the R4
block starting at i. Three consecutive occupied starts are forbidden by
(1), so these adjacent pairs do not overlap. Thus

```text
C = #R4,           E = #C3 + 2#R4 = c3,
feedback = #C3 + #R4 = E-C,
H = 3^(#C3) 5^(#R4) = 3^(E-2C) 5^C = 3^E(5/9)^C.    (3)
```

The exponent E-2C is nonnegative. The rational-looking last expression
is always an integer because interactions occur in disjoint occupied
pairs. The formula for feedback uses the minimum feedback-arc number:
C3 and R4 each require one deletion, and no directed cycle crosses SCCs.
The formula for H uses ordered SCC multiplication and local counts1,3,5.

A minimal scalar hostile is at order6:

| SCC word | Occupancy string | E | C | H |
|---|---|---:|---:|---:|
|C3 C3|1001|2|0|9|
|R4 1 1|1100|2|1|5|

The matching size loses C. Retaining C restores H and feedback without
needing to decode the entire word. It still does not reconstruct the word:
positions remain necessary for the bijection and marked-edge questions.

## 2. Four-state transfer, including n0, n1, n2

Give each occupied start weight u and each adjacent occupied pair an
additional weight v. The weight of a matching is `u^E v^C`. Store the last
two occupancy bits in the order00,01,10,11, with two virtual zero bits
before the first actual bit. Appending1 is forbidden if the older bit is1;
otherwise it contributes u, and an extra v when the immediately preceding
bit is1. Appending0 has weight1. The row-to-column transfer matrix is

```text
             00  01  10  11
T(u,v) = [   1   u   0   0  ]
         [   0   0   1   uv ]
         [   1   0   0   0  ]
         [   0   0   1   0  ].                         (4)
```

Let Z_n(u,v) sum these weights over order-n golden words. For n>=2 there
are n-2 actual edge positions, so

```text
Z_n = e_00 T^(n-2) 1,       n>=2;
Z_0=Z_1=Z_2=1.                                         (5)
```

The two small initial cases must be stated separately; negative matrix
powers are not used. The empty tournament has one empty combinatorial
word/path by the generating-function convention. The atom weights are
`1,u,u^2 v` at vertex lengths1,3,4, giving

```text
sum_(n>=0) Z_n z^n
 = 1+z+z^2 e_00 (I-zT)^(-1) 1
 = 1/(1-z-u z^3-u^2 v z^4).                             (6)
```

Direct determinant expansion gives the same denominator:
`det(I-zT)=1-z-u z^3-u^2 v z^4`. The equality of the full generating
functions also follows from the bijective atom grammar, so it does not
mistake a characteristic polynomial alone for a generating function.

Set v=1 to recover the independent parity-row triangle enumerator. For
Hamiltonian weight, set `u=3,v=5/9`; then

```text
sum_word H(word) z^|word| = 1/(1-z-3z^3-5z^4),
h_n=h_(n-1)+3h_(n-3)+5h_(n-4),
h_0=h_1=h_2=1,       h_j=0 for j<0.                     (7)
```

This h_n is the **sum** of H over all ordinary golden words of order n.
It is not an individual tournament's H, the number of words, or a
recurrence for unrestricted tournament classes. More generally,
`u=t x,v=x^(-1)` grades c3 by t and feedback by x; expression(3) shows why
the resulting powers are nonnegative despite that substitution.

One can avoid rational transfer entries altogether for H. Give a completed
isolated occupied bit weight3 and a completed adjacent pair weight5:

```text
T_H = [[1,1,0,0], [0,0,3,5], [1,0,0,0], [0,0,1,0]],
terminal vector b=(1,3,1,1)^T,
h_n=e_00 T_H^(n-2) b,              n>=2.                 (7a)
```

The terminal3 pays for a last isolated occupied bit not yet followed by
a zero. Omitting that terminal correction would fail already at n3.
This integer form has the same four states and is checked separately.

## 3. What the actual moves do to these coordinates

The [protected triangle creation](inventory.md) changes the endpoint word
`A 1 1 1 C3` into `A C3 C3`. When A is also a golden word, the matching
description adds exactly one selected edge in the new triple and retains
the guard edge. Therefore E and feedback increase by1, C is unchanged,
and endpoint H is multiplied by3. The all-order theorem permits arbitrary
source A, for which SCC multiplication still gives that same H factor;
the occupancy description is restricted to golden A.

The inherited R4 dissolution `R4 -> C3 1` instead changes `(E,C)` from
`(2,1)` to `(1,0)` within that block. It preserves endpoint feedback and
multiplies endpoint H by3/5. Swapping adjacent singleton/C3 atoms preserves
both E and C. Thus H is not a black-path invariant even when both endpoint
tournaments belong to the golden language. Intermediate states of these
black paths need not lie in this language, so no monotonicity or transfer
rule for their intermediate states is asserted.

The [regular triangle-reversal theorem](regular.md) is a different
operation: it preserves every vertex degree and therefore c3, while
allowing H175 to change to171 at order7. General regular tournaments
are not encoded by this golden occupancy model. The comparison concerns
the common observable c3, not an identification of the move systems or
of their state spaces.

## 4. The incoming shear address equals the half-plane matrix

The incoming arithmetic note uses positive shears

```text
L=[[1,1],[0,1]],       R=[[1,0],[1,1]].
```

Our two-step Fibonacci matrix is exactly

```text
M=L R=[[2,1],[1,1]],
(q,p) -> (2q+p,q+p).                                   (8)
```

The rightmost shear acts first. Equivalently, if S swaps the coordinates,
one Fibonacci step is `LS`, and `(LS)^2=LR`. Thus the upper-half-plane
action in our arithmetic note is a specific alternating-shear address in
the incoming rational-pair tree. It advances two Fibonacci steps and
four tournament vertices. It retains the marked coordinate basis and
operation order; it is not a move of the fixed-order black metagraph.

## Reproduction

```text
python3 05-knowledge/results/spectrum_basins_20261020e/arithmetic_interaction.py
python3 -O 05-knowledge/results/spectrum_basins_20261020e/arithmetic_interaction.py
```

The script independently generates words and four-state strings, compares
every `(E,C)` coefficient through order18, checks literal Hamiltonian
counts on all words through order8, expands the4-by-4 determinant exactly
using integer polynomial dictionaries, and checks the shear address.
No main research checker is imported. All validity checks use explicit
exceptions. The JSON contains source hash, full coefficient distributions,
and example strings; normal and optimized receipts agree, with JSON SHA256
`7d665224c7f8aa17209e3b726123655f9338ce0f3620b3b4869b7ae3888d490b`.
