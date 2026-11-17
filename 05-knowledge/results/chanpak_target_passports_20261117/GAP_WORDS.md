# Sevenfold Fermat counts cannot come from homogeneous apex gaps

Fix a transitive spine and one apex. A bit **1** means the apex
beats that spine vertex, in increasing spine order. For a finite
word W let L(W) be its ordinary binary integer, R(W) the binary
integer of its reverse-complement, and H(W) the total Hamiltonian
path count of the tournament. Leading zeros and trailing ones
are ordinal padding and leave H unchanged. Write F(a) for H of
the ordinary binary word of nonnegative a, including F(0)=1.
Converse reverses and complements the word, so H(W)=F(L(W))=F(R(W)).

**Scope.** Here a growing gap is a monochromatic block of positive
length affine in the integer parameter k, with strictly positive
integer slope. The finite words between these gaps are fixed.
The target exponent is allowed to be *any* positive-affine
integer Lk+t. Counts are required to match at infinitely many k.
The theorems do not concern arbitrary contact words, variable gap
types, varying intervening words, or more than one apex.

## Exact read and concatenation

A directed odd cycle must use the apex, an outcontact at position
i and an incontact at j>i, and an even subset of intervening vertices.
Two such cycles cannot coexist in a pack because they meet at the
apex. The odd-cycle packing identity therefore gives

    H(W) = 1 + Σ_(i<j, W_i=1,W_j=0) (2^(j-i-1) + [j=i+1]).

If W=UV with both parts nonempty, cross inversions factor:

    H(UV) = H(U)+H(V)-1+L(U)R(V)+[last U=1, first V=0],
    L(UV) = 2^|V| L(U)+L(V).

An empty part contributes no boundary term. In particular F(2t+1)=F(t)
for nonnegative t, since appending a 1 adds no inversions. Another
useful exact identity is

    G(a) := a-F(a)+1-[a odd]
          = number of adjacent 11 pairs
            + Σ_(all pairs of 1 positions i<j) 2^(j-i-1).

This follows by first counting the outgoing weight from a 1 to *all*
later positions, then subtracting those later positions that are also 1.
In particular G≥0, G vanishes on powers of two, and G(7·2^t)=6.

## Arbitrarily many gaps of the same symbol: no plus ray

By converse it suffices to consider zero gaps. Write

    W(k)=P₀ 0^(m₁(k)) P₁ ··· 0^(m_r(k)) P_r,   r≥1.

Discard irrelevant leading zeros until the first fixed prefix contains
a 1; if none of the growing gaps follows a 1, the H read is fixed
and the theorem is immediate. Set a_i=L(P_i), b=R(P_r), and
x_i=2^(m_i(k)). For a single appended stretch the concatenation law gives

    H(U 0^m S) = H(U)+F(R(S))-1-L(U)+[last U=1]
                         + L(U)(R(S)+1)·2^m.

It follows that H(W) is an integer polynomial in the x_i. Its constant
coefficient and unique top multidegree monomial are respectively

    C = F(b) - Σ_(i<r) G(a_i),
    A·x₁···x_r,       A=a₀ (b+1) 2^(|P₁|+···+|P_(r−1)|).

Indeed when a gap variable is set symbolically to 0, the current
left word's integer read becomes its last fixed piece a_i (or 0
when that piece is empty), and its true last bit is the parity of
a_i. Hence the constant terms accrue −G(a_i). In the last stretch
the top monomial comes from the leading contribution to L(U) of the
very first piece; no previously computed H term contains the final
growth variable. These polynomial evaluations are formal and do not
identify x_i=0 with any actual empty gap or tournament.

After substituting x_i=2^(s_i) (2^k)^(r_i) with r_i>0, this top
monomial stays strictly highest, even if a few s_i are negative.
The substituted constant remains C. Infinite equality with the
polynomial target 7(2^(Lk+t)+1) would force the leading coefficient
to have odd part 7 and the constant to be 7. But the leading
coefficient forces b+1 to be a power of two or 7 times a power of
two. In these cases b has, respectively, F(b)=1 or F(b)=F(6)=5
(strip trailing ones). Thus C≤5, a contradiction.

For one gap, when the leading coefficient has odd part 7 the
constants are exactly +5 and −5: if a₀ is a power of two then
G(a₀)=0 and F(b)=5; if a₀=7 times a power of two then G(a₀)=6
and F(b)=1. For multiple gaps, intercept −7 forces ΣG=8 in the
second case or ΣG=12 in the first case. In the second case exactly
one intermediate piece contributes 2 and all others contribute
zero, since 2 is the least positive G. In the first case at most
six pieces contribute positively. These are restrictions on the
minus sign, not a claim that arbitrary many same-sign gaps miss
that sign as well.

## Two gaps of either symbol: no plus or minus ray

Take W=P σ^m R τ^n S and set a=L(P), b=R(S), r=L(R),
r*=ordinary bit reversal of R (preserving its length), p=R(R),
and l=|R|. Here the function R(·) still denotes reverse-complement
integer, so p=2^l−1−r*. Put x=2^m and y=2^n.
By the concatenation rule H=Axy+Bx+Dy+C, with

| (σ,τ) | A | B | D |
|---|---:|---:|---:|
| (0,0) | a·2^l·(b+1) | −a r* | r(b+1) |
| (0,1) | ab·2^l | a(p+1) | b(r+1) |
| (1,0) | (a+1)(b+1)2^l | (a+1)(p−2^l) | (b+1)(r−2^l) |

The (1,1) case is the converse of (0,0). The constant C may
incorporate a fixed-piece adjacency indicator and is not needed
for cancellation when A≠0. The table also holds when R is empty,
by using its length 0 and integer reads zero.

In the mixed cases, if A>0 the lower two coefficients have the
same strict sign and cannot cancel after substitution of any
positive-affine lengths. For (0,1), if a=0 or b=0 then the first
or second growing gap respectively is irrelevant ordinal padding,
reducing to one gap. The other mixed order has A>0 automatically.

For (0,0) assume a>0 after trimming the irrelevant initial gap.
If r=0, then R is all zero and the two gaps combine to one.
If r>0, both lower coefficients are nonzero. Unequal growth rates
put them at different exponents so they cannot cancel. Equal rates
would require a r* and r(b+1) to differ by a power of two only.
If the leading coefficient has odd part 7, the ratio of the odd
parts of a and b+1 is 7 or 1/7. But the odd parts of r and r* are
the two reversals of the same trimmed binary string with first and
last bits 1. Their bit lengths agree, so their ratio is strictly
between 1/2 and 2. This rules out cancellation. The residual
one-gap constants ±5 already exclude ±7.

Run python replay_gap_words.py for exact coefficient checks through all fixed
pieces of length at most three and a few independent literal
tournament counts. These checks validate the formulas; the
all-order restrictions follow from the algebra above, not a
search over gap instances.
