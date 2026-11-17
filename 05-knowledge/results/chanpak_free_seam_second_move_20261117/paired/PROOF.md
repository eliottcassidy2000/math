# Common-prefix trades: an exact +512 move and a sharp barrier

Use the [known free reflected seam](../../chanpak_wheel_square_selector_20261117/FREE_SEAM.md).
Write u=Ya+r, w=Ya+s, with Y=2^t, a>0, 0<=r,s<Y,
r odd and s even. This guarantees that swapping u,w still gives
the exact +2 pair. We can now also change (r,s) on the lower side and
carry the same change to the upper side.

## An actual second move independent of the high prefix

The physical moves on the lower side are literal changes of the
apex contacts to the two binary low words. The upper side uses the
reflected contact changes. Let rev be t-bit reversal. The suffix law
shows that the dependence of H on a, apart from the common quadratic
and the common F(a) terms, uses the single linear coefficient

    P Y (r+s) + h(2Y-2-rev(r)-rev(s))

and the parity coefficient h times the number of low words with most
significant bit zero. The fixed core has (h,P)=(5,256). Thus keeping
the column sums of r,s fixed keeps the whole variable part fixed.
The remaining difference is precisely the difference of

    P r s + h(F(r)+F(s)).

For t=3 choose (r,s)=(5,2) or (3,4), i.e. 101/010 or 011/100.
The sum of the two F reads is unchanged. The products are 10 and
12, giving a common exact +512 at *every* a. This makes a physical
four-corner commuting square with one separation +2 and the other
+512. At t=2, trade (3,0) for (1,2): the products grow 0→2 and
the F sum grows by 2, giving +522 at every a. The two displacement
sizes have gcd 2, but these are different tiles; there is no claim
of independent iterable composition or of a Bezout move of length 2.

## Why this entire kind of trade cannot give anything finer

Suppose two low-word pairs of the same width give a constant read
difference for all a (even if a is restricted to one parity).
Equality of the linear coefficients forces both ordinary and
reverse sums equal: |h Delta(rev(r)+rev(s))|<2hY<=PY, so the
integer ordinary-sum difference must be zero. These two sums
hold exactly when the binary-column sums agree. To see the
nontrivial implication, take a first differing low column. Its
difference must be even by equality of direct sums, hence +/-2.
In the reversed-sum equation the earlier columns contribute
nothing and the later terms have magnitude strictly less than
this column's contribution, contradiction. This also ensures
that the parity coefficient is the same for unrestricted a.

Write e_i=r_i-s_i at bit position i, counted from the least
significant bit. Only split columns have e_i=+/-1. To get the
second allocation, reverse the signs on a set A; split columns
not reversed form B. Bit 0 is in B, since the lower orientation
always has r odd and s even. Let D_A=sum_(i in A) e_i 2^i
and define D_B likewise. A nonempty change gives |D_A|>=2 and
|D_B|>=1: these signed sums cannot vanish because the largest
power strictly exceeds the sum of the smaller ones.

The product change is D_A D_B, by subtracting squares of the
old and new differences r-s. The exact F-sum change is

    DeltaF = sum_(i in A, j in B) e_i e_j (2^(|i-j|-1)+[|i-j|=1]).

For the latter identity, F(r)-1 sums the binary prefix and
previous parity each time a zero digit is appended. For i<j
the contribution in the two words is the coefficient
2^(j-i-1)+[j-i=1] times (1-r_i)r_j+(1-s_i)s_j.
The parenthesis equals t_j-(t_i t_j+e_i e_j)/2 for fixed
t_i=r_i+s_i. Reversing just one of its signs changes this by
e_i e_j, as claimed. Leading zeros do not alter F.

For any i in B, the portion D_A^(>i) above i differs from D_A
by at most 2^i-1. It is a multiple of 2^(i+1), and is therefore
zero if 2^i>|D_A|; otherwise its magnitude is at most 2|D_A|.
Summing the powers 2^(-i) over i in B consequently bounds the
non-adjacency kernel from pairs whose lower index is in B by
2|D_A|. Pairs whose lower index is in A similarly cost at most
|D_B|, since 0 is not in A. An adjacency with lower index i in
B and upper index i+1 in A necessarily has D_A^(>i) nonzero
(the largest signed power cannot cancel). There are at most
floor(log2|D_A|)+1<=|D_A| such indices, and at most |D_B| for
the reverse situation. In total

    |DeltaF| <= 3 |D_A| + 2 |D_B| <= 4 |D_A D_B|.

If |D_A D_B|>=4, the full displacement's magnitude is at least
(256-4*5)*4>512. In the remaining case |D_A|=2, |D_B|=1,
bit 1 is in A and B is just bit 0: the higher B contribution
is a multiple of 4 and cannot be nonzero while leaving magnitude 1.
Then DeltaF=e_0(D_A/2+e_1), which has the sign of D_A D_B
or is zero. The magnitude is at least 512, attained above.

Consequently a second fixed-width parallel move is real, but every
nontrivial common-prefix-preserving trade in this complete class
skips at least 512. For any analogous balanced core with P>=8h,
the same proof gives the lower bound 2P. To get exact small errors
one must change the common-prefix architecture, add a new kind of
physical contact, or arrange two trade effects to interact without
forcing this single common-prefix constraint.
