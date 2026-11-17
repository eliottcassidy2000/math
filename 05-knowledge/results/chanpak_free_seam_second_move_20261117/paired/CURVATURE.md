# Why two balanced digit trades do not add independently

The +512 and +522 moves are not a pair of additive odometers.
Already within one fixed binary-column frame there is an
all-width obstruction to creating independent signed controls.

Use the notation of [the paired low-word trade](PROOF.md).
Fix the per-column digit sums and flip signs independently
on two nonempty disjoint sets A,C of split columns, neither
containing bit 0. All four resulting words keep the required
flank parities and therefore still have the exact reflected
+2 move. Write D_A and D_C for the signed binary sums on
the two sets. These are nonzero even integers.

The mixed difference (both flipped minus the two single
flips plus the initial read) is exactly

    -2(P D_A D_C + h J),
    J = sum_(i in A, j in C) e_i e_j
                              (2^(|i-j|-1)+[|i-j|=1]).

The product part follows by writing rs=((r+s)^2-(r-s)^2)/4
and observing that the two flips change r-s by -2D_A and
-2D_C. In the binary F-sum only pairs with one index in each
flip set survive that mixed difference. Each such pair gives
minus twice its single-flip term.

In fact |J|<=2|D_A D_C|. For a fixed i in C, the high
truncation D_A^(>i) (the sum in A of terms above i) vanishes if
2^i>|D_A| and otherwise has magnitude at most 2|D_A|, exactly
as in the main proof. The contribution of 2^(j-i-1) over
all pairs with lower index in C is at most
|D_A| sum_(i>=1) 2^(-i)=|D_A|. For the +1 adjacency term,
nonzero contributions require 2^i<=|D_A|, so there are at
most |D_A|. Do the same for lower indices in A. Then
|J|<=2(|D_A|+|D_C|)<=2|D_A D_C|, since both magnitudes are
at least two.

Thus whenever P>2h, in particular for our core (h,P)=(5,256),
the mixed difference is strictly nonzero and has the sign of
-D_A D_C. Its magnitude is at least 2(P-2h)|D_A D_C|.
This holds at every word size and every common high prefix.
It even holds with distinct, but fixed, high prefixes for the
two words if the designated split columns are below them:
interactions with fixed high columns cancel in a mixed
difference. The result is a scoped constraint, not a universal
prohibition on two physical moves. One could seek deliberately
interacting edits, different core contacts, or a change of reader.
It rules out simply invoking coprimality of observed low-word
steps in this common fixed-column grammar.
