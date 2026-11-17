# Changing the physical flank removes the fixed-core macroscopic holes

The octal holes in [the fixed-core shadow](PROOF.md) are real but
not intrinsic to the exact +2 seam. Here changing one actual
strong flank removes every fixed relative gap in the joint
source/shift family. The conclusion is relative approximation
with a bounded error in the requested displacement, not
integer placement.

Use [the Farey nontransitive flank](../../chanpak_free_seam_second_move_20261117/extra/PROOF.md):
replace the odd address's bare transitive flank by a transitive
prefix with the same address followed by R_m marked only at 0.
Leave the other address w even and transitive. The Farey determinant
with the even flank is one, so reflection still moves the count
by exactly 2; the fact that R_m is nontransitive is essential
to changing the coefficients. Its counts are

    h_m=2^m+1, P_m=4·2^m+2(3^m-1).
    h=5h_m, P=256P_m, A=147P_m-10h_m, B=142h_m+256(2^(m+1)-1).

The lower read is K_m + h(F(u)+F(w))+A u+B w+Puw, with K_m
constant at fixed m and u odd, w even. This is the previous
Farey read after collecting its coefficients.

Use the same three-contact words, impose z=x-y=1, and pick an
odd r within 1 of G/2. The product and suffix calculations give

    Delta = 3h y -(64P G²-h) z - B_r^(m),
    B_r^(m)=20P G²-4PG(s-r)+2h-4h(C_r-C_s)-4G(A-B).

Here x is required even, so choose y odd. For any requested even
shift delta, choose the odd y nearest to
(delta+(64P G²-h)+B_r^(m))/(3h). The actual shift is even and

    |Delta-delta| <= 3h.

All counts here are physical, and both heights use exactly
three apex reversals on each track. At fixed m, for fixed delta
and G tending to infinity, y/G² tends to 28P/h. Therefore

    sqrt(H)/G³ -> c_m := 448 P^(3/2)/h.

The constants are independent of dyadic level G. For even m,
the phases log_8(c_m) mod 1 are dense. In fact writing m=2n,

    log_8(c_(2n)) = c + n log_8(27/4) + o(1),

because P=(512+o(1))3^m and h=(5+o(1))2^m. The slope is
irrational by unique factorization (a positive power of 27/4
has prime factor 3 and cannot be a power of 8). Irrational
rotation is dense even after removing any finite set of indices,
so the vanishing phase error does not affect density.

**Relative joint placement theorem.** For every epsilon>0 there
is a finite palette of positive even R_m heights and a constant
C_epsilon such that for every sufficiently large positive N,
and any fixed prescribed even delta, the palette produces a
physical four-corner repair with

    |H/N - 1| < epsilon,  |Delta-delta| <= C_epsilon.

The sufficiently-large threshold may also depend on delta.
Choose an epsilon-sized phase net from the dense even-height
phases (with margin for squaring); then let N select the nearest
phase and the appropriate dyadic scale G. The finite palette
makes the asymptotic estimate uniform as G grows, and take
C_epsilon=max 3h over the palette. In particular one can choose
a palette to avoid any proposed fixed relative hole. If the
desired correction must be exact, the separate [arithmetic
gate](../physical/PROOF.md) matters; the rounding here bypasses
that gate only by allowing a bounded discrepancy.

This is a qualitative escape from the original recurring bands,
not a claim of O(N^(1/3)) absolute accuracy with this palette,
nor a joint exact hitting theorem. The number and sizes of
available phases at the much finer resolution N^(-2/3) are
the next issue. The paper's filtration analogy is useful here:
the terminal obstruction in the fixed coefficient layer need
not survive when the coefficient representation itself varies.

## A quantitative version, without a Diophantine approximation black box

For fixed even delta the relative conclusion can be made
uniformly vanishing:

    |H-N| = O(N/log log N),   |Delta-delta| <= exp(O(sqrt(log N))).

This is not a better absolute approximation than the earlier
one-third selector when no repair is required. It gives a
quantitative way to choose an almost specified shift and
a relatively close source **in the same physical family**.

Here is an elementary phase bound. Put alpha=log_8(27/4).
For any integer q>=1 and p nearest to q alpha, the distinct
integers 27^q and 2^(2q+3p) are comparable within a constant,
so |q alpha-p|>= c*27^(-q), by taking their logarithmic ratio.
Pigeonholing the Q+1 fractional parts of multiples of alpha
in Q intervals gives a rational p/q with q<=Q and
|q alpha-p|<=1/Q. Reduce the fraction; its denominator q'
still satisfies that bound, so q'>=c' log Q. The first q'
multiples of the reduced rational form an equally spaced
grid, and the corresponding multiples of alpha are within
1/Q of their grid points. Thus any Q consecutive multiples
of alpha form an O(1/log Q) phase net, including those with
indices between Q and 2Q. Taking n in that range makes the
o(1) coefficient correction exponentially small in Q.

Choose Q comparable to sqrt(log N), fix m=2n from the net,
and choose G by the dyadic normalization. All chosen m
are O(sqrt(log N)), so G=N^(1/6+o(1)) for this palette.
The count estimate used above has relative error
O(1/G+|delta|/(PG²)) uniformly in m: the subleading pieces
in B_r^(m) are O(PG), A+B=O(P), and rounding y costs at
most 1, while the principal y is comparable to PG²/h.
The lower-order terms of H are O(P(u+w+1)). This relative
error is smaller than the phase mesh. The rounding error
in the shift stays O(h)=exp(O(sqrt(log N))). This proves
the displayed quantitative version without asserting a
particular rate for exact-offset phases restricted by gcd.
