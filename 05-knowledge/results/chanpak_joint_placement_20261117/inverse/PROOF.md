# Exact inversion of the three-contact joint source/shift decoder

Use the three-contact family in `free_second_move/three_contacts/PROOF.md`. Its lower source count is

    N = 81 + 142(u+w) + 256uw + 5(F(u)+F(w)),
    u=G(16x+5)+r,  w=G(16y+4)+(G-1-r),

where G is a power of two at least 4, r is odd with 1<=r<G, x>0 is even, and y>0. Reflection gives N+2, and the three-arc move sends both tracks to N+Delta and N+Delta+2. We solve the **joint** membership query with both N and the signed Delta prescribed, not just displacement. For each G there is at most one physical solution, and over all G there are O(log N) candidates which can be recovered in time polynomial in the input bit lengths. This is not universal placement.

## Eliminate the suffix without forgetting it

Let S=u+w and d=u-w. Put z=x-y and t=rev(r), where the reversal uses exactly log2(G) bits. We have

    S+1 = G(32y+16z+10),      d=16Gz+2r+1.
    Delta = 15y - (16384G^2-5)z - B_r,
    B_r = 4096G^2+1004G+30+2048Gr+40t.

The last expansion follows since the two suffixes are complementary, so their reverse-complements differ by G-1-2t. By eliminating y and z, this displacement equation becomes

    a*d = b*S+c-40t + 5(2r+1)/(32G),
    a = 1024G + 5/(32G),   b = 15/(32G),
    c = -Delta -4096G^2+20G-30 + b(1-10G).

Thus d is the **unique** nearest integer to p(S)=(b*S+c)/a whenever a physical solution exists. Indeed |d-p(S)|<1/25, since 0<=t<G and 1<=r<G; the last term is positive and at most 5/16. The true suffix is not being projected away: Euclidean division of d by 16G recovers z as quotient and 2r+1 as remainder. One can recover the small reversal correction exactly after recovering d. An integer rounded to the right apparent displacement can still fail the exact check.

The remaining inversion is one quadratic. Set

    Phi(S)=64(S^2-p(S)^2)+142S+81.

On physical solutions, |d|<S and 1<=F(v)<=v+1 for v>0. The dual bound is therefore

    |p(S)|<= S+1/25,     |Phi(S)-N| <= 21S        (S>=1).

For the second inequality, replace the exact 64d^2 by 64p(S)^2, costing at most 64(1/25)(2S+1/25), then add no more than 5(S+2) for the F term. These constraints leave at most one **integer** S for each G, even before any congruence test. On the interval cut out by the first inequality, the derivative of Phi is

    128(S-k*p(S))+142,    k=b/a < 1/1000,

and in particular is at least 120S. If S1<=S2 both satisfy the second inequality, integrating gives 60(S2-S1)(S2+S1)<=21(S2+S1), hence S2-S1<=7/20.

## Bit-length algorithm

Because the addresses are positive and the formula contains 256uw, we have S<=N and G<=sqrt(N/256). Enumerate the O(log N) dyadic G. If p(S)=kS+l, the first inequality is just the lower bounds

    S >= (l-1/25)/(1-k),    S >= (-l-1/25)/(1+k).

On their common interval with [1,N], both Phi(S)+21S and Phi(S)-21S increase strictly. Binary-search the first integer S for which Phi(S)+21S>=N, and discard it if Phi(S)-21S>N or the first search has no result. Then recover d by rounding p(S), divide by 16G, require an allowed odd r, and recover y and x from the displayed sum. Verify x,y and both original count equations exactly using the F recurrence. The decoder runs in O(log^2 N) bounded-size rational comparisons plus O(log N) exact word checks. The number of bits of each rational expression is polynomial in log(N+1)+log(|Delta|+2); this claim concerns arithmetic with integer inputs, not a unit-cost oracle for DP on the giant physical graph. If given the upper reflected count or a destination count, subtract 2 or Delta respectively before entering the corresponding branch.

This was suggested by the shape of the argument in Sa'ar Zehavi's [Chabauty–Kim dimension inequality paper](https://arxiv.org/abs/2611.18218), specifically Proposition 2.21 and the surrounding kernel comparison: change the operator/readout while controlling the exceptional part, before bounding the actual residual kernel. The analogy is structural only. Here the discarded bit-reversal term is not a zero module or exact kernel; its normalized magnitude below half a lattice spacing means it cannot give a second integer lift. The exact residual check remains necessary.

Run `python decoder.py` for direct independent prefix count evaluation and joint inversion, including arbitrarily large, highly unbalanced addresses and nearby ghosts. The checks compare the original F-based count at both physical endpoints, not merely the rearranged linear shift.

## A bounded-shift inverse without prescribing the shift

At the native repair scale one can also leave Delta unknown. Fix positive constants c1,c2,c3. Consider only G with c1*G^6<=N<=c2*G^6 and require |Delta|<=c3*G^2. For a physical source, the count bounds imply xy=O(G^4), and the displacement law implies z=O(1) uniformly: for z>=0, y=O(G^2) from x>=y; for z<0, positivity of y in the shift equation bounds |z|. Hence S=Theta(G^3), d=O(G). For the lower bound on S, use N<=64*S^2+147*S+91 and N>=c1*G^6. For the upper bound, bounded z and xy=O(G^4) imply x,y=O(G^2). These bounds have constants depending only on c1,c2,c3.

Now |N-(64*S^2+142*S+81)|<=5(S+2)+O(G^2). For sufficiently large G this is at most 6S, as S is Theta(G^3). The positive real solution to 64*S^2+142*S+81=N determines S uniquely by rounding (the derivative is >=128S on the relevant scale, and the true root differs by <1/2). The finitely many smaller dyadic G can be handled by direct tests. The still unknown integer z is bounded, and is restricted by (S+1)/G == 10-16z mod 64, since x even gives y==z mod 2. For each z, the sum fixes x,y.

Let A0=16x+5, B0=16y+4, and write t=rev(r), s=G-1-r. By the suffix expansion for F, exactly

    F(u)+F(w) =
      F(A0)+F(B0)+(G-1)*A0-2
      -(16z+1)*t+[r<G/2]+F(r)+F(s).

The bracket occurs because A0 is odd and B0 is even. The fluctuating second line is bounded in magnitude by O(G(|z|+1)). Consequently, once S,z are fixed, the target determines d^2 within an interval of length O(G(|z|+1)): subtract the fixed F baseline along with 64*S^2+142*S+81 and use N=64*(S^2-d^2)+142*S+81+5(F(u)+F(w)). For z nonzero, all allowed d=16Gz+2r+1 have |d|>=14G|z|. Thus this squared interval leaves only O(1) possible r, found by two integer square roots and exact formula tests, including the shift-bound check. This proves a polynomial bit-time inverse **with no prescribed Delta** for native-scale, nonzero-z witnesses; the explicit universal shift selector in the three-contact proof already uses 1<=z<=15. At z=0 the same squared interval leaves at most O(sqrt(G)) possibilities and we make no polylogarithmic claim for that branch.

There is an independent absolute-scale obstruction on this native-width scheme: for N in [X,2X], native G is comparable to X^(1/6), so only O(1) dyadic scales occur. S ranges over O(sqrt(X)/G) values because S == -1 mod G, each has O(1) possible z and O(G) possible suffixes. Thus even allowing **all signed shifts in the native-width bound simultaneously**, there are only O(sqrt(X)) source counts at this scale. This is a dimension bound on the exact placement family, not an exclusion for genuinely different templates.
