# When a single arc transmits a prescribed star-seam gap

**PROVED.** A five-vertex one-apex pair turns the signed intercept of any right star arm into the exact half-gap of a single arc. A different construction, valid on *every* two-arm star, has minimum nonzero width-independent gap six at intercept one. The distinction is whether the added vertices carry their own tags. Both conclusions concern actual arc reversals.

## A target-first determinant test

Write the left arm's data as \((h,k,C)\), where \(h=H(A)\), \(k\) includes the tag-reading apex, and \(C=h+P_A\). For the right arm write \((g,J,P)\), where \(g=H(B)\), \(P\) counts paths starting tagged and \(J\) counts even paths starting and ending tagged. Each selected path is weighted by the complementary Hamilton count. At width \(m\ge1\) the established star formula is
\[
                      H_m=gk-2Jh+CP2^m .
\]
We take \(P>0\) for a growing ray. An empty left arm has \(h=k=C=1\). [Star product](../chanpak_moonshine_prime_blocks_20261118/STAR_SEAM_PRODUCT.md) · [Signed endpoint identity](../chanpak_moonshine_prime_blocks_20261118/STAR_INVERSE_4096.md)

Suppose only one internal arc of the right arm flips, with data \((g_i,J_i,P_i)\), \(i=0,1\). For the same **signed** difference at two different widths, equality \(P_0=P_1=P\) is necessary and sufficient. To prescribe both counts at one width with the growing term fixed, put their desired residuals \(r_0,r_1\) equal to the two targets minus \(CP2^m\). With \(D=g_0J_1-g_1J_0\ne0\) the only possible left reads are
\[
 k=(r_0J_1-r_1J_0)/D,\qquad
 h=(g_1r_0-g_0r_1)/(2D).
\]
This follows by inverting the two intercept equations. When \(D=0\), keep the two original linear equations. Positive odd integer reads still need one tagged left tournament realizing them together with the requested \(C\).

For a unit seam centered exactly at the growing term, take the unordered residuals \(\{-1,1\}\). Parallel rows are impossible since both \(g_i\) are positive and their requested residuals have opposite signs. Choosing the order that makes \(h\) positive gives
\[
                  h=(g_0+g_1)/(2|D|),\qquad
                  k=(J_0+J_1)/|D|.
\]
In particular, \(k/h\) is the reduced mediant of \(2J_0/g_0\) and \(2J_1/g_1\), and its determinant with each is one in absolute value. Conversely those determinant conditions with opposite signs and a common \(P\) suffice as soon as the actual left arm is supplied. Thus a Farey test selects the *two* required path counts from the right arms before any search on the left. The unsigned slope of that left arm is the remaining coordinate needed to hit the designated absolute center. This complements an earlier Farey test on a different graft. [Near-twin graft](../least_hole_near_twins_20261117/README.md)

Two small certificates illustrate the inverse. In the table a listed arc is reversed, all tags stay fixed and the left arm is empty. Rows are outgoing bitmasks on vertices starting at zero.

| Right rows | Tag mask, arc | \((g_0,J_0),(g_1,J_1)\) | Shared \(P\) | Exact two counts |
|---|---|---|---:|---|
| 34,136,19,133,43,78,31,117 | 164, (5,7) | (431,215), (435,218) | 2230 | \(2230\,2^m\pm1\) |
| 8,5,9,18,7 | 3, (0,1) | (11,4), (13,4) | 32 | \(2^{m+5}+3,\ 2^{m+5}+5\) |

For the first row the fractions \(430/431\) and \(436/435\) select \(1/1\). For the second row, *every* Fermat candidate \(F_b=2^{2^b}+1\), \(b\ge3\), has a one-arc edge at the prescribed counts \(F_b+2,F_b+4\), taking \(m=2^b-5\). This holds independently of primality tests.

## A single-apex flip that reads a signed right arm

In the ascending spine convention with apex vertex zero, take contact words 0101 and 1101; a 1 means that the apex beats the corresponding spine vertex. Give vertices 1 and 3 the star tags in either tournament. The left reads \((h,k,P_A)\) are \((3,29,30)\) and \((5,31,28)\). The flip changes \(h\) and \(k\) by two while changing \(P_A\) by minus two. Thus \(C=33\) stays fixed.

For **every** right arm and **every** width the resulting physical endpoints are exactly
\[
                33P2^m+26g+3s,\quad 33P2^m+26g+5s,
                \qquad s=g-2J.
\]
In particular, a unit seam occurs precisely at \(s=\pm1\), the old signed one-arm selector. The center is then \(33P2^m+26g+4s\). Given a prescribed center \(Z\) and sign \(s\), a necessary congruence is \(26g\equiv Z-4s\pmod{33}\). Given the physical \(g,J,P\), it is sufficient that \((Z-26g-4s)/(33P)\) be an integral power of two at least two. The three inherited right arms with \(g=11,J=5,P=29,37,47\) give actual unit seams centered at \(290+33P2^m\). Their earlier scalar reads alone did not specify a physical edge. [One-arm data](../chanpak_moonshine_prime_blocks_20261118/STAR_INVERSE_4096.md)

## The sharp obstruction for a uniformly tagged block

Append an arbitrary nonempty tournament \(U\) in ordinal order *after* the left arm \(A\) and tag every vertex of \(U\). Write \(u=H(U)\), \(Q(U)=\sum_{S\subseteq U} H(S)H(U-S)\), with \(H(\varnothing)=1\). The new left triple is
\[
                        (hu,ku,CQ(U)).
\]
Indeed the new block is an ordinal suffix of both the deletion tournament and the tournament with its apex, giving the first two products. In the third read, allow the selected path to be empty. If its old part is nonempty and starts tagged, its new part is any first slot of a two-path cover of \(U\). If the old part is empty, the old tournament supplies an unrestricted complementary Hamilton path and every possible new start is tagged. The complementary paths concatenate in ordinal order. This proves the last product including both empty-slot cases.

Consequently **any** old ray \(I+S2^m\) becomes \(uI+Q(U)S2^m\). A physical flip internal to \(U\) yields a width-independent signed gap at two distinct widths of a growing ray exactly when \(Q(U')=Q(U)\). The gap is then \(I(H(U')-H(U))\). The known congruence \(Q(U)\equiv(-1)^{|U|}H(U)\pmod3\) and the oddness of both path counts show that every such gap is divisible by \(6I\) for \(I\ne0\). If \(I=0\), the gap vanishes. [Two-slot congruence](../chanpak_moonshine_prime_blocks_20261118/ORDINAL_RECONTACT_DETERMINANT.md)

The bound is sharp **simultaneously for every star ray**. Use \(U=(6,68,136,131,15,31,61,115)\) and reverse arc \(0,3\). Both have \(Q=5008\), and their Hamilton counts are 313 and 307. Thus every certified intercept-one star, including the earlier bank of 2,842 slope classes, immediately gives an actual edge between
\[
                       5008S2^m+307,\qquad 5008S2^m+313 .
\]
Since \(5008=16\cdot313\), the second count is exactly \(313(1+S2^{m+4})\). This identifies an actual arc leading six below a prescribed multiple, universally on that bank. In the five-vertex construction the tags record changes to endpoint counts, allowing a gap of two. A separate two-port graft uses seven middle vertices and independently branching ends to produce a six-count change. [Star bank](../chanpak_moonshine_prime_blocks_20261118/STAR_INVERSE_4096.md) · [Independent two-flank graft](../chanpak_branching_inverse_20261122/TWO-FLANKS.md)

The finite readings and literal glued controls can be reproduced with `python3 check_star_arms.py` in this directory. The unbounded statements follow from the displayed identities. The existence of a matching physical port for every prescribed odd target remains **OPEN**.

## The upper neighbor at the next perfect exponent

The upper neighbor at exponent 31 also has a short scalar certificate. It predates the [60-vertex single-arc construction](../chanpak_e31_physical_20261122/README.md). Put \(E=2^{30}(2^{31}-1)\). Direct multiplication gives
\[
 E+1=457256543\cdot5042777503=2305843008139952129.
\]
In the physical two-apex join \(H(b,e)=be+F(b-1)F(e-1)-[b,e\text{ both odd}]\), take \((b,e)=(2,159198469)\) and \((2,1996165791)\) for the two factors. The pair \((17,113331563)\) gives \(3200795801=7\cdot457256543\). Here \(F=1+p\) for the binary recurrence in the main note. Ordinal joining therefore supplies both \(E+1\) and \(7(E+1)\). Repeated joins with a tournament of count 49 give \(7^a(E+1)\) for every integer \(a\ge0\). [Two-apex join formula and constructor](../chanpak_prime_frontier_selector_20261118/targeted/PROOF.md) · [Count 49](../least_hole_near_twins_20261117/README.md)

In particular both ordinary neighbors of \(E_p\) occur for each \(p=3,5,7,13,17,19,31\). For \(p\le13\) in this list both lie inside the corrected finite prefix. At 17 and 19 the witnesses share a literal flipped arc. At 31 the lower neighbor follows also from the all-exponent product for \(E_p-1\). The neighboring counts at a general Mersenne-prime exponent remain a sharper question than these finite realizations. [Finite prefix](../../../session203/coverage-recertification/README.md) · [Exact physical flips at 17 and 19](PERFECT-TARGETS.md) · [Perfect-number families](../chanpak_perfect_wieferich_seams_20261122/PROOF.md)
