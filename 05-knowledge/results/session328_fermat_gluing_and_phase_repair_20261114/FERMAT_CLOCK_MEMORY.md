# Fermat seams observe two extra clock coordinates

Let A be a tournament with root v and let (C=A-v). For S⊆C, (P_v(S)) counts paths starting at v and then visiting precisely S, with (P_v(\varnothing)=1). Let (Q_X(z)) be the optional-slot polynomial on X. Define
$$
R^\pm_A(z)=\sum_{S\subseteq C}(\pm1)^{|S|}P_v(S)Q_{C-S}(z).
$$
Gluing the Busch (R_j) seam means add a transitive chain (x_1,\dots,x_j,y) which loses to all old vertices except for the single edge (y\to v). In particular v beats each (x_i). This adds j+1 vertices (j≥1), and outputs a strong tournament whenever v can reach all old vertices. For any A and v the full polynomial is
$$
Q_{\mathrm{out}}(z)=z^{j+1}Q_A(z)
 + z\big[(z+1)^j R_A^+(z)-(z-1)^j R_A^-(z)\big].\tag{1}
$$
At z=1 this gives the established (H_{\mathrm{out}}=H(A)+2^j L) with (L=R_A^+(1)). At z=2 it yields the four-case split also obtained independently by tracking the optional slots directly.

To prove (1), use the odd-cycle packing identity. Each cycle using a new vertex must use the only return edge (y\to v), hence all packings have at most one such nontrivial cycle. If it uses S in the old deleted core and k of the j earlier tail vertices, it has length (|S|+k+2). Its initial old subpath is counted by (P_v(S)); its new portion is ordered and can choose any k-subset. It is retained exactly when (|S|+k) is odd. The cycle weight is (2z); unused j−k tail vertices contribute (z^{j-k}), and old unused vertices contribute (Q_{C-S}(z)). Projection onto odd total cardinality gives (1). With no such cycle all j+1 new vertices are singletons. Reachability suffices for strongness since y can reach v, all other new vertices can reach y, and all old vertices except v have access to y.

## Unsigned interpretation of the negative response

Let (P_{\mathrm{end},v}(S)) instead count paths ending at v with exactly S before it. Then
$$
R^-_A(z)=\sum_{S\subseteq C}P_{\mathrm{end},v}(S)Q_{C-S}(z-1).\tag{2}
$$
This is not a positivity approximation. In the commutative square-zero variables (x_i) for C, let (D_x=\mathrm{diag}(x_i)), G its adjacency, e the all-ones column, and a the indicator column of out-neighbors of v. For a set function identify its square-free generating polynomial with the corresponding element of this ring. Writing (H_C(x)) for the polynomial with coefficient (H(C[S])) at (x_S), and empty coefficient 1, the three path sums are
$$
H_C=1+e^T D_x(I-GD_x)^{-1}e,\quad
\widetilde P_v=1-a^T D_x(I+GD_x)^{-1}e,\quad
P_{\mathrm{end},v}=1+e^T D_x(I-GD_x)^{-1}(e-a).
$$
The tilde introduces ((-1)^{|S|}). Since (G+G^T=ee^T-I), inside the middle expression (I+GD_x) may be replaced by (I-G^TD_x+ee^TD_x): removing a term −(D_x) changes no product after left multiplication by (D_x), because every intervening transition ends in (D_x) and (D_x^2=0). Rank-one inversion now divides (D_x(I-G^TD_x)^{-1}e) by (H_C). Its scalar denominator is (H_C) because reversing all arcs bijects all paths on every subset. After multiplying by (H_C) and transposing the other resolvent, (\widetilde P_v H_C=P_{\mathrm{end},v}). Multiplication is the set convolution; applying the one-slot translation (Q_C(z)=H_C Q_C(z-1)) in this ring proves (2). Hence at z=2 (R_A^-(2)) is exactly the optional two-slot cover with one slot allowed to be empty or end in an in-neighbor of v, with no sign cancellations in the final value. Similarly (R_A^+(2)) is an optional three-slot cover with one slot allowed to be empty or start in an out-neighbor of v.

## Memory obstruction, including the dyadic switches

The triple ((Q_A(2),R_A^+(2),R_A^-(2))) is a complete response for *this profile at a common single seam*, and any three distinct lengths j distinguish unequal triples. This follows since the difference is
$$
2\big[2^j\Delta Q+3^j\Delta R^+-\Delta R^-\big],
$$
and the functions (1,x,x^{\log_2 3}) are linearly independent on any three distinct positive x by strict convexity. When Q already matches, equality at merely two lengths forces both pointed responses to match. A single length can still give an accidental cancellation; the case is not eliminated without a further hypothesis.

For a one-apex source rooted at its apex, with chain of length n and outgoing bit word w, one gets directly
$$
R^+(2)=2^n+\sum_i w_i 2^i 3^{n-1-i},\qquad
R^-(2)=2^n-\sum_i w_i2^i.
$$
The minimum selected index must be a one bit for a rooted nonempty path. The remaining selected vertices are optional and have increasing order, which gives the powers of 3 or 1 in the two respective sums. In particular the negative response alone recovers w at fixed n. Thus distinct same-order one-apex words with equal Q can share the clock under **at most one** common Fermat seam length. At every j with (3^j>2^n-1) they necessarily separate, since equality requires (3^j\Delta R^+=\Delta R^-) and the right side is nonzero of absolute value at most (2^n-1).

For the explicit S326 pair (1^kw_\pm0^k), two old bits change from 0 to 1 and none change back. With plus minus subtraction, (\Delta R^+>0) and (\Delta R^-<0), so they separate under **every** common j. Their powerful exact-budget binary phase spin is an ordinal construction and cannot simply replace an apex gadget in this rooted Fermat seam. This provides an exact structural reason why “same clock” is contextual, rather than a transferable equivalence. A genuine rooted alternative requires matching the new two coordinates, or a deliberately chosen root/unequal seam length that compensates them.

The replay uses independent subset DP for the full optional slots and explicitly counts the rooted marked path part on several literal tournaments; it also checks the growing one-apex separation symbolically via cut evaluation. The arbitrary-j statements follow from (1).

## The fourth coordinate becomes visible on a second seam

For same-root iteration, the complete state is not just the displayed triple. Write (B=Q_C(z)). Appending the chain gives
$$
B'=z^{j+1}B,qquad
R^{\pm\prime}=(z\pm1)^{j+1}R^\pm\mp z^jB.\tag{3}
$$
Indeed the root-started path visits its selected old vertices first and may then use any increasing new subset. The sole excluded configuration is empty old support and new subset ({y}), because the root loses to y. The other chains' leftover vertices factor independently in (Q_{C'-S}). Thus (1),(3) close the full four-coordinate state ((Q_A,R^+,R^-,Q_{A-v})(z)) at any fixed z.

At z=2, if equal-Q sources tie after every first common seam, then (1) forces equality of R⁺ and R⁻. If their deleted-root slots (B) differ, every second common seam separates their Q values: for first length j and second length k the difference is exactly
$$
-2^{j+1}(3^k+1)\Delta B.
$$
So the four-coordinate profile is exactly the equivalence relation for identical Q at every stage of every same-root seam word. One can test just two distinct first lengths and one second step, provided the first Q values match. At z=1 the same recurrence specializes to (L'=2^{j+1}L-H(A-v)).

When v has a unique out-neighbor, the fourth coordinate is already determined by the first three. In the square-zero algebra of (2), write the single-slot v splitting as
$$
Q_A(z)/z = [P_{\mathrm{end},v}(x)P_v(x)Q_C(z-1)]_{x_C}
          =[P_v(-x)P_v(x)Q_C(z)]_{x_C}.
$$
A nonempty root-started path must contain the same out-neighbor, so products of two nonempty path monomials vanish. Consequently (P_v(-x)P_v(x)=P_v(-x)+P_v(x)-1), giving the relation
$$
Q_A(z)/z=R_A^+(z)+R_A^-(z)-Q_{A-v}(z).
$$
This can simplify the search for transferable equal clocks without jettisoning the fourth coordinate for general roots.

## Literal strong equal-clock seeds and every ordered right context

There are strong order-nine seeds at root 8 of outdegree one:

| outgoing masks | H | H(A−v) | R⁺(1) |
|---|---:|---:|---:|
| (24,417,307,486,394,273,439,293,1) | 525 | 263 | 568 |
| (36,273,314,259,361,394,303,351,1) | 561 | 287 | 607 |

They have identical four-state ((Q_A(2),R_A^+(2),R_A^-(2),Q_{A-v}(2))=(13152,7270,3922,4616)). The subset DP in the replay checks all values directly. After any common repeated Fermat seam, they remain strong, have the same order and exact Q, and have distinct H: the initial signed differences ((\Delta H,\Delta h,\Delta L)=(36,24,39)) obey (\Delta L'=2^{j+1}\Delta L-24>\Delta L) and (\Delta H'=\Delta H+2^j\Delta L>\Delta H). In particular the length j=2 gives H=2797 and 2989, both prime to 99 at the new common Q.

In fact the four-state is stable under **every** common ordered right attachment, not just a Fermat seam. Let B be any second rooted tournament sharing the root v, write (D=B-v), and orient all arcs from C to D. Let (E_B^\pm) denote the signed versions of (R_B^\pm) with the marked path *ending* at v instead of starting there. At any fixed z,
$$
\begin{split}
Q_\mathrm{new}={}&Q_A Q_D+Q_C Q_B-zQ_C Q_D\\
 &+z\sum_{\epsilon=\pm1}(R_A^\epsilon-Q_C)(E_B^\epsilon-Q_D),\\
R_\mathrm{new}^\epsilon={}&(R_A^\epsilon-Q_C)Q_D(z+\epsilon)+Q_C R_B^\epsilon,\\
Q_{\mathrm{new}-v}={}&Q_C Q_D.
\end{split}\tag{4}
$$
Terms without arguments are at z. Rooted-start decomposition in the second line first follows a selected C part then any D part if C is nonempty; if C is empty it uses the B-root rule instead. The weighted subset convolution proves the evaluation at (z+\epsilon). For the first line, a mixed cycle must use v. It consists of a nonempty C part, then a nonempty D part, so its length is odd exactly when their combined nonroot cardinality is even. This gives the sum of the two signed terms, unlike the single-(y) parity convention in (1). Internal-only packings are the first three terms, subtracting the twice-counted case v singleton. The last line is ordinal multiplicativity. Equations (1),(3) already showed that four equal coordinates are necessary for equality in all right contexts. Equation (4) proves sufficiency across arbitrary ordered right words at z=2. It does not promise that Hamilton counts stay different under *all* arbitrary right words.

## A cheaper source of exact binary phase bits

Use the literal order-nine pair above at one common seam j. Its H values are
$$
a_j=525+568\,2^j,\qquad b_j=561+607\,2^j,
$$
and its exact shared (q_j=2^{j+1}13152+2(7270\cdot3^j-3922)=\Theta(3^j)). For ℓ>j, using the ratio (r_j=b_j/a_j), one gets the sharp physical bit
$$
v_2(r_j/r_\ell-1)=j,
$$
since the cross determinant (525\cdot607-561\cdot568=27) is odd. Take opposing signs at j and ℓ as in the dyadic construction. Begin at j=2, seeding the two H units mod4 with the order-eight 373/379 exact-Q switch, then use successive j=2,…,e−1 with ℓ=j+1. This gives every H unit mod (2^e) at one exact clock and order with (\log Q=(e^2-4)\log3+O(e)), an improvement in the leading constant over the dyadic-word bit family. To keep H prime to 99 exclude (j\equiv1\pmod{10}) and use up to three copies per stage to bridge the isolated skips, just as before: both H values are units mod3, (b_j) is always an 11-unit, and (a_j=0\pmod{11}) iff (2^j=2\pmod{11}). There is no claim all those j are powers of two. Each bit is built from a strong tied pair that can still be continued at its root, but the multi-bit clock source itself is an ordinal product, not a strong tournament.

## Fermat-only complementary phase rotation at the exact same strong clock

A second literal pair at the same type of root (vertex 8, sole out-neighbor 0) is:

| outgoing masks | H | H(A−v) | R⁺(1) |
|---|---:|---:|---:|
| (168,441,379,368,449,336,387,300,1) | 829 | 349 | 877 |
| (82,316,353,469,484,329,258,359,1) | 841 | 349 | 892 |

Their shared clock state is ((Q,R^+,R^-,B)=(17732,9583,4771,5488)). Strongness and every entry are literal subset-DP controls. The common sequence of seam lengths (j=(1,1,16,2)) uses only powers of two. The exact H outputs are 7,247,174,993 and 7,388,732,915, at common order 33 and exact (Q=1,909,836,390,404,740). Modulo 396 the resulting H values are 377 and 179; they differ by 198 and both are prime to 99. Thus the complementary 99-bank phases can occur at one actual shared clock in **strong** sources, with a ready root for any further ordered seam.

There is an unbounded phase rotation from this pair. On both tracks repeat j=2. Modulo 99 we have (L'=8L-349), (H'=H+4L), and the fixed point (L_*=349/7) is a unit. The multiplicative order of 8 modulo 99 is ten, and (\sum_{i=0}^9 8^i=0\pmod{99}), since 7 is a unit. Each ten-step block resets both L's and moves both H values by the **same unit** (40L_*pmod{99}). Modulo four H never changes, and L is always reset to (-349). Therefore the H difference stays 198 mod 396 after each complete ten-step block, while the common H residue modulo 99 goes through every class. Every residue coprime to 99 consequently has a strong complementary phase pair, and recurs at unbounded larger common clocks by 990-step repetition. This does not keep Q fixed *across* the different rotations; it does keep it exactly fixed *between the two tracks at every common time*.

For contrast, the first 525/561 pair cannot accomplish this particular modulo-9 synchronization in any common Fermat word while H stays prime to 3. In fact its first track ((H,h,L)) is congruent coordinatewise to ((\Delta H,\Delta h,\Delta L)/3pmod3), including at the seed ((525,263,568)\equiv(12,8,13)). The linear seam updates preserve the relation. Thus whenever the H difference is divisible by 9, the first H is divisible by 3. This obstruction is exactly why the second, differently aligned strong seed matters.

The large phase-init word is convenient but not cheapest in this **specific seed-and-seam grammar**. Best-first search on the literal common clock gives:

| allowed lengths j | first successful word by Q | H values | common Q |
|---|---|---|---:|
| all positive integers | (1,2,5) | 812755, 828397 | 1090389448 |
| only powers of two | (1,4,2,1) | 2067209, 2106611 | 4682656948 |

Both preserve the complementary phases; the latter yields H=89 and 287 mod396 and can replace the large initializer in the ten-step rotation. This minimality is a small finite certificate solely for the second rooted pair and common words, not for all tournaments. Along a word the clock strictly increases, and at each node the child clocks increase strictly with j by (1). Thus the replay can discard every child past the published candidate cost and enumerate the remaining tree in best-first order. It pops 92 and 80 words, respectively, through the first success.
