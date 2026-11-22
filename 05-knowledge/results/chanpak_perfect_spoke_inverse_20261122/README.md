# Which unequal flanks can place a marked spoke at a perfect center?

**PROVED.** Fix a tagged tournament, a vertex in it, and two desired directed Hamiltonian-path counts. Add transitive binary flanks of arbitrary positive lengths on either side and an apex. When the edited arc joins that apex to the chosen vertex, there is at most one address pair per parity class yielding the prescribed ordered counts. For the centers \(E_p=2^{p-1}(2^p-1)\), rationally aligned address ratios give an all-exponent dichotomy below. The successful exponent-31 construction reverses an arc *inside* a marked core and thereby retains an extra digit freedom. [Exponent-31 realization](../chanpak_e31_physical_20261122/README.md)

## An inverse with at most four exact tests

Take a nonempty tournament \(D\), its tag set \(T\), and an initially untagged vertex \(u\). The apex beats exactly the tags in \(D\). The left flank has ordinary binary address \(x>0\); the right contact word is the reversed complemented binary representation of \(y>0\). Holding both flanks fixed, add \(u\) to the tags. This reverses one arc. The completed tournaments are strong, since the first left vertex is beaten by the apex, the last right vertex beats the apex, and all arcs between earlier and later blocks go forward. The same inverse works for a lower endpoint on either orientation by trying both orders of the requested counts. [Unrestricted marked flank law](../chanpak_decimal_core_selection_20261117/FLANKS.md) · [Earlier single-apex inverse](../chanpak_mixed_defects_20261115/SINGLE_SITE.md)

Let \(h=H(D)\), \(Q=\sum_S H(S)H(D-S)\), and let \(A,B\) count nonempty subset paths starting at \(u\), respectively ending at \(u\), weighted by the complementary path count. Both are positive. For fixed parities \(i=x\bmod2,j=y\bmod2\) the old height is
\[
 H_0=K+Qxy+Rx+Ly+h(f(x)+f(y)),
\]
where \(L\) is the tagged-start mass and \(R\) the untagged-end mass. Here \(f(0)=0,\ f(2v+1)=f(v),\ f(2v)=f(v)+v+[v\text{ odd}]\); thus \(0\le f(v)\le v\). The parity-adjusted constant \(K\) includes both signed endpoint corrections in the cited flank law. The new tag raises \(L\) by \(A\), lowers \(R\) by \(B\), and replaces \(K\) by \(K'\). Therefore
\[
                    H_1-H_0=K'-K+Ay-Bx .
\]
For a prescribed signed change \(d\), take the positive integer solutions of this equation in the chosen parities. They are empty or form one ray
\[
 (x,y)=(x_0,y_0)+t(2A/g,2B/g),\qquad t=0,1,\ldots,\quad g=\gcd(A,B).
\]
Indeed after subtracting the two parity bits this is the ordinary two-variable linear Diophantine equation, and both increments are positive.

Set \(G(t)=K+Qx(t)y(t)+Rx(t)+Ly(t)\). The exact height lies in \([G(t),G(t)+h(x(t)+y(t))]\). These intervals are strictly disjoint in their natural order. In fact, writing \(d_x,d_y\) for the increments,
\[
 G(t+1)-G(t)\ge Q(d_x y+d_y x+d_xd_y)> h(x+y)
\]
since \(Q\ge2h\) by its two empty-slot terms. For a desired old height \(N\), choose the last \(t\) with \(G(t)\le N\) by one integer square root. Only this candidate needs the binary recurrence. This proof neither assumes monotonicity of \(f\) nor imposes a bound on the flank lengths.

## A rational-slope dichotomy without a preselected word template

Fix the entire marked core, parity class, signed change and an integer offset \(c\). Set
\[
 \lambda_x=\sqrt{A/(2QB)},\qquad \lambda_y=\sqrt{B/(2QA)}.
\]
Suppose these leading coefficients are rational (one is rational iff the other is). If either is nondyadic, the old count can equal \(E_p+c\) for only finitely many exponents, even if an entirely new word pair is chosen at each exponent. If both are dyadic, infinitely many hits are possible in this class **if and only if** some fixed integer offsets \(a,b\) give
\[
              x=\lambda_x 2^p+a,\qquad y=\lambda_y 2^p+b
\]
realizing the same two target counts for **every sufficiently large integer** \(p\). This provides a finite coefficient test in that case. The irrational-slope case is outside this classification.

**Proof.** Put \(Y=2^p\) and \(\rho=B/A\). The gap equation gives \(y=\rho x+\eta\) with fixed rational \(\eta\). The height equation and the bounds on \(f\) give \(Q\rho x^2=Y^2/2+O(x+Y+1)\), hence \(x=\lambda_xY+O(1)\) and \(y=\lambda_yY+O(1)\) by division by \(x+\lambda_xY\). The bounds are effective. For example, with
\[
 C_1=|Q\eta|+R+\rho L+h(1+\rho),\quad
 C_0=|K|+(L+h)|\eta|+|c|,
\]
we have \(|x-\lambda_xY|\le (C_1+(1/2+C_0)/\lambda_x)/(Q\rho)\) at every hit; use \(Y\ge1\). The affine relation then bounds the other offset.

On an infinite set of exponents we may fix both offsets from the floors and fix the exponent modulo the periods of the two rational binary expansions. A nondyadic rational has a recurring binary block containing both digits. Its rounded address, after a fixed integer adjustment, has a fixed prefix, a growing repetition of this block, and a fixed suffix. Indeed carries from that adjustment are confined to a bounded suffix because both digits recur. The positive-drift lemma shows that \(f\) of this address is affine in \(2^p\) plus a strictly positive coefficient times \(p\). For a dyadic coefficient, a bounded adjustment instead leaves a growing run of zeros or ones with zero drift. Both addresses themselves are affine in \(2^p\) on this progression. Thus the linear term in \(p\) arising from even one nondyadic side cannot cancel: the two \(f\) terms have the same positive weight \(h\), and the remaining terms are polynomials in \(2^p\). Such an expression cannot equal \(E_p+c\) infinitely often, by successive division by its dominant exponential and then by \(p\). [Positive drift for a recurring mixed word](../chanpak_exponent31_seams_20261122/README.md)

If both coefficients are dyadic, their product \(1/(2Q)\) forces \(Q\) to be a power of two and each coefficient to be an integral power of two. Enumerate the bounded integer offsets just obtained. For each offset, every sufficiently long address is a fixed suffix following a growing run of zeros or ones, with a fixed leading digit. Its binary read is exactly affine in \(2^p\): append the suffix using the binary split law or use \(f(2^s)=2^s\) for \(s\ge1\) and \(f(2^s-1)=0\). Its parity is constant. Consequently the height identity is a polynomial identity of degree at most two in \(2^p\); infinitely many solutions force coefficient equality. The gap equation is affine. These identities also suffice. Earlier work excludes dyadic-polynomial prescriptions in another fixed unit chart. Here the linear equation imposed by a spoke forces the address approximation for *all choices* of flanks. [Dyadic perfect targets](../chanpak_absolute_liouville_leech_20261122/DYADIC-PERFECT.md)

**FINITE-EXACT.** Among all stored isomorphism classes of deletions through eight vertices, the only marked pairs \((D,u)\) with both leading coefficients dyadic have transitive \(D\). The numbers of classes at orders 1–8 are 1, 1, 2, 4, 12, 56, 456, 6880. At order \(n\) exactly \(n\) choices of representative and mark pass, all in its transitive class. Compute \(Q,A,B\) directly and put \(a=A/g,b=B/g\); the criterion is that \(Q,a,b\) be powers of two and \(2Qab\) be square. On transitive deletions the old sharp one-site gap excludes the exact perfect center at every \(p\ge3\). Hence no rationally aligned marked-spoke chart through eight vertices straddles infinitely many perfect-form centers, even if the addresses can vary freely. This says nothing about the order of an arbitrary tournament witness. [Sharp one-site gap](../chanpak_exponent31_seams_20261122/README.md)

## A common divisor visible from either end of a flipped arc

More generally let the deletion be \(L_0\to K_0\to R_0\) in ordinal order, with \(u\in K_0\). The outer blocks may be empty, with empty path count one. For **every** choice of the apex tags other than \(u\), the half-change from adding \(u\) is divisible by \(\gcd(H(L_0),H(R_0))\). In particular the gcd must be one for a physical unit flip. Deleting the other endpoint of the same arc gives a second necessary condition. Arbitrary exterior transitive flanks preserve the two outer counts.

To see this, write the apex count as \(h+2(e-J)\), where \(e\) counts even subset paths starting tagged and \(J\) counts those both starting and ending tagged, with complementary path weights. All new terms start or end at \(u\). If they start there, they avoid \(L_0\), so their complementary counts have \(H(L_0)\) as an ordinal factor. If they end there, they avoid \(R_0\) giving the other factor. [Even double-tag identity](../chanpak_moonshine_prime_blocks_20261118/STAR_INVERSE_4096.md)

For singleton \(K_0=\{u\}\) the complete set of possible changes is explicit. Write \(O_v^+(L_0)\) for the mass of odd subset paths starting at \(v\) and \(O_w^-(R_0)\) for those ending at \(w\). Then
\[
 \tfrac12(H_1-H_0)=H(L_0)\sum_{w\in R_0\setminus T} O_w^-(R_0)
                   -H(R_0)\sum_{v\in L_0\cap T} O_v^+(L_0).
\]
Deleting the initial or final \(u\) from each new even path gives exactly the indicated odd path and leaves the complementary weight unchanged up to its outer factor. The two selections are independent. Hence the gcd of *all* half-changes is exactly the gcd of the displayed individual summands, with an empty list contributing zero. An individual target still has to be a difference of the two bounded subset sums. For instance take \(L_0\) to be the directed triangle and \(R_0\) to be the converse of the rows \((0,1,11,35,15,23)\). Their ordinary counts are 3 and 5, and their respective odd endpoint weights are \((2,2,2)\) and \((5,5,12,10,18,14)\). The gcd of the individual summands is one; their smallest nonzero absolute difference is five.

The first *nontransitive* dyadically aligned example is at order nine. With \(U=(2,8,3,5)\), take \(D=U\to u\to U\). Here \(H(U)=5,Q(U)=32\), whence \(A=B=160,Q(D)=2048\) by ordinal multiplication of the two-slot read. The corresponding leading coefficients are both \(1/64\). Every half-change at this marked vertex is divisible by five, so this alignment produces no unit seam at any height. As a contrasting literal example of the general inverse, on rows \((4,65,50,39,75,147,45,95)\) with tag mask 244, adding vertex 3 with flanks \((25441035,21176353)\) gives a unit seam centered at \(E_{31}-78364\). The exact exponent-31 result uses the independent receiver construction linked above.

Run `python3 spokes.py --catalogue` in this directory for the dyadic census and the physical controls. The inverse itself accepts arbitrary Python integers. Rows are outgoing bitmasks and vertex zero is the low bit.
