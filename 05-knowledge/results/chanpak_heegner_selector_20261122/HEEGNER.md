# Exact target selection on extensions of the Heegner-prime tournament

**PROVED with finite-exact frame data.** Keep the quadratic-residue orientations among the seven odd Heegner primes. Add a singleton with arbitrary contacts, then a vertex with arbitrary contacts, and replace this last vertex by a transitive tournament of any positive width. This gives 28,922 distinct growing Hamilton-path count functions. Every prescribed positive integer has an exact inverse test on their union. The five large Heegner discriminant primes together with 37 suffice to recover every admissible width by modular lifting and an absolute check. Seven different growing count functions have a fixed arc reversal adding two at every width.

## A five-bit contact code

The nine positive squarefree radicands of imaginary quadratic fields of class number one are \(1,2,3,7,11,19,43,67,163\). Their negative fundamental discriminants use \(-4,-8\) for the first two. For the five odd Euler lucky primes, the identities \(d=4p-1=8b+3\) align the three lists in the question. The polynomial \(x^2+x+p=((2x+1)^2+d)/4\) is the norm of \(x+(1+\sqrt{-d})/2\). E.g. \(p=41,b=20,d=163\). Its value at \(x=p-1\) is \(p^2\). At the CM point with discriminant \(-163\) the value of \(j\) is exactly \(-640320^3\); the constant coefficient 744 of the Laurent expansion then gives the integer near \(e^{\pi\sqrt{163}}\). [Stark, class number one](https://doi.org/10.1073/pnas.57.2.216) · [Green, Theorem 58](https://people.maths.ox.ac.uk/greenbj/papers/ramanujanconstant.pdf) · [Earlier lucky-prime bridge, HYP-2225](../../hypotheses/HYP-2225-prime-heegner-tournament-boundary.md) · [Earlier CM and splitting calculations](../ramanujan_heegner7_trivial_cycle_20261006.md)

Put \(D=(11,19,43,67,163)\). For each \(d\) in this list and integer \(u\) coprime to the product, use the sign \(\chi_d(u)=(u/d)\). By reciprocity this is the splitting character of discriminant \(-d\). On the [five known Fermat primes](../chanpak_e61_fibonacci_target_20261122/README.md), the signs are

| \(u\) | 11 | 19 | 43 | 67 | 163 |
|---:|:---:|:---:|:---:|:---:|:---:|
|3|+|−|−|−|−|
|5|+|+|−|−|−|
|17|−|+|+|+|−|
|257|+|−|−|+|−|
|65537|−|+|−|−|−|

**PROVED.** Each of the 32 possible sign words is obtained from a unique positive divisor of \(2^{32}-1\). The five rows form a basis over \(\mathbf F_2\) with minus encoded by one, and their prime product is \(2^{32}-1\). Thus the desired five orientations specify the squarefree representative through binary linear algebra. By comparison, the odd-part and powers of two in \(1,2,5,8,20\) generate just four sign words: only the classes of 2 and 5 occur independently. The five Platonic vertex counts \(4,8,6,20,12\) generate eight words, since they supply the independent classes of 2, 3 and 5. This character rank measures the actual contact information each pentad can supply. [Earlier Fermat and Platonic constructions](../chanpak_fermat_solids_seams_20261122/README.md)

In the Legendre tournament on primes congruent to 3 modulo 4, orient \(q\to p\) when \(q\) is a square modulo \(p\). For a new prime, take the representative divisor above modulo \(\prod_{d\in D}d\), add the chosen nonzero square or nonsquare residues at 3, 7, and any further outside primes, and impose \(3\bmod4\). The Chinese remainder theorem and Dirichlet's theorem give infinitely many primes meeting the specification. Iterating realizes any requested contact pattern and every transitive replacement with actual distinct prime labels. [THM-5125, Legendre universality and the seven-prime tournament](../../../01-canon/theorems/THM-5125-residue-tournament-model-heegner-domination-and-cycle-gas-identities.md)

There are two different nine-vertex questions here. If all nine fields themselves are retained, adjoin the fields of discriminants \(-4,-8\) to the seven-prime frame. Orient each new field toward an old prime when that prime is inert there; for two old primes this is exactly the preceding orientation. The \(-4\) vertex beats all seven, the \(-8\) vertex beats only 7. The pair above 2 ramifies in both fields, leaving their mutual arc unspecified. The two completions have heights 161 and 983. If instead the two new vertices are unrestricted primes as in the theorem, there are \(36-\binom72=15\) free arcs. On the 32,768 labeled extensions the heights have 1,214 distinct values, including **every odd height from 405 through 2297**. Exactly 441 different even midpoints admit a difference-two reversal of a free arc while the frame stays fixed. These are counts of values, without isomorphism multiplicity.

The upper endpoint of that consecutive interval is the Euler norm \(47^2+47+41=2297\). The next odd integer \(2299=11^2\cdot19\) is absent from these nine-vertex extensions. More generally write \(f_p(x)=x^2+x+p\) for a row \(d=4p-1\). Given an odd integer \(m\ge1\), there is an integer \(x\) with \(m\mid f_p(x)\) exactly when \(d^2\nmid m\) and every prime \(q\mid m\) other than \(d\) has \(\chi_d(q)=+1\). Indeed \((2x+1)^2\equiv-d\pmod q\) has two simple roots at split primes, each lifting uniquely to every prime power. At \(d\) it has one root modulo \(d\) and none modulo \(d^2\). The Chinese remainder theorem joins the choices. For \(q\equiv3\pmod4\) distinct from \(d\), this divisibility condition is precisely the directed contact \(q\to d\) in the Legendre tournament. Thus the five bits selected by the Fermat divisors also encode the local permissions for prime divisors of the five Euler norms.

In particular the Fourier coefficient **744** is a literal nine-vertex midpoint. In order \((3,7,11,19,43,67,163,2683,2311)\), the Legendre tournament has 743 Hamilton paths. Replacing the last prime by 751 gives 745. The only arc that changes joins the last two vertices. The singleton 2683 has outmask 27 toward the frame, and 2311 and 751 have respective outmasks 167 and 39 toward the first eight vertices.

## Where the equal Heegner and fork heights diverge

The six frame primes other than 7 give height 29; adjoining 7 gives 71. The [physical tripling fork](../chanpak_ternary_buckets_20261122/DIRECT-FORK.md) also starts at heights 29 and 71. Replace 7 by any nonempty tournament \(M\) with the same external incidences. In the physical construction use its growing port filled by one ordinal singleton followed by \(M\). If \(P_j(M)\) counts ordered covers of all its vertices by \(j\) nonempty directed paths, then
\[
 H_{\rm Heegner}(M)-H_{\rm fork\ source}(M)=9P_2(M)+10P_3(M).
\]
Indeed 7 has two out-neighbors in the frame, so at most three visits to the module occur. The Heegner totals at transitive widths 1, 2, 3 are 71, 173, 437. The triangular path-cover table then gives \(71P_1+51P_2+10P_3\) for every \(M\). The physical source gives \(21Q_2-13P_1\) on its filled port, where \(Q_2=2P_1+P_2\) allows empty slots. Appending the singleton doubles \(Q_2\) by choosing its slot and leaves \(P_1\) unchanged. Subtraction proves the formula. At width zero the deletion has height 29 on both sides. For nonempty \(M\), equality holds exactly at one vertex, since splitting a Hamilton path gives a two-path cover whenever \(|M|\ge2\).

At transitive width \(s\ge1\), the two heights are \(10\cdot3^s+21\cdot2^s-1\) and \(42\cdot2^s-13\). Their difference counts nine copies of the length-\(s\) words on three colors that use each of two specified colors, plus one copy of the words that use all three. The nonempty-slot exponential extrapolates to 30 at zero, whereas deletion gives 29; these two evaluations have different domains. The first disagreement is 173 versus 155.

## The integer inverse on the whole growing bank

Write \(A\) for the mask of the singleton toward the seven primes in increasing order and \(B\) for the module outmask toward those eight vertices. Delete all module vertices from an expanded Hamilton path, marking the gaps they occupied. Each gap contains one nonempty path through the module, and between gaps the outside vertices must have the requisite incoming and outgoing contacts. Every gap except a possible first one requires a different incoming contact, and every gap except a possible last one needs a different outgoing contact. Among eight outside vertices one part has size at most four. Thus at most five gaps occur. Let \(b_i\) count the legal outside orders with \(i\) occupied gaps. For a transitive block of width \(s\), a partition into \(i\) ordered nonempty slots determines each internal path uniquely. Consequently
\[
  H(s)=\sum_{i=1}^5 b_i i!{s\brace i}=\sum_{j=1}^5 a_j j^s,
  \qquad a_j=\sum_{i\ge j}(-1)^{i-j}\binom ij b_i .
\]
The C++ recurrence appends each unused outside vertex directly, or through an occupied gap when both contacts permit this. Its integer coefficients are the \(b_i\), including \(b_0\) for the deleted case only. The sequences \(1^s,\ldots,5^s\) are independent by Vandermonde, so distinct coefficient vectors represent distinct functions. There are 29,006 functions, of which 28,922 grow. On a growing profile each target \(T\) can be attained only at \(s\le\lfloor\log_2 T\rfloor+5\). Indeed for an occupied count \(i\ge2\) and \(s\ge i\), first assign separate vertices to all \(i\) slots and distribute remaining vertices in at least two ways apiece. Smaller widths satisfy the bound directly.

For a prime \(p>5\), initial phase \(r\bmod(p-1)\), and \(F(s)=\sum a_j j^s\), define
\[
  L_p(r)=\sum_{j=1}^5 a_j j^r\frac{j^{p-1}-1}{p}\pmod p.
\]
To lift a solution \(s\equiv r\pmod{p-1}\) with \(F(s)\equiv T\pmod{p^k}\), replace \(s\) by \(s+t(p-1)p^{k-1}\). The new height modulo \(p^{k+1}\) is \(F(s)+t p^kL_p(r)\) by the binomial theorem. A nonzero \(L_p(r)\) gives exactly one next digit. As soon as the exponent period \((p-1)p^{k-1}\) exceeds the absolute width bound, only one positive integer exponent in that phase can remain. Its exact height decides the target. This takes \(O(\log\log T)\) digits per surviving phase.

**FINITE-EXACT phase certificate.** Exactly 26,209 growing profiles admit at least one \(p\in D\) that works on every phase. Permitting \(p\) to depend on the phase leaves only eleven singular classes modulo \(\operatorname{lcm}_{p\in D}(p-1)=62370\), on nine of the 28,922 profiles. Adding 37 resolves them all, even after all lifts to the joint phase modulus. Here the additional Fermat period is precisely the 36 arcs of a nine-vertex tournament. The certificate tests the prime from its derivative; a cardinality comparison alone would give no inverse. Thus iterating the inverse over these six primes is complete across arbitrary widths. To keep the whole prime set Heegner instead, 7 resolves ten of the eleven. The residual coefficient vector is \((285,-638,-130,612,48)\) at \(s\equiv34064\pmod{62370}\). On this even phase the total derivative from the unit bases \(2,4,5\) equals 1 at 3. At precision \(3^k\), check all short widths \(s<k\) directly, and omit the term \(-130\cdot3^s\) for the remaining widths. This last observation also gives an inverse using the seven odd Heegner primes.

## Seven physical single-arc families and an exact perfect-number obstruction

Keep the module fixed across a reversal, its seven original prime contacts fixed, and allow the new singleton to reverse one of its seven contacts with the frame. On **every** positive width, the difference is two exactly when the two nonempty slot vectors differ by two in their first component only. The complete bank has eleven distinct such lower profiles, seven growing. The following rows give one labeled realization per growing profile. The flipped pair joins the singleton to the displayed prime. The coefficients refer to the lower height, the midpoint is one larger.

| \(A\) | \(B\) | Prime | Nonzero-prefix slot coefficients \((b_1,b_2,\ldots)\) |
|---:|---:|---:|---|
|11|1|67|309, 134|
|11|5|67|411, 238|
|55|63|163|497, 351, 72|
|69|223|11|627, 297|
|73|3|19|883, 732, 178|
|119|254|163|1175, 894|
|78|141|11|1329, 1697, 856, 144|

For an **arbitrary nonempty tournament** \(M\) substituted in the same place, these marked reversals change the height by exactly \(2P_1(M)=2H(M)\), since the exterior coefficient discrepancy is confined to the first nonempty slot. Thus their unit-gap substitutions are precisely the transitive tournaments. To see necessity, a nontransitive tournament has a nontrivial strong component. Its directed Hamilton cycle gives at least three distinct paths by rotation, each of which extends by concatenating paths on the other ordered components. Sufficiency follows from the unique transitive order. Consequently the perfect-number exclusion below also exhausts every nonempty replacement tournament on these seven marked arcs.

Prime 19 has nonzero derivative on every phase of the first six, and 67 works on every phase of the seventh. Thus two primes select the width for these persistent physical reversals. At target 3744 the selected masks are 11, 5, width four and prime 67, giving a single flip from 3743 to 3745.

For every odd \(p\ge3\), the Euclid–Euler number \(N_p=2^{p-1}(2^p-1)\) equals 28 modulo 36. Indeed it is zero modulo four and one modulo nine. At any width at least two the seven midpoint functions are periodic with period dividing six modulo 36: reduce the powers of 2, 3 and 4 modulo four and nine. The width-one values miss. Exactly two rows survive, rows 2 and 7, each at \(s\equiv5\pmod6\). Row 2 gives the midpoint \(-64+238\cdot2^s\), equal to 0 or 2 modulo five at odd widths, whereas \(N_p=2u^2-u\) with \(u=2^{p-1}\equiv\pm1\pmod5\) is 1 or 3. Row 7 gives
\[
       N_{s+4}<144\cdot4^s+280\cdot3^s-7\cdot2^s-72<N_{s+5}.
\]
The first difference is \(16\cdot4^s+280\cdot3^s+2^s-72>0\); the second is \(368\cdot4^s-280\cdot3^s-9\cdot2^s+72>0\) at all \(s\ge1\). Thus these seven families miss every perfect-number target with odd exponent, including 61. This argument uses the ordinary size of the surviving expression after modular selection.

From this directory the replay uses Python's standard library and a C++17 compiler.

~~~sh
c++ -std=c++17 -O3 heegner_slots.cpp -o /tmp/heegner-slots
/tmp/heegner-slots > /tmp/heegner-slots.txt
python3 -B heegner_decode.py --bank /tmp/heegner-slots.txt --census
python3 -B heegner_decode.py --midpoint 3744
python3 -B heegner_decode.py --bank /tmp/heegner-slots.txt --height 437
~~~
