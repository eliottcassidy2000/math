# A Lorentz boost on divisor splits, and the split that can still lift

**Proved necessary restrictions, not a prime coverage theorem.** Use the physical two-arm identity \(H=de+2uv\) for \(C_3[A,B,1]\), where the arms have ordinary reads \(u,v\), two-slot reads \(R,S\), and gaps \(d=R-u,e=S-v\). The right arm is **7-carried**, meaning \(7\mid v,S\). Keep the scope of the earlier [shifted-divisor and primitive-gap chart](../blue_gate_spectrum_20261117/DIVISORS.md): demand \(d=s^2,e=t^2\), with odd positive coprime \(s,t\), and \(H=7p\) for a prime \(p>3\). In the principal theorem both arms are nonsingletons.

## The arithmetic intertwiner and where it loses physical meaning

For fixed odd leg \(W=st\), primitive choices give
\[
 C=(s^2+t^2)/2,\quad E=(t^2-s^2)/2,\qquad C^2-E^2=W^2.
\]
Transferring a full prime power \(P\mid t\) from \(t\) to \(s\) changes the null coordinates \((C-E,C+E)\) by \(\operatorname{diag}(P^2,P^{-2})\). Thus it is an exact integral-on-its-domain Lorentz boost, preserving the odd leg and conjugating divisor complement to \(E\mapsto-E\). Only full prime powers can move between primitive roots. This yields the full labelled \(2^{\omega(W)}\)-cube; the exponent-one layer for a squared prime in the divisor diagram is not in the primitive picture.

For \(F=\#\) proper nontrivial divisors, \(U=\#\) proper prime divisors, and \(S_0=\#\) proper nontrivial squarefree divisors, the nontrivial solutions of \(F=U+S_0\) are \(W=\ell^3\) or \(W=\ell^2 q r\) with distinct primes in the second case. Indeed a squarefree composite fails because \(F=S_0<U+S_0\). For nonsquarefree \(W=\prod_{i=1}^a p_i^{e_i}\) the identity is \(\prod(e_i+1)=2^a+a+1\). For \(a\ge4\) the lower bound \(3\cdot2^{a-1}\) is larger; the cases \(a=1,2,3\) give a cube, no solution, and exponent pattern \((2,1,1)\). The prime case is the trivial zero equality.

For any tournament with \(n\) vertices, the odd-cycle packing formula and complementary-subset convolution give
\[
 Q_2(T)=\sum_{\mathcal C}2^{\,n-\sum_{D\in\mathcal C}|D|+2|\mathcal C|}.
\]
To see this, color each packed odd cycle wholly in one of the two slots, independently color each uncovered vertex, and include the factor \(2^{|\mathcal C|}\) from the Hamilton reads of the two slots. Thus \(Q_2\equiv(-1)^n H\bmod3\), \(4\mid Q_2\) for \(n\ge2\), and \(8\mid Q_2\) for even \(n\ge4\) (a single odd cycle leaves at least one uncovered vertex; two cycles already contribute a factor 16).

Now \(3\nmid uv\), or \(3\mid H\). For a prime quotient the two arms cannot both be odd, since then \(d\equiv u,e\equiv v\bmod3\), giving \(3\mid H\). They cannot both be even since that would put 3 on both roots. Hence the unique even arm holds the entire 3-power, and the right arm holds the entire 7-power. Moving either of these powers is a genuine **nonlift** on the fixed physical parity and carrier data. Even before other constraints the boost cube falls to at most \(2^{\omega(W)-2}\) assignments. In particular the prime cubes and trivial prime case in \(F=U+S_0\) cannot occur.

## Only one squared-prime placement survives with two nonsingletons

Write \(W=3^\alpha7^\beta\ell^\gamma\), where \(\ell\ne3,7\) is an odd prime and \((\alpha,\beta,\gamma)\) is a permutation of \((2,1,1)\). If the right arm is even, it carries 3 and7. It cannot also get \(\ell\), since then \(d=1\) and \(d\ge |A|u\) would force the left arm to be a singleton. The only candidate right gaps are accordingly \(441,3969,21609\). If the right arm is odd, it cannot receive \(\ell\) either: the resulting left even gap would be 9 or81. In even order at least four a square gap forces \(u\equiv7\bmod8\), and \(d\ge|A|u\); at gap 81 this leaves \(u=7\) (forbidden) or \(u=15\) only at order four (impossible), and at gap9 none. Order two has gap3. Thus \(\ell\) lies left in the right-odd case; the right gap then is 49 or2401.

Since \(v\) is a nonzero multiple of 7 other than 7 or21, it is at least 35. The pair-swap bound \(H(T)\le n!/2^{\lfloor n/2\rfloor}\) excludes orders at most five for \(v\ge35\) (the bound at five is30), and \(e\ge |B|v\) excludes gap49. The bound follows by independently swapping each disjoint consecutive pair in a permutation, with at most one orientation per pair compatible with a path. If the right-even gap is441, then \(v\equiv7\bmod8\), \(7\mid v\), and \(35\le v\le441/6\); only \(v=63\) meets these displayed tests, and it is excluded by \(3\mid v\). This excludes the exponent-square placement \(\gamma=2\).

Here are complete finite exclusions for gaps 2401 and 3969. Put \(v=7w\). Square gaps give \(w\equiv1\bmod4\); even right order improves this to \(w\equiv1\bmod8\). Also \(3\nmid w\). For 2401 the right order is odd and at least seven. Hence \(e\equiv v\equiv w\pmod3\), and the nonzero square \(e=2401\) forces \(w\equiv1\pmod3\); together with \(w\equiv1\pmod4\) and the size bound, this gives \(w\in\{13,25,37,49\}\) and
\[
 (v,Q_2)\in\{(91,2492),(175,2576),(259,2660),(343,2744)\}.
\]
All have just one nontrivial strong component, if present: the only extra possible factorization is \(175=5\cdot35\), which would force \(Q_2/v\ge5\cdot6>2576/175\) by multiplication of the split-cover ratios of the components. (Counts 7 and21 are unavailable.) The strong component has order at least seven by the pair-swap bound \(H\le90\) at six. Singleton padding is impossible: for \(v=91,259\) the exact valuation \(v_2(Q_2)=2\) allows none since a nontrivial component has \(4\mid Q_2\); for \(v=175,343\) padding would make \(Q_2/v\ge16\). A strong order-nine tournament has \(Q_2\ge2676\): using nested strong deletions and vertex-pancyclicity gives at least \(n-j+1\) cycles of each odd length \(j\), contributing \(512+7\cdot256+5\cdot64+3\cdot16+4\). Strong order at least eleven already exceeds the four reads by the empty and triangle packings. And \(v=343\) has \(Q_2/v=8\), requiring order at most seven. Therefore only strong order seven need be checked for gap2401.

For gap3969, \(w\ge17\), so \(v\ge119>90\), and the even order is at least eight. This leaves \(w\in\{17,25,41,49,65\}\). Again there is only one nontrivial component: \(175=5\cdot35\) would force ratio at least30, and any split of455 into nontrivial factors forces \(Q_2/v\ge16\); the other candidates cannot split without an unavailable count 7. For \(v=287,343,455\), the ratio is less than16, so no padding is possible. For \(v=119,175\) at most one singleton can be padded, since two force ratio at least32. A strong even right arm of order ten or larger already has \(Q_2\ge6032\) by the same odd-cycle lower bounds. If there is one padding vertex and the strong component has odd order at least nine, then \(Q_2\ge2\cdot2676\), also too large. Thus the only checks for gap3969 are strong order eight and order-seven strong cores with one singleton.

The supplied [census replay](fus-census.cpp) fixes a spanning directed cycle in the strong core, enumerating the other \(2^{14}\) or \(2^{20}\) arc orientations. It finds no order-seven profile for either condition. At order eight the only profile carrying 7 with \(Q_2-H=3969\) is \((231,4200)\), already rejected by \(3\mid H\). Fixing the cycle loses no strong tournament, by Camion's theorem (also implicit in [Moon’s vertex-pancyclicity theorem](https://doi.org/10.4153/CMB-1966-038-7)).

Consequently the sole remaining necessary locus for two nonsingleton arms is
\[
 W=3\cdot7^2\ell,\quad (s,t)=(\ell,147),\quad
 d=\ell^2,\ e=21609,\quad |A|\text{ odd},\ |B|\text{ even}.
\]
This does not assert any physical profile with gap21609. Indeed the right order is at most 14 because \(2^{|B|}\le Q_2\le(1+1/8)21609<2^{15}\), with even order at least eight as above. The target-first test on this remaining bank is \(p=3087\ell^2+2uw\), with actual arm reads and fixed right gap. A singleton left arm is separate; it puts every prime power on the right, and gives \(p\equiv1\bmod8\) whereas nonsingletons give \(p\equiv5\bmod8\). No broader primitive Pythagorean or sevenfold-prime lane is excluded.


Replay with a standard C++17 compiler:

    c++ -O3 -std=c++17 fus-census.cpp -o /tmp/fus-census
    /tmp/fus-census 7
    /tmp/fus-census 8

The [recorded transcript](fus-census.out) gives both exact totals. At order eight, 1,232 of the Hamilton-framed masks yield the sole excluded profile, not 1,232 isomorphism types.

The replay is `c++ -O3 -std=c++17 fus-census.cpp -o fus-census` followed by `./fus-census 7` and `./fus-census 8`.
