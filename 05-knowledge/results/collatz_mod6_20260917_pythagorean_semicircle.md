# Pythagorean semicircle lane: odd-square identity, Thales chart, Gaussian doubling, and the 29 coincidence

**Status: PROVED elementary identities (Theorems 1-9, most classical, no novelty claim except the all-PPT doubling bijection of Theorem 8) + FINITE-EXACT censuses + CITED closure of reading (b1) (Fermat's X^4-Y^4=Z^2) + three REFUTED readings (literal (b1) on 3-4-5; 29 as a map; the plus sign in x^2+{0,1,2}) + SCOPE items (no braids2 inheritance; user-wording matches unverifiable by auditors). Nothing remains OPEN inside the lane; the Ljunggren-type boundary of the first draft is retired.** No Collatz, Goldbach, or LRC claim. Session collatz-mod6-20260917 (machine mac-mini), lane `pythagorean_semicircle` (wildcard: geometry); finalized 2026-09-21 after two independent audits (recompute and proof-audit), whose corrections are all applied below.

## Inheritance and concept board

Inherited and not re-derived: [collatz](arithmetic_braids_20260917_collatz.md) (three-row typing; "the order of 2 modulo 9 is 6" row-exponent law with rows 9j+2, 9j+5, 9j+8), [summand](arithmetic_braids_20260917_summand.md) (doubling forest; section 2 squaring forest x->x^2 with the hostile "square-root ancestry is not Pythagorean primitivity"), [geometry](arithmetic_braids_20260917_geometry.md) (section 5: Gaussian squaring as fresh leg primes and content, via [THM-3336](../../01-canon/theorems/THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md)), [divisors](arithmetic_braids_20260917_divisors.md). From canon: [THM-3333](../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md) (3), the Gaussian-square lift Phi(m,n)=(m^2-n^2,2mn,m^2+n^2); [THM-3341](../../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md) (6)-(12) (negative Pell selector, recursions m_(k+2)=6m_(k+1)-m_k and t_(k+2)=6t_(k+1)-t_k+2, odd-index Pell hypotenuses), (23) Gamma(a,b,c)=(|a^2-b^2|,2ab,c^2) = Gaussian squaring, (24) the Berggren A_+ ray (3,4,5),(21,20,29),(119,120,169),(697,696,985),..., (26) (3,4,5)->(7,24,25), (21,20,29)->(41,840,841), (31) the count 2^(omega(m)-1); [THM-3334](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md) (36), Tripathi's fixed-hypotenuse count |X_c|=2^(omega(c)-1) for c>=3 with no prime 3 mod 4; [THM-3357](../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md) (22), the same ray; [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md) (17) the 3-4-5 specialization (t+1,t)=(2,1), (18) 29=5^2+2^2, (20) B^3=-I, B^6=I; [THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md) (29)-(33a): G_D(y)=(y^2-D)/b cycles -(a+b)->h->-(b-a) iff 2a(4a-3b)=0, i.e. (a,b,h)=(3k,4k,5k), D=29k^2, and (33a) already types 29=5^2+2^2=2^2+3^2+4^2 as a non-analogy; [THM-3335](../../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md) Thm 2.1, the hypotenuse = even square + 1 family; [THM-4057](../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md) header, Pell/Fibonacci cycle families in the depth gauge; [THM-2142](../../01-canon/theorems/THM-2142-the-half-angle-bridge-ab-monoid-is-the-ctu-cyclotomic-skeleton.md) section 1, b(a(cos theta))=(1+cos theta)/2=cos^2(theta/2), with the exact map stated in Theorem 1 below (the first draft cited it without a map).

Concept board. Closest proved mechanism: Gaussian squaring Gamma of THM-3341 (23), which is the angle doubling of this lane (Theorem 7, inherited). Canonical hostile: the summand note's "square-root ancestry is not Pythagorean primitivity" for arbitrary integers; on hypotenuses it is answered by the bijection of Theorem 8. Corrected near miss: the first draft's "t/s runs through the convergents of sqrt2-1" (REFUTED by both audits; t/s are the intermediate fractions, n/m the convergents, and the smaller-angle half-tangent alternates between them). Least-used sidecar: tan(theta/2), the coordinate the continuous (e/d,l) chart loses and which alone decides PPT-ness and the leg parity (Theorem 5).

SCOPE: the arithmetic_braids2_20260917_* notes contain no Pythagorean, Pell, or Gaussian-squaring content, so no braids2 inheritance applies to this lane (no map found). SCOPE: every statement of the form "matches the user's wording" (readings (b2), the parity clause, the x^2+{0,1,2} graphs) is HEURISTIC prose interpretation; neither auditor had the user's text, and the arithmetic is what is claimed.

Throughout, a PPT is (a,b,c) with a=m^2-n^2 odd, b=2mn even, c=m^2+n^2, m>n>0 coprime of opposite parity; phi is the angle opposite the even leg b and psi the angle opposite the odd leg a; in Theorems 3-5 theta is the smaller angle. Universe of Theorems 1, 3, 5: all 158 PPTs with c<=1000.

## 1. Hypotenuse plus the even leg is an odd square; (s,t) is the identity

**Theorem 1 (PROVED; classical, no novelty claim).** c+b=(m+n)^2 and c-b=(m-n)^2 are odd squares; c+a=2m^2 and c-a=2n^2 are never squares. With s=m+n, t=m-n, s>t>=1 are odd and coprime and

    a=st,   b=(s^2-t^2)/2,   c=(s^2+t^2)/2,   s^2=c+b.

*Proof.* Direct expansion. m,n of opposite parity make s,t odd. gcd(s,t) divides s+t=2m and s-t=2n, hence divides 2gcd(m,n)=2, and is odd, so equals 1. 2m^2 is not a square because sqrt2 is irrational. Conversely any odd coprime s>t gives coprime opposite-parity (m,n)=((s+t)/2,(s-t)/2). QED. Checked on the universe.

So the user's claim (a) holds for **exactly one** leg, the even one. "Half of their identity" becomes exact:

**Theorem 2 (PROVED).** For odd s>=3 the number of PPTs with c+b=s^2 is phi(s)/2.

*Proof.* The fibre is {t: 1<=t<s, t odd, gcd(t,s)=1}. Since s is odd, t->s-t is an involution on the units mod s that swaps parity, so exactly half of the phi(s) units are odd. QED. Checked s=3..201 (100 values). First fibres:

| s | c+b=s^2 | fibre |
|---|---|---|
| 3 | 9 | (3,4,5) |
| 5 | 25 | (5,12,13), (15,8,17) |
| 7 | 49 | (7,24,25), (21,20,29), (35,12,37) |
| 9 | 81 | (9,40,41), (45,28,53), (63,16,65) |
| 15 | 225 | (15,112,113), (105,88,137), (165,52,173), (195,28,197) |

Boundaries: s=1 has empty fibre; t=s would give the degenerate (s^2,0,s^2).

**Half-angle reading and the THM-2142 map (PROVED).** tan(phi/2)=n/m and tan(psi/2)=t/s, by tan(2u)=2u/(1-u^2); verified as exact Fractions on the universe. The (n,m,sqrt c) triangle is the half-angle triangle of phi **only** (its other acute angle is pi/2-phi/2, and m/n != t/s on every PPT); psi is halved by the (t,s,sqrt(2c)) triangle. The map to THM-2142 section 1 is exact: on a PPT,

    cos^2(phi/2)=(c+a)/(2c)=m^2/c,   sin^2(phi/2)=(c-a)/(2c)=n^2/c,
    cos^2(psi/2)=(c+b)/(2c)=s^2/(2c), sin^2(psi/2)=(c-b)/(2c)=t^2/(2c),

so "hypotenuse plus or minus a leg" is 2c times the half-angle cosine or sine squared, and Theorem 1 is THM-2142's functional b(a(cos theta))=cos^2(theta/2) evaluated at PPT angles (checked on the 158 PPTs, together with cos^2=1/(1+tan^2) for both half-angles).

## 2. Readings of "hypotenuse one more than a square, altitude its square root"

**(b1) literal, k a positive integer: hypotenuse C=k^2+1, altitude to the hypotenuse H=sqrt k. REFUTED as a description of 3-4-5; PROVED (modulo a classical citation) to have no rational-leg member at all.**

Legs a,b satisfy a^2+b^2=C^2 and ab=CH, hence (a+b)^2=C^2+2CH and (a-b)^2=C^2-2CH; a real right triangle exists iff C>=2H. With u=sqrt k,

    C-2H = u^4-2u+1 = (u-1)(u^3+u^2+u-1),

so a triangle exists iff k>=1 or k<=k0=u0^2=0.295597742522, u0=0.543689012692 the real root of the cubic. At k=1 the legs are sqrt2, sqrt2 (hypotenuse 2). The 3-4-5 has hypotenuse 5=2^2+1 but altitude ab/c=12/5, not sqrt2: **minimal witness**, the triple itself. (For rational k the argument below needs k a rational square; k is taken a positive integer.)

Rational legs: ab/C=sqrt k rational forces k=j^2, C=j^4+1, ab=j(j^4+1), (a+-b)^2=(j^4+1)(j^4+1+-2j). Direct search j<=2000 finds no rational-leg member. For j even, gcd(j^4+1,2j)=1, so j^4+1 must be a square, impossible since j^4<j^4+1<(j^2+1)^2 (PROVED). For j odd the gcd is 2 and (j^4+1)/2 must be a square: j^4+1=2w^2, i.e. ((j^2-1)/2,(j^2+1)/2,w) is a near-isosceles PPT whose leg sum j^2 is a square. Among the numerators x of x^2-2y^2=-1 through Pell index 400 only x=1 is a square (FINITE-EXACT), and j=1 is not a member either: (j^4+1+2j)/2=2 is not a square, its legs are sqrt2, sqrt2.

**Fermat closure (CITED).** If j^4+1=2w^2 then w^4-j^4=((j^4-1)/2)^2, by the identity (j^4+1)^2-4j^4=(j^4-1)^2 (checked symbolically). For j>1 this is a solution of X^4-Y^4=Z^2 with XYZ != 0, which Fermat's descent theorem excludes (classical; e.g. Hardy-Wright ch. XIII, Mordell, *Diophantine Equations*, ch. 4; the theorem numbers were not verified here, hence CITED and not PROVED-in-file). So j=1, and reading (b1) has **no rational-leg member for any positive integer k**. Corroboration: odd j<=100000 with (j^4+1)/2 a square: only j=1. The first draft's "UNCITED-RECOLLECTION (Ljunggren-type)" boundary is retired; Ljunggren's equation x^2-2y^4=-1 is a different equation, and the draft's open question on x^4+1=2w^2 is answered (only x=w=1).

**(b2) Euclid column n=1: (k^2-1, 2k, k^2+1), k>=2 (k=1 is the degenerate (0,2,2)). HEURISTIC reading match; PROVED family facts.** Hypotenuse is exactly k^2+1; primitive iff k even, since gcd(k^2-1,2k)=gcd(k^2-1,2); altitude 2k(k^2-1)/(k^2+1); the angle opposite the even leg is 2arctan(1/k)->0, so the normalized vertex converges to the diameter endpoint.

| k | triple | primitive | altitude |
|---|---|---|---|
| 2 | (3,4,5) | yes | 12/5 |
| 3 | (8,6,10) | gcd 2 | 24/5 |
| 4 | (15,8,17) | yes | 120/17 |
| 5 | (24,10,26) | gcd 2 | 120/13 |
| 6 | (35,12,37) | yes | 420/37 |
| 8 | (63,16,65) | yes | 1008/65 |
| 10 | (99,20,101) | yes | 1980/101 |

The only length equal to "the square root of that same number" is the leg k of the half-angle triangle (k,1,sqrt(k^2+1)). Which reading the user intended is HEURISTIC and was unverifiable by the auditors.

**(b6) THM-3335 family (CITED, not re-proved).** Hypotenuse = even square + 1 with the even leg equal to that square: (3,4,5),(17,144,145),(99,4900,4901),(577,166464,166465) (rows k=1..4 of THM-3335's table). It also starts at 3-4-5, but nothing in it equals sqrt(4)=2.

## 3. The Thales chart (e/d, l, theta) and its rational lattice

**Correction (PROVED, Thales).** The right-angle vertex lies on the circle whose *diameter* is the hypotenuse (radius 1/2 after normalizing c=1), not on a circle of radius 1; by the inscribed-angle theorem its central angle is 2theta.

**Theorem 3 (PROVED).** With theta the smaller angle, theta=arcsin(a1/c) in (0,pi/4) for a1 the shorter leg: e=a1^2/c^2=sin^2theta, d=b1^2/c^2=cos^2theta, l=a1b1/c^2=sin(2theta)/2, so

    e+d=1,   l^2=ed,   e/d=(a1/b1)^2=tan^2theta,
    vertex = ((d-e)/2, l) = (b1^2-a1^2, 2a1b1)/(2c^2),

i.e. the semicircle vertex is one half of the **angle-doubled** triple (geometric-mean theorem a1^2=ce, b1^2=cd, a1b1=cl). Verified exactly on the universe, including the vertex on x^2+y^2=1/4.

**Theorem 4 (PROVED).** (e/d, l, theta)=(tan^2theta, sin(2theta)/2, theta) is strictly increasing on (0,pi/4]: d(tan^2)/dtheta=2tan/cos^2>0, dl/dtheta=cos2theta>0. The endpoint (1,1/2,pi/4) is the isosceles triangle, never a PPT. The user's "vary {e/d,l,theta} from epsilon to {1,1/2,pi/4}" is one parameter, theta.

**Theorem 5 (PROVED).** e/d and l depend on tan theta only, so are rational for every rational-slope right triangle; witness legs (1,2): e/d=1/4, l=2/5, hypotenuse sqrt5. PPT angles are exactly the theta in (0,pi/4) with tan(theta/2) rational (a rational point on the circle iff a rational half-angle tangent), dense in (0,pi/4). If tan(theta/2)=p/q reduced: p,q of opposite parity iff theta is opposite the even leg ((p,q)=(n,m)); both odd iff theta is opposite the odd leg ((q,p)=(s,t)). FINITE-EXACT: over the 158 PPTs with c<=1000 the smaller angle is opposite the odd leg in 80 cases and opposite the even leg in 78, so both charts occur (a statement that only opposite-parity fractions arise would be false). Lost coordinate of the (e/d,l) chart: rationality of c; sidecar: tan(theta/2).

**Extremal families.** theta->0: (k^2-1,2k,k^2+1) (k even) and (s,(s^2-1)/2,(s^2+1)/2) (s odd). theta->pi/4: near-isosceles |a-b|=1.

**Theorem 6 (PROVED; INHERITED family, THM-3341 (6)-(12), (24); THM-3357 (22)).** |a-b|=1 iff (m-n)^2-2n^2=+-1; all solutions form the Pell orbit (m,n)->(2m+n,m) from (2,1), with c_{i+1}=6c_i-c_{i-1}, x_{i+1}=6x_i-x_{i-1}+2 for the short leg, (m+n)^2-2m^2=+-1, and hypotenuses the odd-index Pell numbers pell_3, pell_5, ... = 5, 29, 169, 985, 5741, 33461, 195025, 1136689 (THM-3341 (12) also lists the degenerate m_0=pell_1=1).

| (a,b,c) | (m,n) | (s,t) | tan(phi/2)=n/m | tan(psi/2)=t/s | smaller angle |
|---|---|---|---|---|---|
| (3,4,5) | (2,1) | (3,1) | 1/2 | 1/3 | psi |
| (21,20,29) | (5,2) | (7,3) | 2/5 | 3/7 | phi |
| (119,120,169) | (12,5) | (17,7) | 5/12 | 7/17 | psi |
| (697,696,985) | (29,12) | (41,17) | 12/29 | 17/41 | phi |
| (4059,4060,5741) | (70,29) | (99,41) | 29/70 | 41/99 | psi |
| (23661,23660,33461) | (169,70) | (239,99) | 70/169 | 99/239 | phi |
| (137903,137904,195025) | (408,169) | (577,239) | 169/408 | 239/577 | psi |
| (803761,803760,1136689) | (985,408) | (1393,577) | 408/985 | 577/1393 | phi |

**Corrected near miss.** tan(phi/2)=n/m = 1/2, 2/5, 5/12, 12/29, ... are the **convergents** of sqrt2-1=[0;2,2,...]; tan(psi/2)=t/s = 1/3, 3/7, 7/17, 17/41, ... are the **intermediate fractions** (mediants of consecutive convergents; ratios of consecutive companion Pell numbers 1,3,7,17,41,99,...) and are never convergents (checked exactly for all eight rows). The first draft's "t/s runs through the convergents" is REFUTED (minimal witness: 3/7 is not a convergent, 2/5 is). With theta the smaller angle, tan(theta/2) alternates between the two sequences: 1/3, 2/5, 7/17, 12/29, 41/99, 70/169, 239/577, 408/985, monotone and below sqrt2-1. FINITE-EXACT: the near-isosceles hypotenuses with c<=10^6 are exactly {5, 29, 169, 985, 5741, 33461, 195025}.

**The 29 test (REFUTED as a map; PROVED locus).** THM-4146 has D=a^2+2b^2-ab=29 at 3-4-5 and 29=5^2+2^2; the Pell pair (5,2) gives c=29 for (20,21,29). Exact iteration of G_D from -(a+b)=-41 fails to close under **either** labeling: (a,b)=(20,21): D=862, orbit -41 -> 39 -> 659/21 -> 54139/9261; (a,b)=(21,20): D=821, orbit -41 -> 43 -> 257/5 -> 11381/125; neither passes through 29, and THM-4146 (31) needs 3a-b=h, which fails for both (60-21 and 63-20 are not 29). As polynomials in Euclid (m,n),

    D = m^4-2m^3n+6m^2n^2+2mn^3+n^4 = c^2+b(b-a),   c_next = (2m+n)^2+m^2 = 5m^2+4mn+n^2,
    at n=1:  D-c_next = m(m-2)(m^2+1),
    n>=2:    D-c_next = m^2((m-n)^2-1) + m^2(5n^2-4) + 2mn(n^2-2) + n^2(n^2-1) > 0

(every term >=0 for m>n>=2 and the last is >0), so D(m,n)=c_next(m,n) only at (m,n)=(2,1) for **all** m>n>=1 (the first draft's search 300>=m>n>=1, which found only (2,1), is now corroboration). Also b(b-a)=4=m^2 only at 3-4-5. THM-4146 (33a) already types 29=5^2+2^2=2^2+3^2+4^2 as a non-analogy; the new content here is only the negative (20,21,29) test.

## 4. Angle doubling is Gaussian squaring; the PPT doubling forest

**Theorem 7 (PROVED; INHERITED = THM-3341 (23), no novelty claim).** (a,b,c)->(2ab,|b^2-a^2|,c^2) is (b+ai)^2 with norm c^2, i.e. THM-3341's Gamma; in y=2cos theta it is y->y^2-2 (2(2x^2-1)=(2x)^2-2). THM-3333's Phi(m,n)=(m+ni)^2 is the same map, so half-angle triangle -> PPT is one doubling step. THM-3341 (26) already records (3,4,5)->(7,24,25) and (21,20,29)->(41,840,841).

**Theorem 8 (PROVED; the all-PPT bijection is the only new statement).** A PPT is the double of a PPT iff its hypotenuse is a perfect square iff its Euclid parameters (m,n) are the legs of a PPT; the half is unique; for every hypotenuse c, doubling is a bijection PPT(c)->PPT(c^2), both of size 2^(omega(c)-1) (for c a PPT hypotenuse, i.e. every prime factor 1 mod 4: THM-3334 (36), Tripathi; THM-3341 (31) uses the same count). THM-3341 section 2 proves "square hypotenuse iff Gamma-image" on the U-spine only.

*Proof.* If (a,b,c) is a PPT then gcd(2ab,b^2-a^2)=1 (b^2-a^2 odd, gcd(a,b)=1), so the double is a PPT with hypotenuse c^2. Conversely if (a',b',c^2) is a PPT with Euclid (m,n), then m^2+n^2=c^2 with gcd(m,n)=1, so (n,m,c) is a PPT and its double (2mn,m^2-n^2,c^2)=(b',a',c^2); uniqueness of Euclid parameters gives uniqueness of the half. Injectivity of doubling on PPT(c) follows from the same uniqueness; surjectivity onto PPT(c^2) is the converse. QED

FINITE-EXACT census, c<=10^4 (1593 PPTs): exactly 16 have a square hypotenuse, each with a unique half, and #PPT(c<=100)=16:

| double | half | (m,n) | sqrt c |
|---|---|---|---|
| (7,24,25) | (3,4,5) | (4,3) | 5 |
| (119,120,169) | (5,12,13) | (12,5) | 13 |
| (161,240,289) | (15,8,17) | (15,8) | 17 |
| (41,840,841) | (21,20,29) | (21,20) | 29 |
| (527,336,625) | (7,24,25) | (24,7) | 25 |
| (1081,840,1369) | (35,12,37) | (35,12) | 37 |
| (1519,720,1681) | (9,40,41) | (40,9) | 41 |
| (1241,2520,2809) | (45,28,53) | (45,28) | 53 |
| (721,5280,5329) | (55,48,73) | (55,48) | 73 |
| (2047,3696,4225), (3713,2016,4225) | (33,56,65), (63,16,65) | (56,33), (63,16) | 65, omega=2 |
| (3479,1320,3721) | (11,60,61) | (60,11) | 61 |
| (959,9360,9409) | (65,72,97) | (72,65) | 97 |
| (4633,5544,7225), (6887,2184,7225) | (77,36,85), (13,84,85) | (77,36), (84,13) | 85, omega=2 |
| (4879,6240,7921) | (39,80,89) | (80,39) | 89 |

Chain: (3,4,5)->(24,7,25)->(336,527,625)->(354144,164833,390625). Hence the hypotenuse ray c,c^2,c^4,... is a ray of the summand note's squaring forest (section 2), and on hypotenuses square-root ancestry lifts **bijectively** to PPT halving; the summand hostile concerns non-hypotenuse integers, and the geometry note's Gaussian squaring (fresh leg primes, content) is a different question from this bijection. Roots of the PPT forest: PPTs with non-square hypotenuse.

**Theorem 9 (PROVED).** PrePer(y^2-2,Q)={-2,-1,0,1,2}: |y|>=3 escapes since |y^2-2|>|y|, and a non-integer rational never repeats because every denominator prime has its valuation doubled. Graph 2->2, -2->2, 0->-2, 1->-1, -1->-1 = angle doubling on {0, pi, pi/2, pi/3, 2pi/3}. y^2-2 has exact period-3 orbits at 2cos(2pi k/9), gcd(k,9)=1 (minimal polynomial y^3-3y+1; k=3 gives the fixed point -1) and at 2cos(2pi k/7), gcd(k,7)=1 (y^3+y^2-2y-1). For N=9 the angle orbit {1,2,4,8,7,5} has order 6=ord_9(2), the same fact as the inherited row-exponent law, and folds by -1 to the 3-cycle {1,2,4}; the Collatz row residues {2,5,8} meet each pair {k,-k} exactly once (since -1 swaps the classes 1,2 mod 3). N=7 has ord_7(2)=3 already, so its 3-cycle needs no -1 fold: the "order six over three via a central -1" shape is N=9-specific. It is the shape of THM-4139 (20), B^3=-I, but a shared shape, not a map (cubic irrationals versus quarter-integers; x^2-29/16 and x^2-2 are not conjugate).

## 5. Sign audit of the user's x^2+{0,1,2} (REFUTED as written)

Exact integer graphs inside [-3,3]: x^2 has fixed points 0,1 with -1->1; x^2-1 has the cycle -1<->0 with 1->0; x^2-2 has fixed points -1,2 with 1->-1 and 0->-2->2. These are the described graphs. x^2+1 and x^2+2 have no integer periodic points (x^2+c>x for c>=1). Survivor: the minus-sign family, consistent with the user's x^2-7/4 and x^2-29/16. Typing: x^2 is z->z^2 on the unit circle (angle doubling on the circle itself), x^2-2 its trace projection y=z+1/z (Theorem 7); x^2-1 has no angle reading (not conjugate to a power or Chebyshev map). The match to the user's prose is HEURISTIC (SCOPE above).

## 6. Typed analogies

| source | target | map | preserved | lost | sidecar | decisive test | verdict |
|---|---|---|---|---|---|---|---|
| (n,m,sqrt c) half-angle triangle of phi | PPT | (m+ni)^2 = THM-3341 Gamma | phi doubles, norm squares | psi (halved by (t,s,sqrt(2c)) instead) | parity of tan(theta/2) | Theorem 8 census | PROVED bridge, inherited map |
| THM-2142 b(a(cos theta))=cos^2(theta/2) | Theorem 1 | cos^2(phi/2)=m^2/c, cos^2(psi/2)=s^2/(2c) | half-angle functional | nothing (exact) | which leg | 158-PPT check | PROVED identity |
| hypotenuse squaring ray (summand section 2) | PPT doubling ray | Gaussian norm | c->c^2 | argument theta | the Gaussian integer | bijection PPT(c)->PPT(c^2) | PROVED bijective lift (new) |
| Thales semicircle vertex | doubled PPT | central angle 2theta | rational point | none | none | vertex=(b^2-a^2,2ab)/2c^2 | PROVED identity |
| user's continuous {e/d,l,theta} | PPT lattice | theta | monotone order | rationality of c | tan(theta/2) | legs (1,2) | PROVED, lattice dense |
| 29 in THM-4146 | 29 = hypotenuse of (20,21,29) | none | numeral only | everything | none | G_D orbits, D=862 and D=821; locus {(2,1)} | REFUTED |
| ord_9(2)=6 folded by -1 | THM-4139 (20) B^3=-I | none | order 6 -> 3 shape | fields differ; N=7 lacks the fold | none | no common object | shape only |
| (b1) hypotenuse k^2+1, altitude sqrt k | 3-4-5 | none | none | none | none | altitude 12/5; Fermat closure | REFUTED, family empty |
| braids2 notes | this lane | none | none | none | none | grep | SCOPE, no map found |

## 7. Reproduction

    cd /tmp/math-wt-collatz-mod6-b
    python3 04-computation/experiments/collatz_mod6_20260917_pythagorean_semicircle.py > 05-knowledge/results/collatz_mod6_20260917_pythagorean_semicircle.out
    python3 -O 04-computation/experiments/collatz_mod6_20260917_pythagorean_semicircle.py | cmp - 05-knowledge/results/collatz_mod6_20260917_pythagorean_semicircle.out

Final script sha256 020b0c1f9356e5d4ec563541dad8fb2d7879f6eeb98e466047f2e402537f2891; output sha256 0c0ea35fd5569759fa7310851b4ca06ad8ed86459fca2a0a44a53e5d2704d40b (217 lines); normal and `-O` runs byte-identical (the `require` helper raises, so `-O` disables nothing); run metadata about 20 s, under 100 MB resident. Universes: PPTs c<=1000 (158), c<=10^4 (1593), c<=10^6 (near-isosceles census), s<=201, j<=2000 and odd j<=100000, Pell index<=400, (m,n) agreement 300>=m>n>=1 as corroboration of the all-(m,n) proof. Section S6 of the output carries the audit-driven checks (THM-2142 map, parity split 80/78, half-angle-triangle precision). Independent audit scripts of the same session: `04-computation/experiments/collatz_mod6_20260917_pythagorean_semicircle_audit_recompute.py` (Euclid-free brute-force PPT enumeration) and `..._audit_proof-audit.py` ((s,t)-chart enumeration, own polynomial code).

Provenance: the lane script was recovered from agent transcripts on 2026-09-21 after the worktree was pruned; the recovered draft passed unchanged, and the 2026-09-21 finalization added the corrections above (one false sentence in the draft output, the convergents claim, was replaced by the checked true statement).

## 8. Stopping boundary / next question

The lane is closed as geometry: every reading of the user's claims is either PROVED, CITED-closed, or REFUTED with a minimal witness, and the only residual novelty is the all-PPT bijection PPT(c)->PPT(c^2) of Theorem 8 (THM-3341 has the U-spine case). Retired: the quartic-Pell open question (Fermat). Not attempted, and not claimed: whether the additive doubling forest (odd core u->2u) of the summand note admits an integer-geometric object whose norm is n and whose lost coordinate is an angle, as Gaussian integers do for the multiplicative forest (no such object found, SCOPE); and any single object realizing both the 2cos(2pi/9) cubic 3-cycle and THM-4139's quarter-integer 3-cycle (none known; the two live in different fields). The next question worth a lane is the additive counterpart, not further censuses.
