# Square-sum four-cycles: exact activation and the 96 layer

Let Q(N) join distinct positive integers at most N whose sum is a square. The repository's [earlier four-cycle study](../crossroads_family_20260926_squares.md#3-square-labels-and-local-switches) derives the alternating sum and proves that the ordinary first four-cycle occurs at 46. The refinements here are the exact height rule, the repeated-label threshold, and the synchronized activation layers at 99 and 30,013.

## Exact heights of a square identity

Every simple four-cycle has cyclic edge sums (A,B,D,C) up to orientation and origin, where A < B ≤ C < D are squares and A + D = B + C. Opposite sums are equal in total because both cover the four vertices. Equal adjacent sums force a repeated vertex. Thus any repeated label is at the middle (B = C), not an extreme.

The vertices for the ordered pattern are (t, A−t, t+B−A, C−t). The exact height-N condition is:

    max(1,C−N) ≤ t ≤ min(A−1,N−B+A);
    2t not in {A, 2A−B, C, A+C−B}.

These exclusions cover all collisions between the two alternating vertex pairs; within either pair the differences B−A, C−A are positive. If B < C the orientation with edge sequence A,B,D,C is unique, hence one cycle per admissible t. If B = C, the reversal at the smallest edge identifies t with A−t; choose 2t < A.

In the repeated case put (A,B,D) = (a²,b²,c²) and d = B−A, so a²+c² = 2b². The canonical parameter spans 1 ≤ t ≤ floor((A−1)/2), with only 2t=A−d possibly excluded. Its height is B−t. The first height is **floor(c²/2)+1** whenever a ≥ 2: the largest allowed t has A−2t in {1,2}, whereas d ≥ 5. The case a=1 cannot split the smallest sum into two positive integers. For primitive odd-square progressions the collision 2t=A−d is impossible by parity, and the entire height interval [(c²+1)/2, b²−1] occurs.

Only roots c≤14 can give a repeated template by height 99; checking these roots gives just (1,5,7) and (2,10,14). The first fails positivity; the second first appears uniquely at **99** as (1,3,97,99). The height rule gives (24,25,144,145) at 145 and (24,25,264,265) at 265 from the root triples (7,13,17) and (7,17,23).

## All four cycles born at 99 are one factorization layer

The common difference is 96 = 2⁵·3 = 100−4 = 121−25 = 196−100. Factoring y²−x² sends the three same-parity factor pairs (8,12), (6,16), (4,24) to root pairs (2,10), (5,11), (10,14). The only other same-parity pair (2,48) is too high. Pairing any two of the three square differences gives these identities and cycles:

| Identity | New cycles |
| --- | --- |
| 4+121 = 25+100 | (1,3,22,99) |
| 4+196 = 100+100 | (1,3,97,99) |
| 25+196 = 100+121 | (3,22,99,97), (1,24,97,99) |

For an independent completeness check, a newly born cycle contains 99, whose possible neighbours are {1,22,45,70,97}. The only pairs with a second common neighbour below 99 are (1,22) via 3, (1,97) via 3 or 24, and (22,97) via 3. The full four-cycle count rises from 11 to 15. The occurrence of root 11 is part of the same factorization layer as the repeated template.

The repeated template supports an actual Hamiltonian two-switch at its first height. The literal cycle in the replay begins 1,3,…,99,97,…,48, then closes back to 1. Reversing the segment from 3 through 99 swaps {1–3,97–99} for {1–99,3–97}. The segment has 55 vertices, so every integer retains its position parity. This is a new repeated-label spanning-switch threshold, distinct from the ordinary threshold 46. Because one cycle loses edge 1–3, they are not a common 1-to-3 port pair; cutting a common edge such as 6–10 gives ports 6 and 10 instead.

## Pythagorean and golden-ratio coordinates

For a primitive triple of roots a²+c²=2b² all three roots are odd (reduce modulo 4 and exclude a common factor 2). Then (c+a)/2 and (c−a)/2 are the coprime legs of a right triangle of hypotenuse b. Conversely a primitive Pythagorean triangle with legs u,v gives a=|u−v|, c=u+v. Scaling multiplies square sums by a square, but can change positivity: the primitive 3–4–5 triangle gives unusable (1,5,7), and scaling its roots by 2 gives the first usable motif at 99.

For consecutive Fibonacci Euclid parameters m=F(k+1), n=F(k), k≥2, put p=mn and ε=(−1)ᵏ. Cassini gives m²−n²=p+ε and the other leg is 2p. The three roots before reducing common factors are (p−ε, F(2k+1), 3p+ε). The cases (m,n)=(2,1),(3,2),(5,3) give (1,5,7),(7,13,17),(7,17,23), the third after dividing out 2. As m/n tends to the golden ratio, the legs tend to ratio 1:2 and the square labels to ratio 1:5:9. Separately, the Fibonacci matrix Q=[[1,1],[1,0]] has tr(Q⁵)=11 and det(Q⁵)=−1, so F(k+10)=11F(k+5)+F(k). That is an exact recurrence involving 11, not a premise of the height cutoff.

## Synchronizing 7 and 49 at the exact target 30,013

For odd largest root c the two top vertices at first activation are (c²−1)/2 and (c²+1)/2. The odd root c=245 has exactly two positive repeated templates, with root triples

    (35,175,245)  = 35·(1,5,7),
    (155,205,245) =  5·(31,41,49).

The primitive right triangles are 3–4–5 and 9–40–41. At their first height **30,013**, the two cycles are (612,613,30012,30013) and (12012,12013,30012,30013). The six vertices induce exactly the union of the two four-cycles, a theta with three routes of lengths 1,3,3 between the same top vertices. Its two short exchanges share their top edge; simultaneous availability does not make them independent. The two scaled families synchronize whenever c is a positive multiple of 49, with first synchronized height 1201. For even c the two top vertices are c²/2±1, again common to both families.

More generally, for odd c factor its split-prime part in the real quadratic norm ring: let e(p) be the exponent of each prime p≡±1 mod 8 dividing c. The number of distinct root templates at c, before excluding the unusable a=1, is

    (product over split p of (1+2e(p)) − 1)/2.

For the primitive count, Euclid parameters m>n>0 give c=(m+n)²−2n², with gcd(m,n)=1. The ring Z[√2] is norm-Euclidean: rounding each coordinate gives errors of size ≤1/2, whose norm has absolute value ≤1/2. A primitive odd norm d>1 therefore uses only split primes (the criterion is the quadratic character of 2) and chooses one prime over each split factor, with its entire exponent. An inert prime cannot divide the norm without dividing both coordinates; selecting both conjugate primes at a split factor also violates primitivity. A generator can be made totally positive using the unit 1+√2 of norm −1 and an overall sign. The totally positive norm-one units are powers of (1+√2)², as follows by reducing into [1,(1+√2)²) and checking the two integer coefficients. Writing α=(m+n)+n√2, the condition 0<n/(m+n)<1/2 is equivalent to 1<α/√d<1+√2. Modulo positive units and conjugation, the interval gives exactly one representative; boundary points would require n=0 or m=n, impossible for a primitive odd norm d>1. Conjugation pairs opposite split choices without fixed points. Hence the primitive count is 2^(number of distinct prime factors of d − 1). Sum over the possible primitive roots d|c, multiplying the roots by c/d. This yields the product formula. To count physically activated templates, subtract one if there is an a=1 triple. For c=245 = 5·7² there are two templates; 5 is inert and appears only through scaling.

Run python replay_square_quartets.py. It checks the exact-height formula against an independent common-neighbour graph construction through 180 and verifies both 99-vertex Hamilton cycles and the 30,013 synchronization.
