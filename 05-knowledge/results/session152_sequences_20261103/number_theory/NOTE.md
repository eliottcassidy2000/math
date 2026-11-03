# A semigroup distinction for the two lists

Let A = {15,17,19,23,25,27,29,33,37,41} and B = {2,3,5,7,11,13,17,19,23,27}, taken exactly as given. In particular 27 = 3³ is composite, not a prime. Their intersection is {17,19,23,27}. The finite lists by themselves do not determine any continuation rule.

For a finite positive-integer set C, let S(C) be its additive numerical semigroup (including 0). Since 15 and 17 are coprime, S(A) is cofinite. Its Apéry set with respect to 15, indexed by residue r = 0,...,14, is

(0, 46, 17, 33, 19, 50, 36, 37, 23, 54, 25, 41, 27, 58, 29).

Indeed, every displayed number belongs to S(A); the only entries not already in A or 0 are 46=23+23, 50=25+25, 36=17+19, 54=27+27, and 58=29+29. To establish minimality, enumerate the sums of at most three terms of A below 59: four positive terms are at least 60, so no shorter competitor can have been missed. The verifier does this independently of its Dijkstra computation of the Apéry set.

Consequently, for every n >= 0, n belongs to S(A) iff n >= w_(n mod 15). Thus its Frobenius number is max(w)-15=43, its conductor is 44, and its number of gaps is sum_r (w_r-r)/15=26. Its univariate generating function is

H_A(x) = (sum_{r=0}^{14} x^{w_r})/(1-x^15).

All ten given members of A are minimal generators: each is odd and less than 45, whereas a sum of two members is even and a sum of three or more is at least 45. This gives an induction: the 15 Apéry representatives are base cases and adding 15 reaches exactly every other element of S(A). In particular, the 15 consecutive integers 44,...,58 are in S(A), proving by induction that every integer >=44 is in S(A), while 43 is not.

By contrast, S(B) = S({2,3}) = {0,2,3,4,...}. Indeed each listed element except 2,3 is at least 4 and can be formed from 2 and 3, and every n>=2 is so formed by induction on n using n-2. The minimal generating set is {2,3}, the Frobenius number is 1, the conductor is 2, and the genus is 1. Its generating function is

H_B(x) = 1+x²/(1-x) = (1-x^6)/((1-x²)(1-x³)).

The difference between the two generating functions is a structural fact about additive closure, not a claim that the original lists are prime sequences or that a hidden rule has been identified. `python3 verify.py` independently checks the finite certificate and arithmetic; the Apéry-residue induction proves the infinite assertions.
