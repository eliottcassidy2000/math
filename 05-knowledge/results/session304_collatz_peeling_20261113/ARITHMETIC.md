# Guarded ternary descent: Collatz and tournament peeling

## An exact arithmetic tree avoiding the tournament holes

For positive odd n let D_s(n) be the odd part of 3n+s, with s=±1. Put S equal to the positive odds with 7 and 21 removed. Every n>1 in S has a path in S of at most three such transitions to a smaller m, and m≤(253/299)n. Three steps and this uniform contraction constant are sharp.

Choose s congruent to n modulo 4, so 3n+s=2^a m with a≥2 and m<n. Every D_s value is coprime to 3, so m cannot be 21. If m is 7, write n=(7·2^a+ε)/3 with ε=−s. The opposite sign gives y=7·2^(a−1)+ε exactly (one halving). For a≥5 the same sign ε at the next two moves gives z=21·2^(a−3)+ε and the odd part of w=63·2^(a−5)+ε (the first of these takes two halvings, and w is the value after at least two halvings of the last branch). The intermediates avoid both holes. At a=5 the integral case has ε=1 and ends at odd part(64)=1. For a≥6, w is odd. The ratio w/n tends to 27/32 and is monotone towards that limit from above for ε=+1, below for ε=−1. The first a≥6 with ε=+1 is a=7, giving 253/299. Low cases a=2,3,4 are 9→13→5, 19→29→11, 37→55→41→31. Direct non-bypass contraction is at most 5/6. From 299 the only allowed first step is 449; the second outcomes are 337 and 673; third outcomes are 253,505,1009. This proves sharpness; 37 has no allowed descent in two odd moves. Choosing these parents gives a rooted tree with logarithmic macro depth. This is not ordinary Collatz: both forward signs are essential to the theorem as stated.

## Meaning of the lift obligation

A tournament triangle peel preserves the 3-free part and travels on one ×3 ray. Every D_s endpoint instead has valuation zero at 3, and the bypass moves between rays. It therefore cannot be simulated by any sequence of peels and inverse peels. A physical implementation must change the 3-free part and retain an appropriate marked interface. Bare arithmetic solvability or 3-adic coverage cannot discharge this obligation.

There is another conditional model. The ternary digits of (H−1)/2 are implemented by 3H−2, 3H, 3H+2 respectively. If all three operations were physically available on a common renewable ready class whenever their output avoids 7 and 21, only the seeds 1,19,23,61,63,65 would be required: each other allowed odd target has a unique smaller allowed odd parent by removing its last ternary digit. The five nontrivial seed roots are exactly the allowed children of 7 and 21. This is a proof target, not a statement that such a ready class exists.
