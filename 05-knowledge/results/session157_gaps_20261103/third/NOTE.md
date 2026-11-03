# A one-apex family covering every residue modulo 75

Status: PROVED. This supplies a third strong family and closes the modulo-75 gap of the two previously investigated families. It is a congruence result, not a claim about all exact odd Hamilton-path counts. No novelty claim is made.

For n>=1, take a transitive base with orientation word 1(001)^n0, where 1 means the added apex r beats that base vertex, and 0 means the base vertex beats r. Orient base vertices in their listed order. The resulting tournament is strong: r beats the first base vertex, which reaches every later base vertex; every base vertex reaches the last base vertex, which beats r.

Its Hamilton-path count is

    F(n) = (80*8^n + 28n + 67)/49.

For a direct count, delete r and split the transitive base into paths before and after r. For a 0 at position i and an earlier 1 at j, the choice of the intermediate vertices contributes 2^(i-j-1). A 0 followed immediately by a 1 contributes one further path, and the two exterior placements contribute two paths. Let B_i be the sum of 2^(i-j-1) over earlier 1s j, immediately before processing position i. Then B updates to 2B+1 at a 1 and 2B at a 0, and a 0 contributes the current B. After the initial 1, B_0=1. Each block 001 maps B to 8B+1 and contributes 3B+1, including its 01 adjacency. The final 0 contributes B_n. Thus

    B_n = (8*8^n-1)/7,
    F(n) = 2 + n + 3*sum_{j=0}^{n-1} B_j + B_n,

which simplifies to the formula. This counts all paths, independently of finite experiments.

Every residue modulo 75 is attained by F(4t), t>=1. For p=3 or 5 and k>=1, 8^(4*p^(k-1))=1 mod p^k by the binomial theorem, since 8^4=1 mod 15. Hence

    F(4(t+p^(k-1))) - F(4t)
       = (16/7) p^(k-1) (mod p^k),

and 16/7 is a unit modulo p. It follows inductively that t -> F(4t) is bijective modulo 3 and modulo 25. Therefore CRT chooses t modulo 75 that attains any desired combined residue. The dependency of the 3- and 5-adic conditions is fully accounted for by using this same t; 4t remains positive for positive t. In particular the eight classes 9,12,24,27,54,57,69,72 missed by the two earlier families are all attained. This conclusion does not show that any particular integer in those classes is itself a path count.

The same reasoning, with primes processed in increasing order and L_p=ord_p(8), proves every residue modulo any odd M coprime to 7. At each prime p dividing M, L_p divides p-1, so its prime factors have already been assigned; holding n mod L_p fixed gives a unit lift because the linear coefficient 4/7 is invertible modulo p. The 2-adic coordinates appearing in these orders are fixed arbitrarily before the odd primes. This extension needs no 7-adic claim; when 7 divides M, the denominator 49 requires separate analysis. Indeed F(n)=3 mod 4 for all n>=1, so this family alone does not cover all odd classes modulo 4.

Run `python3 check.py` and `python3 -O check.py`. Independent endpoint/subset DP checks n=1..4; the second control verifies all 75 classes among F(4),...,F(300), including the eight target classes. Finite controls complement, but do not replace, the proofs.

## A companion that also handles powers of 2

Replace the final 0 by 00, making the orientation word 1(001)^n00. The same split count gives

    F_2(n) = (192*8^n + 28n + 53)/49.

Indeed, after the n blocks the two trailing 0s contribute B_n+2B_n=3B_n instead of B_n. This family is also strong. For any M not divisible by 7, the union of the two families realizes every eligible residue: all residues if M is odd, and all odd residues if M is even.

For the 2-part 2^e, first choose F or F_2 according to the target modulo 4: respectively 3 or 1. For e=1 either is fine. For large n, 8^n=0 mod 2^e, so F or F_2 equals (28n+c)/49 mod 2^e, with c=67 or 53. The numerator 28 has 2-adic valuation exactly 2. Thus, for e>=3, varying n mod 2^(e-2) gives every residue in the chosen class mod 4. Fix this coordinate first, and extend it to the precision needed by orders of 8 modulo the odd prime factors of M. Process the odd primes in increasing order as before. For each p, the required n mod ord_p(8) has already been fixed and the coefficient 4/7 is a p-adic unit. CRT assembles the independently selected coordinates, and n may be increased while preserving them. Thus the two choices genuinely supply simultaneous coverage, not just separate prime-power coverage. No coverage claim is made here when 7 divides M.

`python3 check_companion.py` independently checks the companion graph by subset DP and exercises a constructive version of the ordered CRT for all eligible residues modulo M<=250 with 7 not dividing M, as well as higher prime-power and mixed moduli. The same checks run with `python3 -O`.
