# A strong two-apex family covering every 3-power residue

Status: PROVED by an exact finite-state recurrence certificate and p-adic lifting. This is not a claim of exact integer-spectrum coverage or of simultaneous coverage modulo every composite modulus. No literature novelty is claimed.

Take apexes 0 and 1 with 1->0, and base vertices x_1,...,x_(2n+2) in transitive order. For n>=1 use the type word (0,1)^n 2,1, where bit a of a type means that the base vertex beats apex a. The tournament is strong: 0->x_(2n+1)->1->x_(2n+2)->0 is a directed cycle; every other base vertex reaches x_(2n+1), while 1 reaches x_1 and x_1 reaches all later base vertices. Deleting both apexes leaves a transitive tournament.

Its exact Hamilton-path count g(n) satisfies

    225 g(n) = 25 + 150n + (668+60n)4^n + 432*9^n.

The accompanying standard-library checker is an exact proof certificate of the identity. It implements both orders of the two-apex gap automaton from S156, with 27 states per order. Let P be the 54-state transfer for the two-symbol word 0,1, let v be the initial state after one P, and let lambda append 2,1 and accept the terminal paths. Write

    q(X) = (X-1)^2 (X-4)^2 (X-9).

The checker verifies lambda(P^j q(P)v)=0 for every j=0,...,53, in exact integer arithmetic. Cayley-Hamilton applied to the 54x54 integer matrix P makes all subsequent terms zero as well. Hence g(n), n>=1, obeys the scalar recurrence q; the displayed formula obeys the same recurrence, and the checker compares its first five values to the independently implemented automaton. Separately it builds the tournaments and checks strongness and their Hamilton paths by subset/endpoint dynamic programming for n=1,2,3. This is a finite, reproducible recurrence proof for all n, not extrapolation from five values.

For each e>=1, g(n) attains every residue modulo 3^e, with n>=1 and the tournament strong. To see this, take n>=e+2. Because 225=9*25, the formula and 4^n=(1+3)^n give

    g(n) = 2-n (mod 3).

For 1<=k<e, put d=3^k. Since (4^d-1)/d=21 mod 27 and 4^d=1 mod 9, while 432*9^n and its change are divisible by 3^(k+3), the numerator of 225(g(n+d)-g(n)), divided by d and reduced mod 27, is

    150 + ((668+60n)*21+60)*4^n
      = 15 + (21+18n)*4^n = 9 (mod 27).

The last equality follows from 4^n=1+3n+9 binom(n,2) mod 27. Consequently

    (g(n+3^k)-g(n))/3^k = 1 (mod 3).

The base map n -> 2-n is bijective mod 3. Each additional base-3 digit of n then selects the next output digit bijectively, proving all residues mod 3^e. The 3-power valuation calculation is used only in the stated range n>=e+2, and no false claim of periodicity at small n is needed.

Combined with the strong one-apex family (10*4^n+6n-1)/9, this proves strong FVS-at-most-two tournament counts cover every residue modulo p^e for every odd prime p and all odd residues modulo 2^e. It does not follow that either single family or their union covers every eligible residue modulo all composite moduli. For instance, this family misses precisely 2, 8, 11, 14, 17, 20 modulo 21: the checker verifies an exact eventual period of 630 from the formula and checks its complete period. The first family cannot give any multiple of 3. These are limitations of the individual families, not of their union: together the two families cover all residues modulo 21. More decisively, modulo 75 the second family misses the eight classes 9, 12, 24, 27, 54, 57, 69, 72, all divisible by 3; the first family misses them as well. The checker proves this by an exact eventual period of 2250. These are obstructions to the two-family construction only, not to strong tournaments in general.

Reproduce: `python3 check.py` and `python3 -O check.py`. The final loop checks all 3-adic residues for e<=7 as a hostile finite control, separate from the proof.
