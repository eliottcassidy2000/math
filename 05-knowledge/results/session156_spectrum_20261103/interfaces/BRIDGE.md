# An exact 7p obstruction in a strong feedback-two clone family

Source: S155 `05-knowledge/results/session155_spectrum_20261103/primes/README.md` at commit `ea9db36824ee1b3c0730a3c5325826b79a8f3265` supplies both the seven-vertex H=35 witness and the all-k clone-slot identity. This note specializes those existing results and adds an elementary factorization and valuation analysis. No literature-priority claim is made.

Let T have out-neighbor bitmasks `(62,60,120,16,96,8,43)` on vertices 0,...,6. Replace vertex 0 with k >= 1 transitive clones, all having the external orientations of vertex 0. Then the resulting tournament T_k is strong and

    H(T_k) = 18*2^k - 1.

Indeed 0,1,2,5,3,4,6,0 is a directed Hamilton cycle in T, and substituting the ordered clones for 0 gives one in T_k. For the clone-slot formula, count each permutation of 1,...,6 by the number t of valid slots and q of backward adjacent pairs, discarding any permutation with an invalid defect. The only surviving (t,q) multiplicities are

    (1,0): 3, (2,0): 14, (2,1): 4.

They contribute 3 + 14*2^k + 4*(2^k-1), proving the formula. A full standalone exact enumeration is in verify.py, which also independently checks this identity by vertex-subset DP for k=1,...,8. The base T has feedback vertex number exactly two by the cited S155 source; as T_k contains an induced T, and deletion of vertices 5 and 6 leaves a transitive tournament, T_k has feedback vertex number exactly two for all k.

**7p restriction.** Since 18=4 modulo 7 and 2 has order 3 modulo 7,

    7 divides H(T_k) iff k = 1 mod 3.

If k = 1 mod 6 and k >= 7 then

    H(T_k) = (3*2^((k+1)/2)-1)(3*2^((k+1)/2)+1).

The factors are odd, exceed 7, and are relatively prime because they differ by 2. Exactly one is divisible by 7; after removing one factor of 7, both remaining factors exceed 1. Therefore H(T_k)/7 is composite. For k=1, H(T_1)=35=7*5.

If k=1+3t, reduction modulo 49 uses 2^(3t)=(1+7)^t=1+7t (mod 49):

    H(T_k) = 36*8^t-1 = 7*(5+t) (mod 49).

Therefore 49 divides H(T_k) iff k = 7 mod 21 (among the 7-divisible cases). Combining the two restrictions, for prime p and H(T_k)=7p, necessarily k=1, or k=4 mod 6 and k is not 28 mod 42. These conditions are necessary, not claimed sufficient. Explicit examples are k=4: 287=7*41, and k=10: 18431=7*2633; 41 and 2633 are prime by trial division. No inference to any other tournament family, or to infinitude of prime quotients, is made.

Run `python3 verify.py`. The mod 49 condition and the factorization are proofs for every k; the diagnostic scan is finite and is not used as their proof.
