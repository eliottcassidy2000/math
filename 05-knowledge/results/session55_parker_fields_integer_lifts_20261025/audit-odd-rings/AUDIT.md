# Independent audit: complete quotient-ring Parker classification

25 October 2026. PASS, with explicitly identified dependency on the independently audited prime-field classification in `session55-finite-field-geometry`. This audit does not re-run that separate 1730-field census. All rings here are Z/nZ, all squarehood includes zero, and all nine entries must be distinct as residues. No integer-square solution or first-discovery claim follows.

## 1. Prime-power proof, including the small exceptional squares

Use the division-free magic plane M(c,u,v), with entries c plus the nine offsets su+tv for s,t in {-1,0,1}, in the usual magic arrangement. Its universality follows from total=3c and opposite pairs=2c. There is no assumption that 2 or 3 are invertible in this identity.

Let L be Lo Shu, with entries 1 through 9. The integer magic array 1+pL has all entries congruent to 1 modulo an odd prime p. Every such integer has a square root modulo every p^k by ordinary simple-root Hensel induction. Its entries are distinct modulo p^k exactly when p^(k-1)>8. Therefore all p>=11 have non-Parker quotient rings at every k>=2, and p=3,5,7 have them at every k>=3. This disposes of all higher prime powers without needing any lift of a field witness having a zero entry.

Modulo 3, opposite square entries summing to twice the square center must both equal the center: the only squares are 0 and 1. Every square magic array is constant. Modulo 9 there are just four square residues, too few for nine distinct entries.

Modulo 5, Q={0,1,-1}. If the center is nonzero, square scaling normalizes it to 1. The only pair of squares with sum 2 is (1,1), so the array is constant. If the center is 0, then u,v,u+v,u-v lie in {0,1,-1}. Both u,v nonzero would give u+v or u-v equal to 2 or -2, impossible. Thus one is zero. A nonconstant array has exactly three zero entries; a constant zero array has nine.

Modulo 7, Q={0,1,2,4}. If the center is zero, opposite entries sum to zero and -1 is not a square, so every entry is zero. Otherwise normalize the center to 1. The allowable offsets x with 1+x and 1-x both squares are precisely {0,1,-1}. Hence u,v,u+v,u-v all lie in that set, forcing one of u,v to be zero. Unless constant, the array has three entries equal to zero.

A square residue modulo p^2 that is zero modulo p equals zero modulo p^2. Thus the three zero entries in either nonconstant reduction at p=5 or 7 remain equal upstairs. A constant nonzero reduction allows at most p possible lifted entry residues, fewer than nine. A constant zero reduction is all zero upstairs. This proves 25 and 49 Parker, not merely a census observation.

The known prime-field Parker list is 3,5,7,11,13,17,19,23,31,43,47,67. Combined with these elementary arguments, the complete odd prime-power list is that prime list together with 9,25,49. The finite-field input is used only at exponent one.

## 2. CRT construction and odd classification

The three integer squares 1,25,49 form an arithmetic progression. Their differences are 24,24,48, none divisible by any odd m>=5. Place y,z,x / x,y,z / z,x,y, with (x,y,z)=(1,25,49). All lines sum to 3y. Its equal-cell partition in row-major labels is

A={{0,4,8},{2,3,7},{1,5,6}}.

Reflect vertically to obtain

B={{2,4,6},{0,5,7},{1,3,8}}.

Every A block meets every B block in exactly one cell. Over coprime odd factors a,b>=5, CRT-combine these two arrays at fixed cell labels. Equality globally means equality in both components, so all nine combined entries are distinct. Square roots CRT-combine as well. This proves every such product non-Parker.

A non-Parker witness on a complete prime-power component extends across the coprime complementary factor by combining with a constant-zero array. This assertion is not arbitrary divisibility monotonicity: the factors must be coprime.

The mod-3 component of any square magic array is constant. Hence for gcd(3,m)=1, an array modulo 3m is distinct exactly when its component modulo m is distinct; a witness modulo m extends by a constant modulo 3. Thus Parker(3m) iff Parker(m).

If an odd n has two distinct primes >=5, the orthogonal-partition construction applies. If 9 divides n and a prime >=5 also divides n, it applies with factors 3^k and the remaining odd factor. If n is a pure 3-power, only 3,9 survive. All remaining cases are p^k or 3p^k with p>=5, whose classification was just proved. Therefore there are exactly 28 odd Parker moduli n>=2, with largest 201. Modulus 1 is also Parker if the zero ring is admitted, and is included explicitly in the computational convention.

## 3. Complete all-modulus finite reduction

The integer magic array 1+8L has entries 1 modulo 8. Each is a square modulo every 2^k, and its nine entries are distinct modulo 128 and every higher power. Therefore every Parker n has 2-adic exponent <=6. Also its full odd part must be Parker: any non-Parker odd component combines with a constant array modulo its coprime dyadic component.

Thus all Parker n belong to exactly 7*29=203 candidates 2^a*m, with 0<=a<=6 and m=1 or one of the 28 classified odd moduli. The largest candidate is 64*201, but enumeration of that finite universe gives largest actual Parker modulus 3216.

## 4. Independent finite certificate

`audit.py` does not import any producer code. It enumerates arrays using top-left, top-right and center square entries, reconstructing the remaining entries from opposite pairs and line sums. This is a complete enumeration by Section 1. It enumerates each FULL odd candidate modulus, whereas the producer factors those odd moduli into prime-power mask tables. This is an independent CRT aggregation cross-check, not merely the same table reuse.

It computes full fixed-label equality masks, inclusion-minimal masks, and all 203 CRT combinations against the dyadic tables. Removing a mask that contains another realizable mask is safe because intersection is monotone. Each positive CRT witness is checked directly for all eight line sums, nine distinct entries, and square membership at the full candidate modulus.

The independent output agrees with every producer local prime-power antichain and all 203 statuses. Exactly 137 candidates are Parker including modulus 1, or 136 for n>=2. The largest is 3216. Results and independent tables are in `audit_results.json`. Producer certificate checked: `session55-crt-collision-wiring/classification_candidates.json`.

## 5. Literature boundary

Cain's primary source https://arxiv.org/html/1908.03236v2 , Section 9, explicitly labels the odd maximum 201 as Conjecture 9.3 and the unrestricted maximum 3216 as Conjecture 9.4. The present proof reduction plus certified finite exhaustion proves those statements relative to the independently proved prime-field classification; this is not a claim to have solved the integer problem.

The author-maintained repository https://github.com/jiji7879/magic_squares_fields advertises a finite-FIELD classification and arbitrary powers of entries, not a classification of quotient rings Z/nZ. Its contents cannot be treated as establishing priority or absence of priority for this ring result. Targeted searches found no additional primary ring-classification proof. That limited search is not a novelty certificate.
