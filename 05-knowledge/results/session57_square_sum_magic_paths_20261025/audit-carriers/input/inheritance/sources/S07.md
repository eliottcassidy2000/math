# Session 55: finite-field square residues versus no-wrap integer lifts

Status: PROVED reductions and completeness statements; FINITE-EXACT bounded certificates. The existence and least size of a 3-by-3 magic square of nine distinct integer squares remain OPEN. No computational record or novel classification of the classical problem is claimed.

## 1. Main findings

1. The proposed no-wrap route is valid **if the actual nonnegative entries being added are integer squares**. Quadratic residues alone do not suffice.
2. The smallest prime modulus admitting a canonical, distinct, no-wrap magic square of quadratic residues is **59**, in the finite-exact sense certified here. Its entries are

   28   3  17
    5  16  27
   15  29   4

   Every line sums to 48, below 59. Every entry is a square modulo 59. Only 16 and 4 are squares as integers. Even requiring the center to be an integer square does not repair the implication.
3. A particularly transparent hostile is the Lo Shu square

    8  1  6
    3  5  7
    4  9  2

   modulo 311: all entries are quadratic residues and every line sums to 15, but only 1, 4, 9 are integer squares.
4. A complete center-first integer search is simpler than scanning and rescaling finite-field solutions. It searches a four-gap additive pattern in the exact opposite-square pairs. The supplied run excludes center roots 1 through 100,000, with no primitivity filter. This is a reproducibility demonstration, far below older published search bounds.

## 2. The exact no-wrap lifting lemma (PROVED)

Let n>=2 and let A be a 3-by-3 array of nonnegative integers. Assume its eight row, column and diagonal sums are congruent modulo n. If all these sums lie in an interval of length strictly less than n, they are equal as integers: any difference is a multiple of n of absolute value less than n.

In particular, sums in [0,n) suffice. Sums in (0,n] suffice as well. Thus the user's weak upper bound “never exceeds the modulus” works for nine distinct nonnegative entries: every line sum is positive, so the second interval applies. Without positivity, the closed interval [0,n] does not itself make reduction injective.

For a magic square **of integer squares**, additionally verify all nine equations A_ij=r_ij^2 for nonnegative integer roots r_ij, and verify distinctness of the values A_ij, not merely distinctness of chosen signed roots. These conditions, together with the lemma, are sufficient. For canonical representatives A_ij in [0,n), the assertion that A_ij is a square modulo n is strictly weaker than A_ij being an integer square.

The two hostile examples above satisfy every residue, distinctness, sum and size test. They fail precisely the integer-square test. Squaring a modular square root without reduction repairs that entry's squarehood but changes its magnitude and generally changes line sums by different multiples of n. A separate carry sidecar is indispensable.

This argument applies to Z/nZ, hence to the prime field F_p. It does **not** identify the extension field F_(p^r) with Z/(p^r)Z for r>1. Integer reduction into that field lands in its prime subfield and has characteristic p. Numbering its p^r elements does not produce compatible integer arithmetic or a modulus-p^r lift.

## 3. Why the exact version is as hard as the integer problem (PROVED)

A no-wrap prime-field witness whose representatives are nine distinct integer squares is already an integer solution by Section 2. Conversely, given any integer solution with magic sum T, take a prime p>T. Every entry is nonnegative and at most T, so its reduction is distinct and its actual line sums remain T<p. Thus the existence of a genuine no-wrap prime-field witness is equivalent to existence of an integer solution. No field-first search removes the integer-square existence obstruction.

Every integer 3-by-3 magic array with center c has sum 3c and opposite pairs summing to 2c. One proof sums the four center-crossing lines: this equals the sum of all entries plus 3c; the three row sums give the sum of all entries, hence 4T=3T+3c.

Consequently an integer-square witness has c=e^2 and T=3e^2. A nonnegative such array with c=0 is all zero, not distinct. It cannot contain a zero elsewhere either: the opposite entry would be 2e^2, which is not a square for e>0. Every genuine distinct nonnegative witness is therefore positive.

In the exact prime no-wrap problem, 3e^2<=p is automatically strict: equality with prime p would force p=3,e=1, which cannot accommodate nine distinct square entries. A complete negative result through e<=B therefore excludes genuine no-wrap witnesses at every prime p<=3B^2. This is a **bounded exclusion**, not a bound on where a solution must exist.

## 4. Complete center/gap carrier (PROVED)

For e>=1 define

D_e={d>0 : e^2-d and e^2+d are integer squares}.

Every such d is strictly below e^2. All integer magic arrays with center c have the form

 c+x    c-x-y    c+y
 c-x+y  c        c+x-y
 c-y    c+x+y    c-x.

This follows by choosing the two top corners and using opposite-pair sums and the top row/left column. Conversely, direct addition gives all eight sums 3c.

Rotations and reflections act by signed permutations on (x,y). If all nine entries are distinct, the two corner deviations are nonzero and have different absolute values. A unique representative has x>y>0. The nine deviations are exactly 0, ±x, ±y, ±(x+y), ±(x-y).

It follows that a distinct integer magic square of squares with center e^2 exists if and only if there are x>y>0 satisfying

x,y,x+y,x-y in D_e, with four different values.

The last condition is equivalently x!=2y given x>y>0. This is an exact iff, not just an obstruction. The code retains the explicit cardinality test. It checks every unordered pair of distinct positive gaps and uses hash membership for the sum and difference. A center with fewer than four gaps is safely discarded.

### Complete Euclidean sieve

Every positive gap is generated by

e=k(m^2+n^2), d=4mn(m^2-n^2)k^2,

where m>n>=1, gcd(m,n)=1, m-n is odd, and k>=1. Restricting e<=B therefore gives finite loops m<=sqrt(B), n<m, k<=B/(m^2+n^2).

Proof of completeness: write a^2+b^2=2e^2 with 0<a<e<b, divide the common gcd of a,b,e, and obtain a primitive triple. Its roots a,b,e are odd: if e is even then a,b are even, contrary to primitivity; if e is odd their squares sum to 2 modulo 4. Set U=(b+a)/2 and V=(b-a)/2. These are coprime positive integers of opposite parity and U^2+V^2=e^2. The primitive Pythagorean parametrization gives {U,V}={m^2-n^2,2mn}. Restoring the common gcd k gives d=(b^2-a^2)/2=2UV=4mn(m^2-n^2)k^2. Conversely these formulas directly produce squared endpoints. No common-factor assumption on the full nine-cell square, or on different opposite pairs, is introduced.

The sieve is therefore complete, including nonprimitive arrays. It has O(B log B) generated-gap work as an upper bound, since summing B/(m^2+n^2) over the parameter lattice is O(B log B). The additive test costs O(sum_e |D_e|^2), with constant-time expected set lookups. This is an algorithmic improvement over separately enumerating many moduli, not a claim of improving the known research record.

### Independent count check

For e=2^a product(p_i^alpha_i) product(q_j^beta_j), with p_i=1 mod4 and q_j=3 mod4, the standard sum-of-two-squares representation count applied to 2e^2 gives

|D_e|=(product_i(2alpha_i+1)-1)/2.

Indeed, r_2(2e^2)=4 product_i(2alpha_i+1). Four solutions have |a|=|b|=e; every nontrivial positive unordered pair supplies eight signed ordered representations. Subtract four and divide by eight. The code obtains this count independently by a smallest-prime-factor sieve and compares it at **every** searched center; it separately checks each generated gap's two integer square roots. Direct square-root enumeration also checks all centers through 2,000.

## 5. Exact computational packet

Run from this directory:

    python3 verify_and_search.py --bound 100000 --independent-bound 2000
    python3 verify_and_search.py --bound 100000 --independent-bound 2000 --verify
    python3 test_verifier.py

The third command runs both ordinary and optimized Python, requires identical output, and exercises actual optimized CLI rejection of changed content, missing keys, extra keys and a missing certificate file; results are saved in verifier_checks.json.

The second command recomputes everything and compares the entire deterministic JSON object with certificates.json. All checks use explicit RuntimeError guards, not optimization-removable assert statements. The verifier internally tests rejection of changed content, a missing key, and an extra key. Standard-library Python, integer arithmetic only. No probabilistic primality or floating-point square test is used.

Recorded universe and results:

- Center roots: all integers 1<=e<=100000; arbitrary nonprimitive arrays included.
- Primitive Pythagorean parameters: 15,919.
- Generated gap insertions: 161,436.
- Pair tests at centers having at least four gaps: 347,561.
- Maximum gap count: 40; first center attaining it: 32,045.
- Genuine solutions: zero.
- Consequent no-wrap exclusion: every prime modulus p<=30,000,000,000.
- The serialized complete gap-bank digest is in certificates.json.

For the weak predicate too, allowing sums equal to p makes no difference: T=3c=p would force p=3,c=1, whereas nine distinct nonnegative entries have total at least 0+1+...+8=36, greater than 3T=9. Thus the strict and weak bounds have exactly the same witnesses for prime moduli.

The weak residue search enumerates all primes through 59, all admissible centers c with 0<3c<p, all d with c±d quadratic residues, and every required pair of deviations. Each resulting full canonical orbit set is independently compared with direct enumeration of the center and two ordered corner entries from the quadratic-residue set. Earlier primes all return no weak witness; p=59 is the first. This establishes minimality **only for this weakened prime-field predicate**, not for genuine integer-square witnesses, and not for fields of extension degree greater than one.

Controls also include Cain's F_29 example, whose canonical line sums wrap and are not all equal; the classical Parker near miss, whose actual entries are squares but are not all distinct or magic; an all-4 square, which passes squarehood and magic but fails distinctness; and synthetic gap set {1,2,3,4}, which correctly produces Lo Shu.

## 6. Research inheritance and connection ledger

Anchor: exact no-wrap integer-square lifting. Niche: a minimal weak-predicate hostile. Wildcard: replace field-square orbits by a finite additive-pattern search on exact opposite-square gaps.

Concept board: (1) modular magic equations; (2) integer squarehood; (3) archimedean height/carry; (4) dihedral orbit normalization; (5) Pythagorean gap generation. Each certificate separates these predicates. The prior general guardrail “a quotient owes a loss ledger” is instantiated here: reduction preserves equations and square-residue status, but loses integer magnitude, squarehood of canonical representatives and carries. The restoring sidecar is an exact root for every chosen entry plus the eight integer sums. Lo Shu is the decisive hostile. No artificial tournament relation is imposed.

The closest classical mechanism is opposite-square arithmetic progressions; the corrected near miss is Parker's actual-square array; the least-used useful sidecar is the eight unreduced line sums. Repository startup found no dedicated Parker route. The parent owns current synchronization, source routing and publication; this packet modifies no existing repository files.

## 7. Literature and attribution boundary

Primary source: O. M. Cain, *Gaussian Integers, Rings, Finite Fields, and the Magic Square of Squares*, arXiv:1908.03236v2, sections 5–6, https://arxiv.org/html/1908.03236v2 . Cain uses nine distinct squared field elements, proves F_29 the smallest non-Parker finite field, and gives a field search normalized by multiplicative scaling. This packet independently proves its integer-lift criteria and search completeness; it does not import the paper's algorithm as a completeness certificate.

Historical calibration: Pierrat, Thiriet and Zimmermann, *Magic Squares of Squares*, dated March 13, 2015, https://members.loria.fr/PZimmermann/papers/squares.pdf , introductory paragraph and section 2. They report a much larger earlier center-cell bound and discuss the need not to impose coprimality on every arithmetic progression. Our bounded run is only a reproducibility/control exercise; it is not a new lower-bound record. No assertion about the current strongest bound is made.

Copyright-safe source extracts and exact locations are in sources.json. All longer mathematical arguments above are independently written.
