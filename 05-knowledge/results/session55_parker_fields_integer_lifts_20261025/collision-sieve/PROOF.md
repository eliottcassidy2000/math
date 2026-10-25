# Exact center-independent Parker collision sieve

Status: proved necessary condition, with exact finite controls and a bounded performance experiment. No integer existence conclusion, no improved search-bound record, and no literature novelty claim.

## 1. Statement and discriminant route

Write an integer magic square as

    c+u    c-u-v    c+v
    c-u+v  c        c+u-v
    c-v    c+u+v    c-u.

Its nine offsets are precisely su+tv for s,t in {-1,0,1}. If every entry is a square, its reduction modulo a Parker prime is a square-entry magic array. By the definition of a Parker field it cannot have nine distinct reduced entries, even if the original integer entries were distinct. A collision is equivalent, for odd p, to

    C(u,v) = uv(u²-v²)(u²-4v²)(4u²-v²) = 0 mod p.

Indeed the difference of two grid indices has coefficients in {-2,-1,0,1,2}. Removing common factors of 2 yields the eight lines u=0, v=0, u=±v, u=±2v, 2u=±v. Every such line really is a difference of two offsets. This argument uses only a collision, not the center, and remains true when p divides c or either offset.

For an equivalent discriminant interpretation let F(T)=product_(s,t)(T-c-su-tv). Direct pair grouping gives

    disc(F) = 2^16 u^18 v^18 (u²-v²)^10 (u²-4v²)^4 (4u²-v²)^4.

There are nine pairs in each axis direction, five in each diagonal direction, and two in each remaining direction; eight differences have an extra factor of 2. The odd-characteristic zero set of this discriminant is exactly C=0. Thus this is a collision/discriminant sieve, not an independent squarehood certificate.

The proved classification supplied by session55-finite-field-geometry gives the odd Parker primes

    3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 47, 67.

Its proof and census, rather than a guessed list or statistical model, are the upstream classification dependency. Extension Parker fields give no new integer reduction primes: integers map into their prime subfields, not arbitrary extension elements. In particular F_27 is not Z/27Z.

## 2. Exact counts and CRT, with the content qualification

For p>=7 these are eight distinct projective lines. A line has p-1 nonzero affine points and all lines share only the origin, hence exactly 8p-7 affine pairs pass. For p=3 the lines merge to all four projective directions; for p=5 they merge to all six. Therefore all affine pairs pass for p=3,5,7.

For a squarefree product M of selected Parker primes, CRT is a bijection of the full residue-pair universe with the product of local universes. The exact passing fraction is the product of (1+L_p(p-1))/p², where L_3=4, L_5=6 and L_p=8 otherwise. For the complete list:

    M = 468229211845395
    affine fraction = 18129782170164195 / 874207008934460827483.

Condition instead on gcd(u,v,M)=1, meaning the pair is not (0,0) at any selected prime. The exact fraction is product L_p/(p+1), since each projective direction has p-1 representatives. For the complete list this is 32/1590435. The same formula applies to the full product of projective lines P^1(F_p). These are exact uniform finite-residue densities, NOT predictions for gap-bank pairs or arbitrary bounded primitive integer boxes.

Content must not be divided away blindly. Put g=gcd(u,v), u=gU,v=gV. C has degree eight, so C(u,v)=g^8 C(U,V). At primes not dividing g the primitive direction must be one of the collision directions; at primes dividing g the original test is automatic. The example (u,v)=(11,33) passes mod 11 but the divided pair (1,3) fails. It is a test of unsafe normalization, not an integer-square witness. If a *whole square array* is divided by a proven common square factor and actual roots remain integral, the normalized array may be tested afresh. Dividing only offset content does not establish this premise.

For a fixed g, the primitive-projective constraint is therefore imposed only at selected primes not dividing g. There is no missing primitive integer solution caused by a gcd filter in this experiment: none is imposed.

## 3. Conditional exact-gap rule

For a square center c=e², an exact opposite-square gap d satisfies c±d square modulo every prime. If p does not divide e, square scaling sends c to 1 and d to d/e²; if p divides e the center is 0. Hence the two finite sets

    G_p(c) = {d in F_p : c-d and c+d are squares}, c in {0,1},

cover every exact gap bank modulo p. The certificate enumerates both sets and all ordered pairs for every prime. For p=11, the square set is {0,1,3,4,5,9}, G_11(0)={0}, and G_11(1)={0,3,-3}. Every pair in each gap set has a zero coordinate or u=±v, hence is already a collision pair. Thus 11 adds no rejection once the opposite-gap condition is enforced. This is proved by complete 11-element enumeration; it is not inferred solely from the bounded integer run. Primes 3,5,7 were already vacuous. Consequently the reduced mask test may safely skip all four primes for exact gap-bank candidates.

## 4. Exact integer sidecar and bounded universe

For each 1<=e<=100000 the program enumerates every d>0 for which e²-d and e²+d are integer squares. Parameterization is

    e=k(m²+n²), d=4k²mn(m²-n²),
    m>n>0, gcd(m,n)=1, m-n odd, k>=1.

Completeness: writing a²+b²=2e² with a<b gives integers A=(a+b)/2 and B=(b-a)/2 satisfying A²+B²=e². Their primitive Pythagorean parameterization gives the formula. The excluded trivial d=0 cannot contribute nine distinct entries. The possible a=0 would imply b²=2e² and has no positive integer solution. All generated roots are checked with integer isqrt; brute root enumeration independently checks banks through e=1000. The complete serialized bank hash is compared with the upstream producer certificate by test_controls.py.

At each center with at least four gaps, enumerate every pair u>v>0. Four distinct positive gap magnitudes u,v,u-v,u+v are necessary and sufficient for nine distinct entries, and u±v must be in the same exact bank. The finite universe is exactly the upstream 347561 pairs. Keep original integer entries and roots, not least field representatives. Survivor samples store all nine integer entries, exact roots where they exist (null otherwise), gcd(u,v), and Euclidean (quotient,remainder) pairs mod 67. These quotients are the discarded carries; a survivor with null roots is NOT an integer-square magic square. The magic equalities alone hold identically even for signed entries.

## 5. Reproducibility and benchmark interpretation

Run `python verify.py` to rebuild certificate.json and benchmark.json. Run `python verify.py --verify` to recompute and compare the complete deterministic certificate. Timing samples are deliberately excluded from equality comparison. Run `python test_controls.py` for independent formula, bank-hash and mutation controls.

All computations use Python standard-library exact integers. Three wall-clock samples per algorithm measure filtering already-built banks/pairs, excluding bank construction, global lookup-mask construction, and certificate generation. Per-bank signature construction is charged to the signature algorithm. Benchmark numbers are observational, machine/load dependent, and not complexity guarantees. Both a conventional distinctness-first full-gap baseline and a cheaper difference-first full-gap baseline are included; any benefit from test ordering must not be attributed to the sieve when no candidate reaches it.

## 6. Inheritance and novelty boundary

Repository search examined magic/Parker text and collision/discriminant matches in math/01-canon, math/05-knowledge and math/07-reflections. Broad matches overwhelmingly concern unrelated topics; no dedicated pre-session55 integer Parker collision sieve was found in the searched local snapshot. This is not a global literature search and proves no priority. session55-finite-field-geometry already explicitly derives the eight collision directions; session55-parker-lifts already constructs the exact gap-bank universe and reports zero solutions. Those are inherited, not new here. This packet contributes the explicit integer necessary sieve, exact CRT/content bookkeeping, conditional p=11 redundancy, and measured algorithm comparison, including its negative performance result.
