# Two-burst tournament arithmetic: an exact inverse and a sharp thin-image law

Session241, 2026-11-10. **PROVED** statements below are elementary all-order statements about a specified tournament language. Computation is a separate **FINITE-EXACT** check. No literature-priority claim and no general Chan–Pak coverage claim are made.

## Inheritance, concept board, and scope

The closest mechanism is the physical O/Z state calculus in THM-4888, `01-canon/theorems/THM-4888-finite-oz-banks-mod3-law-involution-rigidity-residue-confinement.md`. Its coordinates are `(D,E,U,d)` and moves

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d).

The initial rooted edge has state `(1,0,2,1)` and D is its Hamilton-path count. The root-change extension J and its modular confinement, product identities and local factors are inherited from THM-4916 (`root-change-local-factor-identities-confinement-multipliers-product-family`); those facts do not imply integer coverage. THM-4861 (`the-five-fermat-primes-are-hamiltonian-path-counts...`) supplies Fermat counts but excludes prime counts for nontrivial vertex-transitive tournaments other than C3. THM-4986 (`the-7-21-omission-proof-has-an-exact-spanning-tree-twin...`) and THM-1370 supply the all-order omission of 7 and 21. THM-4885 (`representation mass of the one-apex join families`) records coverage through 2^35 only as finite evidence, not the conjecture.

The live board is: (1) shifted products, (2) binary runs, (3) multiplicative orders and prime valuations, (4) source-language coverage, (5) literal Hamilton paths. The hostile is “complete modular information implies global coverage”; the corrected near miss is the height-weighted/letter-weighted distinction explicitly retained in THM-4916. The least-used sidecar here is the binary run structure of the *shifted* count rather than the residue of H itself. Some historical source directories cited by the canon are absent in this checkout; no unseen source is claimed as inspected. This is a new packaged restriction of the inherited calculus, not a claim to have invented its moves or general exact inverse.

## 1. Physical construction and shifted product

For a,b >= 1, let T(a,b) have a root rho and a transitive spine. In spine order the root-edge labels are

    O O^a Z O^b Z.

The root beats an O vertex and loses to a Z vertex; earlier spine vertices beat later ones. Thus the graph is defined without a cosmetic orientation of arithmetic ties. Its order is a+b+4. It is strong: rho reaches the first spine vertex, every spine vertex reaches the last Z (or is that vertex), and the last Z reaches rho.

Put r=a+2 >= 3, s=b+1 >= 2. Then

    H(a,b)=1+2^(a+1)+2^(a+b+2)-2^b,
    2H(a,b)-3=(2^r-1)(2^s+1).                         (1)

Proof: after O^a Z the state has D=1+2^(a+1), U=2^(a+2)-1, d=1. The nonempty O^b block resets E to zero and multiplies U by 2^b; the last Z adds that U to D. This proves (1). The graph also makes the inherited state formula independently testable by direct subset DP.

More generally, for any physical seed (D,E,U,d), the same nonempty two bursts give

    2(H-D)-d=(2^(a+1)U-d)(2^(b+1)+1).                (2)

The first O erases E. The second factor is a generalized Fermat number (a genuine Fermat number only when b+1 is a power of two). Equation (2) is an exact divisor-based inverse criterion for a *fixed seed*: enumerate positive divisors B=2^(b+1)+1, b>=1, of the shifted target, and require ((2(H-D)-d)/B+d)/U to be a power 2^(a+1), a>=1. No claim of injectivity for arbitrary seeds is needed.

## 2. Root inverse: binary runs uniquely recover the whole word

The map (r,s) -> (2^r-1)(2^s+1) is injective on r>=3, s>=2. More precisely, write M=2H-3 in base two.

- r=s: M is a single run of 2r ones.
- r<s: M has exactly two runs, each of r ones, separated by s-r zeros. Its length is r+s and its leading run length is r>=3.
- r>s: M has leading run length one, total length r+s+1, and trailing run length s. Indeed its high bit is 2^(r+s), and its lower part is 2^r-2^s-1. That lower part has s trailing ones, then a zero, then r-s-1 ones (the last run may be empty).

These types are disjoint. Recover (r,s) as follows: for an all-ones string take half its length for both; otherwise, if its leading run has length >1, take r to be that length and s=length-r; otherwise take s to be the trailing-one run length and r=length-s-1. Check r>=3, s>=2 and equation (1). The final check rejects malformed strings. Thus the algorithm is an exact membership test and inverse for every positive integer H, with O(log H) bit inspection and a constant number of big-integer operations. It needs neither factoring nor a primality oracle.

The inequalities r>=3,s>=2 are essential to this clean parser: extending the domain can merge binary types, so no unrestricted exponential-product theorem is asserted.

## 3. Exact coverage size at every dyadic bound

For every integer k>=4,

    #{H(a,b)<=2^k : a,b>=1}
      = (k-4)(k-3)/2 + floor(k/2)-2.                 (3)

For k<=4 the image below 2^k is empty. In particular, for even k=2m>=4 the count is 2(m-1)(m-2). For arbitrary real X tending to infinity the count is

    (1/2)(log_2 X)^2 + O(log X).

Proof: injectivity permits counting ordered pairs. If r+s<=k, all pairs qualify; if r+s>=k+2, none do. This follows directly from M=2^(r+s)+2^r-2^s-1 and r>=3,s>=2 (the possible correction has magnitude strictly less than the leading power's half). At r+s=k+1, H<=2^k iff 2^r-2^s+2<=0, equivalently r<s. There are (k-4)(k-3)/2 pairs in the first triangle and floor(k/2)-2 pairs on the admitted boundary. Monotonicity between successive powers of two proves the asymptotic.

Thus the *exactly invertible arithmetic family* has zero density even among odd integers. This is a quantified global obstruction, not another modular obstruction. It places no bound on arbitrary O/Z words, arbitrary root changes, or all tournaments. The values 7 and 21 fail the inverse but their universal nonattainment uses the independent inherited theorem, not this restricted language.

## 4. Prime relations preserved by the shifted coordinate

Let p be an odd prime and t=ord_p(2). By (1),

    p | (2H-3)
      iff t | r, or (t is even and s = t/2 mod t).   (4)

Consequently the prime test is periodic separately in r and s with period t, which divides p-1; for odd composite moduli one can use the Carmichael exponent lambda(m), or phi(m), as a common (not necessarily least) period. This is an exact use of group exponents, not the distinct Carmichael totient-fiber analogy of THM-4904.

Prime powers are retained as well. If t|r, then the first factor has valuation

    v_p(2^r-1)=v_p(2^t-1)+v_p(r/t).

If t is even and s=(t/2)u with u odd, the second has valuation

    v_p(2^s+1)=v_p(2^(t/2)+1)+v_p(u).

Otherwise the corresponding valuation is zero. Add the two valuations for the shifted count. These are the elementary odd-prime lifting-the-exponent formulas; exceptional large base valuation (including Wieferich behavior) is retained, not silently set to one.

With g=gcd(r,s), shared prime factors satisfy the exact identity

    gcd(2^r-1,2^s+1) = 2^g+1 if r/g is even, and 1 otherwise. (5)

For completeness, write r=gR,s=gS with coprime R,S. Any common odd divisor makes x=2^g satisfy x^R=1,x^S=-1. If R is odd, Bezout or raising the relations to S and R gives -1=1, impossible. If R is even, S is odd and Bezout gives x=-1 modulo every common divisor, so it divides 2^g+1. The reverse divisibility is immediate. This proves equality, including prime-power multiplicities.

When r=s, M=2^(2r)-1 is a repunit, and H=2^(2r-1)+1 is divisible by 3. If r=s=2^j (j>=2), the classical identity makes M the product of Fermat numbers F_0 through F_j. The graph therefore captures an actual multiplicative natural-number relation in an explicit affine coordinate of its path count.

**Lost data:** (4) is not a primality test for H. H(1,1)=19 is prime with shifted product 35, whereas H(1,2)=33 is composite with shifted product 63. The preserved object is the product and its order/valuation data after H -> 2H-3. The inverse restores the two burst lengths only within this language; it does not identify arbitrary tournaments with that H.

## 5. Reproduction and proposed canon statement

Run `python3 05-knowledge/results/session241_prime_relations_20261110/coverage/verify.py`. It writes `verification.json`. Tests compare: recurrence versus closed form on 10,000 pairs, independent literal path DP on 16 graphs, exhaustive inverse rejection/acceptance on 500,000 odd targets, exact gcd identities, 72,000 prime-period/valuation predicates, 5,112 fixed-seed identities, and the exact dyadic cardinality formula at k=4,...,160. The checks use explicit raises rather than Python assertions, so optimized mode retains them.

Suggested canon statement, without reserving an ID: **The two-burst strong tournament family has an injective Mersenne–generalized-Fermat shifted product, a linear-bit membership inverse, and exactly (k-4)(k-3)/2+floor(k/2)-2 distinct Hamilton counts below 2^k.** Include the exact prime-order sidecar and distinguish this zero-density language from the still-open general Chan–Pak conjecture.

## Exact inherited canon paths

- `01-canon/theorems/THM-4888-finite-oz-banks-mod3-law-involution-rigidity-residue-confinement.md`
- `01-canon/theorems/THM-4916-root-change-local-factor-identities-confinement-multipliers-product-family.md`
- `01-canon/theorems/THM-4861-the-five-fermat-primes-are-hamiltonian-path-counts-of-tournaments-with-least-orders-3-4-6-8-11-a-vertex-transitive-tournament-has-a-prime-count-only-for-c3-and-fermat-primes-are-free-vertices-of-the-prime-square-sum-graph.md`
- `01-canon/theorems/THM-4986-the-7-21-omission-proof-has-an-exact-spanning-tree-twin-block-floor-crevasse-lemma-and-the-first-creator-contrast.md`
- `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`
- `01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md`
- `01-canon/theorems/THM-4904-root-change-word-fibers-carmichael-chamber-parity-law-sqrt-density.md`
