# A constructive bridge from perfect-number deficiency to amicable pairs

Date: 2026-11-04. Status: rigorous classical reduction, explicit modular corollary, reproducible finite certificates. **No new infinitude theorem or claimed resolution of an open problem.**

## The adjacent yes/no question

Are there infinitely many amicable pairs? This remains open. An amicable pair is a pair of distinct positive integers N,M with sigma(N)=sigma(M)=N+M. Thus perfect numbers are fixed points of s(n)=sigma(n)-n, whereas amicable numbers are genuine 2-cycles. This is a closer constructive bridge to the perfect-number equation than an unsupported analogy with twin primes.

The proposal is to turn the shared divisor-sum structure into exact, finite certificates, then isolate precisely what would still be required for infinitude. The central factorization below is classical Euler, not a novel theorem: Euler's Problem 1, §27, uses x=p+1, y=q+1 and the reduced fraction a/(2a-sigma(a))=b/c to obtain (cx-b)(cy-b)=b². Our unnormalized version exposes the deficiency and supports a directly reproducible checker.

## Fixed-kernel certificate theorem

Let a be a positive integer, S=sigma(a), D=2a-S, and E=S-a. Consider distinct primes p<q and a prime r, with gcd(a,pqr)=1. Then

    N = apq,   M = ar

are amicable if and only if D>0 and there are positive divisors d<e with

    de=a²,
    d ≡ e ≡ -E (mod D),
    p=(d+E)/D,  q=(e+E)/D,
    r=pq+p+q,

where the primality, distinctness, and coprimality hypotheses above hold.

### Proof

Multiplicativity gives sigma(N)=S(p+1)(q+1) and sigma(M)=S(r+1). Equality of these sums forces r=pq+p+q. Their common value equals N+M precisely when

    S(p+1)(q+1) = a(2pq+p+q),
    Dpq - E(p+q) - S = 0.

Here E≥0 and S>0, so D must be positive. Completing the product gives

    (Dp-E)(Dq-E) = E²+DS = a².

Both factors are positive: the preceding equation gives q(Dp-E)=Ep+S>0 and symmetrically p(Dq-E)>0. Setting d=Dp-E and e=Dq-E proves necessity, including d<e. Conversely, substitute the asserted divisors, expand the product, and reverse the equalities. Finally r>pq, hence N<M and the cycle is nontrivial. QED.

This theorem excludes every perfect or abundant kernel a from this particular ansatz. It does not exclude perfect numbers themselves, nor amicable pairs of other forms.

## The 12±1 connection: an exact one-sided sieve

If p,q>3 and r=pq+p+q is prime, then r≡11 mod12. Indeed r+1=(p+1)(q+1) is divisible by four. If both p,q were 1 mod3, r would be 0 mod3; since r>3 this would contradict primality. Hence at least one is -1 mod3, so r+1 is divisible by three. Therefore 12 divides r+1.

Thus the single-prime cofactor must be 12k-1, never 12k+1. The two-prime side may use any prime residue classes 1,5,7,11 mod12 except that both factors cannot lie in {1,7}. This is a necessary local condition, **not** a primality guarantee or an infinitude argument.

## Worked certificates

- a=4: S=7, D=1, E=3; choose d=2,e=8. Then p=5,q=11,r=71; N=220,M=284 and both divisor sums are 504.
- a=585=3²·5·13: S=1092, D=78, E=507; choose d=351,e=975. Then p=11,q=19,r=239=12·20-1. N=122265 and M=139815; both divisor sums equal 262080. This explicitly demonstrates that odd–odd amicable pairs exist; that existence question is not open.
- a=819: S=1456,D=182,E=637; d=273,e=2457 gives p=5,q=17,r=107 and the odd pair (69615,87633).

## Reproducible validation and scope

Run `python verify_amicable.py > results.json` with Python's standard library only. It enumerates every positive kernel a≤5000 and every factor pair of a², checking the exact certificate by factorization and verifying both divisor sums again using a separate divisor-pair summation engine. Output: 10 certificates, of which 4 have odd a. It also compares against 9,975 separately enumerated small prime-pair cases for a≤100,p,q≤100, and exhaustively checks all 16 reduced residue pairs mod12.

This is exhaustive for the specified ansatz and kernel bound, not for all amicable pairs up to a numerical size. It uses exact trial-division primality tests, not probabilistic tests. The attached JSON contains full witnesses. Explicit checks remain active under Python -O; optimized and normal runs give identical output.

## The remaining global obstruction

Each fixed a offers only finitely many divisor pairs. Therefore infinitude through this route would require infinitely many different deficient kernels whose divisor pairs produce three simultaneous primes while satisfying coprimality. A congruence sieve can reject candidates but cannot supply that global theorem. No implication from odd-perfect-number nonexistence to amicable infinitude, or conversely, has been established here.

## Primary sources and source audit

1. Euler, *De numeris amicabilibus* (1750), §27; Jonathan David Evans translation, revised September 2025: https://arxiv.org/html/2409.08783v2 . This identifies the exact classical precedent for the factorization.
2. Carl Pomerance, *On amicable numbers*, author manuscript dated March 30, 2014: https://math.dartmouth.edu/~carlp/amicable032914.pdf . Establishes a sparsity bound and explicitly distinguishes the unresolved infinitude question.
3. Paul Pollack, *Amicable numbers*, submitted October 5, 2026: https://arxiv.org/abs/2610.05661 . Reports a substantially stronger asymptotic upper bound for amicable-number density and explicitly states in its introduction that infinitude remains unproved. An upper bound does not prove infinitude.
4. Source-quality caution: Chen–Tang–Zhan, arXiv:2601.07444v1, introduction and §3.8, https://arxiv.org/html/2601.07444v1 , incorrectly treats odd–odd amicable existence as open / proposes their nonexistence. The certificates above refute that assertion. Its other formalization claims were not independently audited and are not needed for this note.
