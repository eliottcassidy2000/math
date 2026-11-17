# A prime-capable exact collision bank and its quadratic inverse

**Status:** exact construction, exact inverse, and infinitely many distinct prime divisors proved. No assertion of infinitely many prime-valued terms, an arithmetic progression of all attainable integers, or Chan–Pak completeness.

## 1. A collision mechanism that really reaches a prime

The zero-shift construction in `joint_placement/arithmetic/PROOF.md`, with G=4 and t=24, gives address pairs

    (406426391,406403152), (406426407,406403136).

Both have lower count

    p=42284279505304223759,

which is prime. Each presentation also has its physical reflected count-p+2 companion. The same literal one-arc continuation produces respectively

    42284279712886532305, 42284279712886540517.

Their difference is exactly8212. Thus prime count does not force the scalar fiber to have one continuation state. In particular, this collision interface is not obtained merely by multiplying two nontrivial path counts through ordinal composition.

`PRIME_CERTIFICATE.json` and `verify.py` supply a standard-library Lucas certificate, not a probable-prime test. For each node n, its complete factorization of n-1 is certified recursively and a base a satisfies a^(n-1)=1 mod n and gcd(a^((n-1)/q)-1,n)=1 for every prime q dividing n-1. Thus a has order n-1 modulo every prime divisor of n, forcing n prime. Leaves at most100 are trial-divided. The certificate has eight nontrivial nodes.

This is a finite prime example, included to establish that the construction escapes a compositeness barrier, not as statistical evidence for coverage.

## 2. An exact infinite recurrence bank with no common prime obstruction

Specialize the earlier dyadic collision family to G=4 and t=2^(22+j), j>=0. Its shared source count is exactly

    Q(2^j),   Q(X)=a X^2+b X+c,
    a=1267674777822284490394719748096,
    b=135388433329826992581696,
    c=3614891648909077.

Every term retains two same-count physical presentations and the distinct continuation law, as well as the reflected +2. This is inherited from the suffix identity, not an interpolation of initial values. The characteristic roots4,2,1 give the exact recurrence stated in the earlier proof.

The discriminant b²-4ac is positive and not square; its integer square root is2824044124272513318. Thus Q is irreducible over the rationals. Also

    gcd(Q(1),Q(2),Q(4))=1.

So no prime divides every term of this exact infinite bank. These statements do not imply infinitely many prime-valued terms: even irreducible polynomials restricted to powers of2 can have difficult primality behavior.

### Infinitely many different primes nevertheless occur as factors

Suppose instead only finitely many primes S divide all the Q(2^j). All terms are odd. For each odd p in S put e_p=v_p(Q(1)) and choose j to be an arbitrarily large multiple of the multiplicative order of2 modulo p^(e_p+1), simultaneously for all p in S. Then

    Q(2^j)=Q(1) mod p^(e_p+1),

so v_p(Q(2^j))=e_p exactly. Under the finite-support assumption this forces Q(2^j)=Q(1) for arbitrarily large j. But Q has positive leading coefficient and tends to infinity, contradiction.

This argument is a useful finite-place/absolute-size separation: make all hypothesized prime valuations locally stationary while the real magnitude grows. It needs Q(1) nonzero and odd Q(2^j); no unproved distribution of primes is invoked. It applies to many other nonconstant positive polynomial dyadic banks, but not to pure monomial banks whose prime support can stay fixed.

## 3. An exact target decoder

Given N, form

    D=b²+4a(N-c).

N belongs to the shared-source dyadic bank if and only if all of the following hold:

1. D is a nonnegative perfect square;
2. sqrt(D)-b is divisible by2a;
3. X=(sqrt(D)-b)/(2a) is a positive power of2.

Then j=log2(X) and the addresses are reconstructed from t=2^(22+j). This is a direct exact inverse, not a scan over levels or binary reversal words. Integer square root and division give polynomial bit complexity. Subtracting2 first handles reflected companions. It is an exact membership criterion, not an existence claim for arbitrary N.

## 4. Constructive implication and unresolved issue

A bank can simultaneously carry exact count collisions, nontrivial continuation differences, no universal prime divisor, irreducible closed recurrence, and infinitely many distinct prime factors. These are compatible properties, rather than a forced choice between flexible interfaces and composite-only construction.

The remaining problem is still absolute placement. This particular bank has only O(log X) shared-source terms up to X. Its inverse sharply separates the arithmetic membership obstruction (a square discriminant followed by a power-of-two constraint) from the already solved physical realization. A future fusion must relax that restriction while preserving the continuation coordinate, not merely add more proofs of coprimality.

## Reproduction

Run `python bank_fusion/arithmetic/verify.py` from the workspace, or invoke the file by absolute path. It independently checks the prime certificate, exact address counts,150 large dyadic targets and their inverse, nearby reflected nonmembership, the nonsquare discriminant and gcd. It uses no third-party package. SymPy was used only to discover a factorization during certificate construction; replay does not require or trust it.
