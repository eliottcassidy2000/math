# Independent arithmetic and code audit

2026-10-29. Reviewed `arithmetic/ARITHMETIC.md` and `codes/SYMMETRY.md`. **PASS within the stated scope.** No blocking mathematical defect found. The audit script is standalone stdlib Python and imports no producer code; normal and optimized execution both pass. The JSON records are reproducible with `python3 check.py` and `python3 -O check.py`.

## Citation and theorem boundary

Independently retrieved [Bennett–Levin, *The Nagell–Ljunggren equation via Runge's method*](https://personal.math.ubc.ca/~bennett/BeLe-Mona-2015.pdf), Monatsh. Math. 177 (2015), 15–31, DOI 10.1007/s00605-015-0748-1. Proposition 1 is on printed page 16 (PDF page 2). Its first condition excludes square exponents outside the displayed known solutions. For positive base, the two square cases give sigma(3^4)=121 and sigma(7^3)=400. Exponent one must be handled separately and gives sigma(3)=4. This supports the complete prime-power classification, but not a general perfect-power classification. The paper aggregates classical results for this proposition; this session does not establish that deep theorem anew. The script's bounded search is only a guard.

## Independently confirmed

- Divisor enumeration, not producer factorization routines, reconstructs the bank's seven distinct nonzero F2^3 columns.
- Kernel/image dimensions and distances are [7,4,3]/[7,3,4]. The zero column gives [8,5,1]/[8,3,4], including the stated kernel weight distribution.
- Exhaustive coordinate permutations yield group orders 168 for both seven-coordinate and zero-extended response codes; adding the constant row changes the group to order 1344. The affine code is [8,4,4], self-dual, and punctures to the original kernel.
- The 19 stabilizer has 24 elements; its action on four avoiding lines has 24 distinct permutations, proving the claimed faithful S4 identification. There are seven lines, 21 flags, and one flag orbit.
- Joint direct/sigma columns for 3,5,19 have rank four.
- sigma(5*19*29)=3600. Six zero-state products preserve sigma squareclass when the residual factor is coprime to 21. Exponent replacement and telescoping check independently; the shared-port cocycle checks on all triples from 1 through 15.
- Coprime multiplicativity and the distinct aliquot identity check on coprime pairs below 60. sigma(9)=13, sigma(243)=364, s(15)=9 and s(81)=40 independently confirm the hostile controls.

## Scope and caveats

The algebraic proofs in the producer notes justify universal claims; finite checks do not replace those proofs. Exponent-state symmetries preserve squareclass, not numerical divisor sum. Neutral attachment requires a coprime base. Fixed prime labels are not freely interchangeable arithmetic identities. The Fano 7/21 count does not identify any unrelated 21 objects. No novelty or convergence inference is warranted. There are no unresolved blockers in the audited arithmetic/code scope.

## Addendum: arithmetic parity completion (2026-10-30)

Reviewed new codes §6 and independently extended `check.py`. The joint observer (prime-count parity, sigma squareclass) has sixteen fibers of sixteen selections; its zero fiber equals the self-dual affine code already checked. All 256 selections are checked using explicitly constructed divisor lists, not producer routines. All 344,064 affine-map/selection identities independently confirm transport `(epsilon,s) -> (epsilon,Ls+epsilon*t)`. Translation takes the singleton prime3 to singleton prime7 in the displayed hostile case, so it does not preserve the old sigma-zero fiber. The larger zero-fiber group and changed observer are correctly distinguished. A dual-character index is a choice of character, not a fresh random bit from the observation prime. No defect found in the completion.

Also checked the arithmetic common-unit caveat independently: sigma(3)=4 has Legendre value +1 modulo11, whereas sigma(81)=121 has value zero. Equal rational squareclass does not preserve raw modular zero/support observations without excluding primes dividing representative values. On the six-state family, excluding prime divisors of sigma(u)·2·5·11 suffices.

## Final focused addendum: Phi6 selection and typed return (2026-10-30)

**PASS.** Independently reviewed final arithmetic §6. Given the cited complete neutral prime-power classification, solving `x²−x+1=n` for n=3,81,343 reduces to discriminants 9,321,1369 and positive integer roots x=2,19. Checking the return relation removes x=2 and leaves x=19. The universal identity `Phi6(x)+3x=(x+1)²` establishes the claimed equivalence to the particular square value, rather than merely an unspecified square. The prime-power guard remains explicit.

Independently checked the prime-cube proof by formal coefficient convolution, not point sampling: both sides of the stated identity have ascending coefficient vector `[7,-1,0,-7,1]`. Thus p>1 forces p=7, and x=19; using s(p³)=p²+p+1 requires and has the stated primality condition. This elementary restricted proof does not purport to replace the cited classification on all prime powers.

The signed polynomial evaluations at 18 and −19 both give343, while independent divisor enumeration gives sigma(324)=847. The typed sequence uses polynomial evaluation, aliquot sum, and division by3, respectively. Direct guards s(19)=1 and s(57)=23 further exclude reading it as an ordinary aliquot cycle. No Collatz conclusion follows or is claimed.

`check.py`, `normal.json`, and `optimized.json` now include this final focused addition. No further extension or publication was performed by this auditor.
