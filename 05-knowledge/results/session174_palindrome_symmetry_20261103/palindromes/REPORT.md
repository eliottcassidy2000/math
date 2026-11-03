# Palindromic quadratic norm sieves and tournament reversal

**Status: PROVED** elementary statements below, supplemented by **FINITE-EXACT** tests. No global novelty or prime-infinitude claim. Restored after /tmp reset, including independent audit scope corrections.

## 1. The record-form quadratic field

For b>=2,m>=1, put P_m(b)=b^(2m+1)-b^m-1. Its base-b digits are m copies of b-1, then b-2, then m copies of b-1. The exact identity is

(2b^(m+1)-1)^2=4b P_m(b)+(4b+1).

Every odd prime divisor q of P_m(b) is coprime to b and has 4b+1 square or zero modulo q. Proof: P=-1 modulo prime divisors of b, and the displayed identity gives the square. Away from the discriminant, such primes split in Q(sqrt(4b+1)); when the discriminant is square, this field language degenerates. This is only a necessary factor restriction.

**Decimal strengthening.** For m>=1 every prime divisor q of P_m(10) is a nonzero quadratic residue modulo 41. The final digit 9 excludes 2,5. Powers of 10 modulo 41 are 1,10,18,16,37; the only root of10x²-x-1 modulo41 is39, so 41 never divides P. Quadratic reciprocity gives (q/41)=(41/q)=1. Hence every positive divisor also belongs to the square subgroup modulo41. This does not mean q=1 modulo41: 989=23*43 is a hostile control.

## 2. Exact exponent automaton and paired residues

For odd prime q not dividing b, let t=ord_q(b), and
R_q={r modulo t: b^(2r+1)-b^r-1=0 modulo q}.
Then q|P_m iff m mod t belongs to R_q, and |R_q|<=2. The state x=b^m advances by multiplication by b; accepting states solve bx²-x-1=0.

If t is odd, |R_q|<=1. If t is even, R_q is invariant under
r -> t/2-1-r modulo t.

Proof: two roots have product -b^(-1). For odd t, -1 is outside <b>, so both roots cannot be inside it. For even t, -1=b^(t/2), yielding the stated partner, including a double root. Discriminant squarehood alone loses the subgroup/phase information and is not sufficient.

Decimal examples (q,t,R): (23,22,{1,9}), (31,15,{12}), (43,21,{1}), (59,58,{9,19}), (61,60,{14,15}), (83,41,{32}). The user's prime 7417 supplies a particularly concrete bridge: ord_7417(10)=2472 and its exact accepting residues are 1312 and 2395. Their sum is 1235 modulo 2472, exactly the reflection constant t/2-1. Thus 7417 is an actual divisor sieve for this palindrome family. A finite-prime sieve is periodic with period dividing the lcm of the orders. Shared order factors prevent naive independent-density multiplication. Divisibility proves compositeness only when P>q.

## 3. Belphegor and symmetric defects

Let B_n be the decimal word 1, n zeros, 666, n zeros,1 (n>=0). For x=10^(n+1),
B_n=100x²+666x+1,
(100x+333)²-110789=100B_n,
110789=7³*17*19=49*2261.
Thus any prime divisor q outside {2,5,7,17,19} satisfies (2261/q)=1. Modulo each q not dividing 100, at most two states of <10> can be roots. Belphegor's 13-zero example is n=13; its primality is not tested or needed here.

For b>=2,m>=d>=1 define
N_(m,d)(b)=b^(2m+1)-b^(m+d)-b^(m-d)-1.
With x=b^(m-d), A=b^(2d+1), C=b^(2d)+1,
N=Ax²-Cx-1,
(2Ax-C)²=4AN+C²+4A.
For fixed d this gives the same discriminant and at-most-two-root mechanism.

**Canonical palindrome boundary:** N has a palindromic base-b expansion if b>=3 OR m>d. At b=2,m=d, both outer digits become0; fixed-width palindromy remains but ordinary canonical palindromy need not. Explicit hostile: N_(1,1)(2)=2 has digits10, not palindrome, although padded010 is symmetric. The algebraic identity remains valid. Also d=0 subtracts 2b^m and is not the single-defect family.

The reported 2,000,005-digit two-defect example has m=1,000,002,d=51,044. This is a symbolic interpretation, not independent prime certification.

## 4. Cyclotomic obstruction and actual tournament bridge

Repeating a digit block W of length ell, r times gives
W(b)(b^(ell*r)-1)/(b^ell-1).
Examples: 1001=Phi_2(10)Phi_6(10)=11*7*13; 2222=2Phi_2(10)Phi_4(10)=2*11*101.
For F(z)=z^d F(1/z) with odd formal degree d, F(-1)=0. Thus all even-digit base-b palindromes are divisible by b+1.

For a tournament T on n vertices, let f_T(p) count forward adjacent arcs along permutation p. Set
F_T(z)=sum_p z^(f_T(p)), S_T(z)=sum_p sgn(p)z^(f_T(p)).
Reversal sends f to n-1-f and permutation sign to epsilon times itself, where epsilon=(-1)^(n(n-1)/2). Hence
F_T(z)=z^(n-1)F_T(1/z),
S_T(z)=epsilon z^(n-1)S_T(1/z).

For even n, (z+1)|F_T(z), so b+1 divides F_T(b).
**All signed-factor claims below assume n>=2.** S_T(1)=0 because permutation signs sum to zero. When epsilon=1, differentiation gives 2S_T'(1)=(n-1)S_T(1)=0, hence a double zero at 1. Substitution at -1 supplies the other forced factor. Guaranteed signed factors by n modulo4:
0: (z-1)²(z+1); 1: (z-1)²; 2: z-1; 3: (z-1)(z+1).
They hold also when S_T=0. At n=1, S_T=1, so the restriction is essential.

Source-target map: formal coefficient reversal to polynomial evaluation at integer b. Preserved: divisibility forced by evaluation at -1. Lost: coefficients can exceed b-1 and carries can destroy literal digit palindromy. Required sidecar: the polynomial and its formal reversal degree. No artificial tournament on integers is constructed.

Repository connection: THM-583-palindromic-count-is-half-system-hampath-count.md reconstructs a path's second half using an involutory anti-automorphism. Here reversal reconstructs a numeral's second half. Retaining the involution makes both compressions lossless. This is not a reduction of prime existence to Hamiltonian paths.

## 5. Signed reciprocity correction

THM-088-signed-F-palindromic.md(C) correctly derives, for Q=S_T/(z-1),
Q(z)=(-1)^(binom(n,2)+1)z^(n-2)Q(1/z),
but then drops the parity-dependent sign and claims Q is always anti-palindromic. That conclusion fails: n=2 gives S_T=±(z-1), so Q is a nonzero constant. The corrected formula follows directly from reversal. No prohibited mistakes ledger was read and no older theorem modified by this worker.

## 6. External scope and research questions

Primary Banks--Hart--Sakata source: https://arxiv.org/abs/math/0405056, Almost all palindromes are composite. Its all-palindrome asymptotic neither settles infinite palindromic primes nor transfers automatically to these exponentially thin families.

A discovery table https://www.stdkmd.net/nrr/prime/prime_n2pl.htm reports the 2,000,005-digit Propper/Batalov example dated 2024-01-06 (t5k 136852). This is not an audited certificate. The 2021 1,888,529-digit number must not be called current largest without primary-record verification; parent source lane owns this issue.

OPEN directions: design central blocks with favorable discriminants but no fixed divisor; solve finite covering systems of exponent progressions with CRT dependence retained; compare one-/two-defect families through fields and orbit hits; use reversal factor signatures to check tournament enumerators. None is an assertion of infinitely many primes or an LRC result.

## Verification

Run python3 verify.py and python3 -O verify.py from this directory. Checks use explicit exceptions, never optimization-disabled assert statements. The test universe covers bases 2..16, m=1..15; bases 2..12, primes<=300 for root pairing; decimal Belphegor n=0..29; symmetric defects bases 2..12, m<=9 with the binary boundary; every labelled tournament n=2..5 and every permutation, plus explicit n=1 hostile. Finite tests supplement proofs, not replace them.
