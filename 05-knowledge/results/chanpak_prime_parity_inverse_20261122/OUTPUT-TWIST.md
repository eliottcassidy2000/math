# A genuine Liouville twist on the affine grammar, and the untouched twin correlation

**Status: PROVED in the scopes below.** The determinant Liouville sign
is controllable and can be held fixed in an unconditional joint
almost-prime theorem. The Liouville sign of the **output** is a
different observable; no cancellation for that sign or prime pair
infinitude is asserted.

## 1. Three coprime physical outputs from one address

Use the positive physical grammar of [the affine-sieve packet](../chanpak_affine_sieve_metagraph_20261122/AFFINE-SIEVE.md),
with chronological matrices

```text
A=[[3,3],[6,8]], B=[[9,9],[30,32]], C=[[3,3],[12,16]],
(h,R) starts at (1,2).
```

For each positive-depth word w, finish with the all-ones address
of width j, for j=1,2,3. The rooted migration formula gives the
three literal Hamiltonian counts

```text
N_j(w)=h+2^j R.       (j=1,2,3)
```

The finish adds a transitive block and its new root. The new root
beats that block, which beats the old tournament, and the old
terminal beats the new root. Every old vertex can reach the old
terminal along a suffix of the certified rooted full path.
Consequently each finish is strong. There are at most 6n+6 vertices
after n letters and any of these three finishes, and
3^n ≤ N_j ≤ 24·41^n.

All counts are odd. After the first letter the state has h=0 mod 3,
R nonzero mod 3, and these properties persist, so every N_j is
coprime to 6. For primes p≥5 the state is primitive, since the seed
is primitive and all three matrices have determinants supported on
2 and 3. The pairwise coefficient differences 2,4,6 are units
mod p. Thus no prime divides two different N_j. In each good
residue field the union of the three forbidden lines has precisely
3(p−1) nonzero vectors, yielding local sieve density 3/(p+1).
Any two output forms recover the state over their exact integer image.

There is no hidden multiple counting of triples at a fixed depth.
Indeed slope R/h is initially 2 and remains at least 2. According
as the last letter was A, B or C, its value belongs to the respective
disjoint intervals

```text
[22/9,8/3), [94/27,32/9), [44/9,16/3).
```

This follows by applying the corresponding fractional-linear
function to the interval [2,∞). The last letter is determined;
invert its rational matrix, requiring positive integral coordinates,
and repeat. The seed is out of all three image intervals, so the
algorithm stops exactly at the correct depth or rejects. From
(N_1,N_2), recover R=(N_2−N_1)/2 and h=2N_1−N_2. This is an exact
recognizer and inverse for the physical pair image, not a claim
that independently chosen output values (prime or composite)
satisfy the inverse.

## 2. Liouville of the determinant, not of the output

The determinants are 6,18,12. Thus complete multiplicativity gives

```text
λ(det M_w) = (−1)^(#B(w)+#C(w)).
```

The affine-sieve proof's two-step comparison extends to the
signed operator with signs +,−,− on A,B,C. More generally use any
fixed abelian letter character. On primitive residue states modulo
q coprime to 6, replace each unitary permutation ρ(U) by
χ(U)ρ(U). The same exact pairwise-square energy applies to
the average of its nine two-letter operators. The relative moves
G=BA^(-1) and D=CA^(-1) acquire phases, but their commutator
J=DGD^(-1)G^(-1), and independently the commutator of the
conjugated relative moves giving AJA^(-1), acquire phase one.
Hence the two integral shears used in the existing
[Bourgain–Varjú comparison](https://arxiv.org/abs/1006.3365)
are still obtained as bounded-length unitaries with no phase.
The same mean-zero operator norm gap is uniform in both q
and the character. For the determinant twist, on the constant
vector the one-step eigenvalue is (1−1−1)/3=−1/3. Thus for
either determinant sign ε, any measurable U in the finite
primitive state set Ω_q has

```text
Pr[state_n in U and λ(det M_w)=ε]
    = |U|/(2|Ω_q|) + O(sqrt(|U|) exp(−η n)),
```

with η>0 independent of q. The signs are exact arithmetic signs
of an integer determinant, not assumptions about independence of
Liouville values.

Using squarefree d prime to 6 for divisibility of N_1 N_2 N_3,
the same CRT and level-of-distribution argument as in the
affine-sieve packet now has local density

```text
β(d)=product_(p|d) 3/(p+1)
```

and yields a power-saving summed remainder up to D=exp(δ n),
for some δ>0. The [combinatorial lower-bound sieve](https://www.ias.edu/sites/default/files/math/Bourgain_Gamburd_Sarnak_Affine.pdf),
with this fixed dimension three and z=D^(1/28) (Theorem 3.1
uses s>9t for dimension t), therefore gives ≫3^n/n³ words in each determinant
sign for which **all three** coprime counts have no prime divisor
at most z. All three counts then have an absolute bound on Ω,
counting multiplicities, since N_j≤24·41^n. These are
≫3^n/n³ distinct ordered triples at depth n by the exact inverse
in §1. Moreover at least two of the three finish types give
≫3^(n/2)/n^(3/2) distinct **scalar** counts at this depth: any pair
of coordinate projections identifies the ordered triple, so the
product of the two smallest projection cardinalities is at least
the number of triples. The finishes are not claimed to be close
together, or prime.

The same proof applies in any fixed compatible primitive
residue phase: mix modulo qd and retain the fixed state modulo q.
For example, at each of finitely many primes p≥5 one can
prescribe a positive exact valuation e to a selected finish j
and keep the other two coprime to p: take
R=1 and h=−2^j+p^e modulo p^(e+1). The same Ω bound then holds
for the remaining parts, after removing the prescribed factors.

The number three is convenient because it needs no further
coprimality phase. More generally fix any J≥2 and use j=1,…,J.
Place in a fixed bad modulus all primes p≥5 dividing some
coefficient difference 2^i−2^j and all such primes p≤J.
At that modulus retain state (1,0). All the finished counts
are then units at every bad prime. At any remaining prime
the J lines are distinct and the forbidden proportion is
J/(p+1)<1; the dimension-J sieve applies. Thus both determinant
sectors contain ≫_J 3^n/n^J distinct ordered J-tuples of
pairwise coprime strong counts with Ω bounded in terms of J.
This observation changes the number of arithmetic slots, not
their gaps.

## 3. The output twist is not the determinant/letter twist

Two literal chronological words, CAABBBB and CBBBBAA, have the
same exact letter histogram (2,4,1), the same determinant, and
the same state (110,87) modulo 385. Their first finishes are

```text
28274210479 = 19 · 1488116341  (both factors prime),  λ=+1;
27107951539 is prime,                            λ=−1.
```

Both outputs are units modulo 385. Thus even exact letter counts
and that joint congruence state cannot recover the output sign.
Already AAB and BAA, which have the same determinant and the same
state (3,8) modulo 11, have outputs 23207=23·1009 and 22811
prime. These finite witnesses disprove an exact reduction to
the determinant character or these finite observers. They are
not statements about cancellation of the true output twist.

## 4. Exact signed barrier on a physical +2 chart

Now use an actual constant-gap branch from the
[unit-switch packet](../chanpak_affine_sieve_metagraph_20261122/UNIT-SWITCHES.md):
its two outputs are N and N+2. This is a different question from
the three simultaneous almost-prime finishes above. Let N+2≤Y,
let z³>Y, and let Q_z select words for which both outputs
have no prime divisor at most z; restrict to outputs >1.
Every z-rough output in this range has at most two prime factors
with multiplicity. On this set the prime indicator is exactly
(1−λ)/2. Consequently for arbitrary nonnegative word weights a(w),

```text
4 T_twins = U − S_0 − S_2 + S_02,
U    = Σ_(Q_z) a(w),
S_i  = Σ_(Q_z) a(w) λ(N+i)   for i=0,2,
S_02 = Σ_(Q_z) a(w) λ(N)λ(N+2),
```

where T_twins is the corresponding mass of actual twin outputs.
This is an exact identity including multiplicities; no error terms
are hidden. An alternative puts Λ(N) in the first slot and only
requires N+2 to be z-rough:

```text
2 Σ_rough a(w) Λ(N) 1_prime(N+2)
  = Σ_rough a(w) Λ(N) − Σ_rough a(w) Λ(N)λ(N+2).
```

Replacing Λ by θ(N)=log(N)·1_prime(N) gives the exact twin
identity with the same algebra and no prime-power ambiguity.
For Λ itself one must account for prime powers: for the explicit
bilinear chart N=Auv+Bu+Cv+d, the number of histories for one N
is at most the divisor-pair count of A(N−d)+BC by the exact inverse.
Thus with unweighted histories the prime-power contribution is
O_ε(Y^(1/2+ε) log Y), an ambient estimate that may be much too large
for a sparse word atlas. It is not automatically negligible.

If z is smaller, parity alone does not isolate primality. The
exact correction to (1−λ(n))/2 on a rough n>1 is the indicator
of Ω(n)=3,5,7,… up to the available finite factor bound; those
layers must be controlled in addition. In particular the weak
sieve scale in §2 cannot be inserted in these identities to get
primes for free. A real-positive unsigned lower bound and the
specific output-signed estimates (or another exact selector) are
still needed on the genuine word measure.

## 5. Why the conditional 180 result does not bridge that step

[Allahyar–Yianni, §§1–2](https://github.com/sahibzada-allahyar/math/blob/main/preprint/prime-gaps-180/prime_gaps_at_most_180.pdf)
uses a 39-shift admissible tuple of diameter 182 and a nonnegative
weight. There is a single pair of shifts at distance greater
than 180, the two endpoints. Its detection criterion subtracts
the expected mass of this bad double hit as well as baseline mass
from the sum of the single-prime masses. The preprint explicitly
depends on A4, an unestablished Type III exponential-sum input.
The argument and its conditional conclusion do not give a
Liouville or von Mangoldt estimate on the physical word measures
above, or select gap exactly two.

There is a sharp finite information statement for that detection
interface. Write M for total mass, b for bad-pair mass, and a_i
for singleton prime-slot marginals. If the only bad pair is
{0,last}, then an abstract nonnegative measure with those
statistics and **no good pair** exists exactly when the natural
constraints hold and Σ_i a_i≤M+b. Indeed put mass b on that
bad pair, mass a_i on each interior singleton, mass a_i−b on
the two endpoint singletons, and put the remaining mass on the
empty set. The necessity is the pointwise bound |S|≤1+1_bad(S)
when there is no good pair. Thus even perfect knowledge of those
statistics gives no better universal detection criterion; a
different choice of physical slots or genuinely new mixed
moments is needed to improve it. This is a deduction about the
interface, not a re-check of the prime-gap proof.

## Reproduction

The companion checker verifies the explicit matrices, the exact
small-depth pair inverse, determinant parity, joint coprimality,
and the literal output-sign counterexamples by deterministic
trial division. Normal and optimized Python give identical results
on all 9,840 nonempty positive words through depth eight.

Exact finite replay: `python3 check_output_twist.py`.
