# The prime 463 is a generator of the entire dyadic conductor clock

2026-10-17, clock_duality_oct15. **PROVED, scoped and self-audited:** the
valuation law, exact conductor-kernel generator, and all-level group structure
below. **FINITE-EXACT:** 256 power checks, independent quadratic-form
composition through conductor256, and the separate discriminant-463 control.
No independent referee or novelty claim. No Collatz convergence, plane-colouring,
or planar Jacobian conclusion.

## Inheritance: a coefficient becomes a marked ideal

The preceding [conductor cutoff D2](escape_renewal_idoneal_20261017d.md)
proved the class numbers and fixed-class counts, but gave the full group
structures only as a finite PARI table. Its first non-self-inverse witness
was (16,8,463). The present calculation upgrades that table to all conductors.
The least-used sidecar is the ideal over463 itself, rather than the scalar
class number or the inverse-merged orbit count.

Other recovered occurrences of the numeral are typed separately below.
In particular THM-463 is a theorem identifier for the two-cube divisor
criterion, not a theorem about the prime463.

## P1. The exact clock

Let m be positive, odd, and squarefree, put theta^2=-2m, and let

    O = Z[theta],  O_k = Z + 2^k O,  alpha=1+theta.

O is the maximal order of Q(theta), with discriminant -8m and only the
units +1,-1. Write alpha^n=A_n+B_n theta for n>=1. Then

\[
                 A_n\text{ is odd},\qquad v_2(B_n)=v_2(n).       \tag{1}
\]

**Proof.** Modulo2, theta^2=0, so (1+theta)^n=1+(n mod2)theta.
In particular A_n is always odd and B_n is odd for odd n. Squaring gives

\[
 A_{2n}=A_n^2-2mB_n^2,\qquad B_{2n}=2A_nB_n.                  \tag{2}
\]

Each squaring therefore raises v2(B) by exactly one. Strip the power of2
from n to obtain (1). There is no probabilistic or asymptotic step.

## P2. A cyclic conductor kernel at every level

Let I_k be the contraction of the O-ideal (alpha) to O_k. Its norm 2m+1 is
odd, so it is coprime to the conductor and is an invertible proper ideal.
Its class lies in the kernel of Cl(O_k)->Cl(O), and

\[
                     \operatorname{ord}[I_k]=2^k.              \tag{3}
\]

Indeed I_k^n is principal exactly when alpha^n lies in O_k. Necessity:
extending a principal generator to O gives the same ideal as alpha^n;
the generators differ only by an O-unit, +1 or -1. Sufficiency uses the
extension/contraction bijection for conductor-coprime ideals and the fact
that alpha^n has odd norm. By (1), membership in O_k is equivalent to
2^k dividing n.

The classical order class-number formula gives h(O_k)=2^k h(O): the
conductor prime2 is ramified and the unit index is one. Thus the entire
kernel, not just a subgroup found in a finite search, is cyclic:

\[
       \ker(\mathrm{Cl}(O_k)\longrightarrow\mathrm{Cl}(O))
                         =\langle[I_k]\rangle\cong C_{2^k}.    \tag{4}
\]

The extension/contraction input and order class-number formula are classical;
[Chao Li's exposition, Ring class fields, Proposition2 and Theorem9](https://www.math.columbia.edu/~chaoli/docs/MinorThesis3.html)
records the precise conductor-coprime restriction. The explicit generator
and valuation calculation above are the mechanism used here.

If Cl(O) is elementary abelian of rank r=omega(m), the genus count from D2
is |Cl(O_k)[2]|=2^(r+1) for k>=1. Its order is2^(k+r), and (3) gives an
element of order2^k. The structure theorem for finite abelian groups now
forces

\[
                  \mathrm{Cl}(O_k)\cong C_{2^k}\times C_2^r,
                  \qquad k\ge1.                              \tag{5}
\]

To see that no hidden larger invariant factor remains: r+1 cyclic factors
consume at least r bits outside the largest one. That largest factor has
at least k bits by (3), exhausting the total k+r bits. This conclusion
requires the stated maximal-order hypothesis; it is not a claim about all
quadratic class groups.

## P3. Why463 and1848 meet

Take m=231, so theta^2=-462 and

\[
                        N(1+\theta)=463.                     \tag{6}
\]

463 is prime. The complete reduced-form list for discriminant-1848 has
eight ambiguous forms, proving Cl(O)=C2^3. Therefore

\[
 \boxed{\mathrm{Cl}(-1848\,4^k)=C_{2^k}\times C_2^3\quad(k\ge1).} \tag{7}
\]

At k=1 the idoneal parameter is1848 and the group is C2^4. At k=2 it is
C4 times C2^3. A prime ideal over463 represents (16,8,463) or its inverse,
according to the chosen prime above463; these are exactly the first
non-self-inverse witness classes from the preceding packet.

The same marked prime generates the full conductor kernel at every later
level. The first squaring steps are

    (1,1) -> (-461,2) -> (210673,-1844),

where the pairs are (A,B). The imaginary coefficients have valuations0,1,2.
This is why merely retaining the inverse orbit would lose the growing clock.

The group counts require no further form enumeration:

\[
 \begin{aligned}
 |\mathrm{Cl}(O_k)|&=8\,2^k,\\
 |\operatorname{Fix}(g\mapsto g^{-1})|&=16,\\
 \#(\text{inverse-merged classes})&=4\,2^k+8 \qquad(k\ge1).
 \end{aligned}                                                \tag{8}
\]

Multiplication by the marked class above463 gives eight directed cycles,
each of length2^k, one in each conductor-kernel coset. This is an intrinsic
binary relation, but not a tournament: most vertex pairs are unrelated.
Inversion reverses these arrows. The all-order converse/fusion theorem for
tournaments supplies the same involution orbit formula, not an identification
of the two source categories. Retaining the marked ideal and its power is the
analogue of retaining a Hamiltonian observer path above a merged metagraph
node; a quotient count alone does not reconstruct either source.

## P4. Other appearances, without conflating their predicates

The earlier [Paley/Heawood note, section6.4](procgen_petersen_20261001_petersen_heawood_paley.md)
lists463 among primes whose fourteen-root induced Paley tournament is
all-odd. Independently, reduced forms give h(-463)=7. The exact finite
character checks are

    sum_(a=1..462) a (a/463) = -3241 = -463*7,
    sum_(a=1..231) (a/463) = 7.

These retain the additive ordering of the quadratic-residue signs, while the
fourteen-root tournament restricts those signs to multiplicative positions.
No implication between its all-odd property and h(-463)=7 is asserted here.
The fields Q(sqrt(-463)) and Q(sqrt(-462)) are different. In the latter463 is
a split norm prime; in the former it is the ramified discriminant prime.

The reusable conclusion is a representation with a proved update law:
(A,B) remembers the conductor depth that a norm or inverse orbit forgets.
For Collatz a valuation law at a fixed centre still does not control changing
centres or replenished excursions; the earlier finite-centre obstructions
remain in force. There is no transfer of (1) to universal integer descent.

## Reproduction

[Checker](../../04-computation/experiments/prime463_conductor_20261017e.py),
[receipt](../../04-computation/experiments/marked_contact_current_20261017e/prime463.json).
It uses exact integer multiplication for n=1..256 and independent PARI
binary-form composition, testing both the proposed order and its half for
k=1..8. Reduced forms independently give the maximal group and h(-463).
Normal and optimized outputs agree byte-for-byte. The all-order assertions
come from (1)--(5), not from extrapolating these controls.
