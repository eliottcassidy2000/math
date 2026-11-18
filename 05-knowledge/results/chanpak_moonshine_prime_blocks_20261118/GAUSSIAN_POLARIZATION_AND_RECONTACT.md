# Gaussian polarization and the physical recontact determinant

**Status:** PROVED structural synthesis and scoped corollaries; literal examples independently replayed. The diagonal Gaussian reader is inherited, not new. No universal prime placement, Chan–Pak proof, or infinitude of prime outputs is claimed.

## Provenance and the question being retained

`pole_prime_placement/NORM_READERS.md` already gives a literal Gaussian reader H²+P₂², alongside readers of discriminants−8,+8,−44,+88 and a bank with no common quadratic splitting obstruction. It explicitly separates arithmetic representation from physical profile existence. `machin_tournament/README.md` uses Gaussian multiplication for exact angular corrections, while `minus_machin/README.md` records that cancellation does not remove all physical source factors.

The present synthesis connects that inherited Gaussian reader with the **signed determinant from selective recontact** in `niche_bridge/ORDINAL_RECONTACT_DETERMINANT.md`, and with the shell identity from `square_reflection_bridge/REFLECTION_SHELL_LAWS.md`. The new useful perspective is a single complex coordinate carrying two different physical reads, not a claim that Gaussian norms alone classify tournament counts.

## 1. Physical polarization of the norm reader

For a nonempty tournament A let

h_A=H(A), p_A=P₂(A), z_A=h_A+i p_A.

All path counts are directed. P₂ counts ordered two-path covers, allowing empty slots. Given two nonempty tournaments A,B, form their ordinal sum A→B, then its cyclic shell: add a,b with b→a, a→A,B, and A,B→b. The resulting tournament is strong. Ordinal composition multiplies slot counts and a nonempty shell has H'=H+P₂. Hence its **actual Hamilton count** is

N(A,B)=h_A h_B+p_A p_B.                              (1)

When A=B this is the inherited Gaussian norm |z_A|². For different sources it is a dot product, **not generally a Gaussian norm of one integer pair**.

Define the signed determinant

D(A,B)=h_A p_B−h_B p_A.

Then ordinary Gaussian multiplication gives

conjugate(z_A) z_B=N+iD,
N²+D²=(h_A²+p_A²)(h_B²+p_B²).                       (2)

The earlier selectively contacted ordinal sources F(A,B) and F(B,A) have counts h_B(p_A+h_A) and h_A(p_B+h_B). Their difference is

H(F(A,B))−H(F(B,A))=−D(A,B).                         (3)

Thus the Gaussian imaginary coordinate is, up to this explicit sign, the physical recontact defect. Swapping A and B preserves N and negates D.

Important distinctions: tournament converse fixes h and p, so it fixes z; it is **not** Gaussian conjugation. Ordinal source composition sends (h,p) to (h_Ah_B,p_Ap_B), not to a Gaussian product. Multiplying Gaussian coordinates cannot be assumed to produce another physical profile.

## 2. The actual profile cone and the mod-four barrier

The two covers with one slot empty contribute2h to P₂; all other contributions are nonnegative. Therefore

p>=2h,

with equality only for a singleton. If the order is at least2, complementary subsets pair in the sum for P₂. Each induced H is odd and there are an even number of complementary pairs, giving

p≡0 (mod4).

Thus for two sources of order at least2:

D≡0 (mod4),  N≡h_Ah_B (mod4).                       (4)

A diagonal norm from such a source is1 modulo8, since h is odd and p is divisible by4. The singleton diagonal gives5. This is stronger than the unrestricted Gaussian norm condition1 modulo4.

The cone also yields

N>=5h_Ah_B,       |D|<N/2.                          (5)

For example if D>0, D<h_Ap_B<=p_Ap_B/2<N/2; exchange A,B for the other sign. These are restrictions on this reader, not on all tournaments. They make the target inverse finite and decreasing in the scalar source counts, but do not construct the required joint profiles.

The classical split primes13=3²+2² and29=5²+2² are not physical diagonal outputs: their even coordinates2 violate p>=2h. More generally all primes5 modulo8 except5 are excluded by(4). Even splitting and the allowed mod-eight class are insufficient:41=5²+4² also violates the cone, although41 is1 modulo8 and is a realizable tournament count by other constructions.

## 3. Split, inert and ramified primes: exact scope

In the Gaussian integers,2 ramifies:2=−i(1+i)². For an odd prime ell, the quotient by ell is governed by X²+1: it splits when ell≡1 mod4 and is a field when ell≡3 mod4. This is the familiar split/inert distinction; it is arithmetic, not a classification of tournament symmetry.

For an inert prime ell≡3 mod4,

v_ell(h²+p²)=2 min(v_ell(h),v_ell(p)).               (6)

Proof: divide out the common ell-power. If the remaining sum of squares vanished modulo ell with one coordinate nonzero,−1 would be a square, a contradiction. Thus a primitive Gaussian norm has no inert prime factors. Nonprimitive norms can contain inert primes, but only to even exponent.

For the mixed read(1), an inert prime may divide N without dividing D. If ell divides **both** N and D, equation conjugate(z_A)z_B=0 in the field F_(ell²) forces either z_A=0 or z_B=0 there: one source has both coordinates divisible by ell.

If N itself is an inert prime, neither source can have that common factor, since0<h_A,h_B<N by(5). Hence D is nonzero modulo N, and in fact0<|D|<N/2. A zero-determinant or perfectly aligned reader therefore cannot produce an inert prime, whereas a mixed, misaligned reader can.

For any prime N in(1), both physical profiles must have gcd(h,p)=1: each gcd divides N and is strictly smaller than N. Primitive profile content is necessary, not sufficient, and is not equivalent to having trivial automorphism group.

The ramified coordinate1+i does not supply a parity-preserving tournament move. It takes h+ip to(h−p)+i(h+p), whose coordinates are both odd and whose real part is negative in the physical cone. Its norm is even. This violates the physical chart rather than constructing an even or shifted odd Hamilton count.

## 4. Physical examples in both prime sectors

The standalone replay gives full adjacency rows and checks strong connectivity and exact counts.

* A=B=TT₂, profile(1,4): N=17,D=0.
* A=B the four-vertex source with rows[8,1,3,6], profile(5,32): N=1049,D=0.
* A singleton, profile(1,2), B the five-vertex source[24,1,3,6,14], profile(11,100): **N=211,D=78**. Here211 is3 modulo4 and is therefore not a Gaussian norm, but it is the real part of the polarized product. Its actual eight-vertex tournament rows are

  [190,176,130,134,140,156,63,64].

  Check:211²+78²=5·10121. The imaginary coordinate is not divisible by4 because one source is a singleton, precisely the boundary excluded in(4).
* A=TT₂ and B the profile(5,32) source give **N=133=7·19,D=12**. Thus an inert prime factor to odd exponent does not obstruct the mixed reader. Check:133²+12²=17·1049.

Taking A=TT_m and the fixed profile(11,100) source yields the all-order strong family

N_m=11+100·2^m,       m>=1,

on m+7 vertices. Every count is3 modulo4. The first two prime examples are211 (m1) and811 (m3); even m forces divisibility by3. No infinitude of primes in this family is asserted. Its sparse dyadic addresses do not select an arbitrary prescribed prime.

## 5. Why source content must be broken physically, not divided out

The earlier rigid source C=[16,1,3,7,14] has profile(9,84), gcd3, and trivial automorphism group (verified by all120 permutations). Its diagonal read is7137=9·13·61, and every mixed dot read using it is divisible by3. Thus rigidity does not eliminate the arithmetic content barrier.

Nevertheless the selective nonuniform returns of `niche_bridge/EQUAL_GAIN_RECONTACT.md` turn a host using this same source into131 and133. The operation genuinely changes the read and breaks the common divisor. Nothing was divided out numerically. The equal-slot square-reflection source subsequently moves the window to141/143.

This is the useful constructive lesson: the Gaussian/profile obstruction tells us **which read needs to change**. It cannot prohibit other physical constructions. In particular7 and21 are not singled out by Gaussian inertness:3,11,19 and the realized count63 also fail the sum-of-two-squares test. Treating norm representation as necessary for Chan–Pak would reject many actual witnesses.

## Remaining selector

For an exact target N, the polarized reader asks for two actual profiles in the cone with h_Ah_B+p_Ap_B=N. Equations(2),(4),(5),(6) prune candidates and show when alignment or source content is incompatible, but do not supply them. If this finite inverse has no physical solution, change the quotient or use a selective contact that reads marked data, as the earlier return experiment does. Congruence success and Gaussian splitting alone are not absolute placement.

## Reproduction

Run `python gaussian_profile_replay.py` and `python -O gaussian_profile_replay.py`; outputs are byte-identical. The script computes induced path tables, P₂, literal shell-of-ordinal tournaments, exact Gaussian identities, strong connectivity, source automorphisms, inert valuation identities on the displayed sources, and six terms of the all-order family. The all-order proofs above do not depend on these finite checks.
