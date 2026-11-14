# Nonlinear joint source selection modulo every modulus coprime to six

**Theorem.** Given M coprime to 6 and arbitrary residues h,q modulo M, there is a literal nonempty tournament T with (H(T),Q_T(2))≡(h,q) modulo M. Both residues may be zero or nonunits. The chart uses ordinal joins of the two explicit one-apex patterns from [Session 322](../session322_target_inverse_prime_channels_20261113/RESIDUE_APERTURE.md) and transitive singleton padding. The ordinal construction is followed below by a single cyclic wrap giving a **strong** realizer for every such residue pair. No assertion is made that (h,q) is realized as an absolute pair.

## Local normalized profiles

Let P_d(n) have transitive chain of length dn and apex cut (10^(d−1))^n. We use d=2,3. The inherited odd-cycle calculation gives
```
H = 1+n+(2^(d−1)−1) Φ_n(2^d),
Q(2)/2^(dn) = 2 + A_d Φ_n((3/2)^d) + B_d Φ_n(2^(−d)),
Φ_n(r)=n/(1−r)−r(1−r^n)/(1−r)^2,     Φ_n(1)=n(n+1)/2,
A₂=B₂=1/2, A₃=5/4, B₃=3/4.
```
Take n in a progression of step D coprime to p, divisible by p−1, and with n≡0 mod p−1. The modulo-p powers are then trivial. At p≥7 use d=2; with prescribed H residue x, its normalized second coordinate f(x)=Q(2)/2 is (x+4)/5 mod p. At p=5 use d=3, with f(x)=(3−x)/2. Also d=2 at p=5 gives f(x)=(x²+x+1)/8, always nonzero modulo 5.

For both main choices of d, the H slope in n is a p-adic unit, and H may be prescribed independently as a parameter x modulo p^e. Indeed a shift of n by j D p^t changes each exponential by zero modulo p^(t+1), while the H linear term changes by a unit times j p^t. The same increment calculation for Q(2)/2, including the 2^(dn) factor, shows its change with x is the slope of the displayed affine f modulo p at every lift level. The d=2 chart at p=5 also permits H to be prescribed: its H formula has unit denominator 9 and increment slope 2/3, even though its Q coordinate is quadratic at the first digit.

## Independent H and Q digits via a nonlinear pair

Use ordinal products, so both coordinates multiply. Work first modulo one prime p and use normalized second coordinate Q(2)/2 per letter. At p≥7 the nonzero normalized Q values obtained from two letters with H product 1 form the set
```
A = { [17+4(x+x^−1)]/25 : x∈F_p*, expression nonzero }.
```
The involution x↦x^−1 has two fixed points, so its image before deleting zero has (p+1)/2 distinct values. Zero appears exactly at x=−4,−1/4. Thus A has half the nonzero field elements. It is not a coset of an index-two subgroup. Indeed the sum of its elements is 17/50: the x+x^−1 terms sum to zero over these inversion classes, since their two fixed values cancel. This sum is nonzero except at p=17, while a nontrivial multiplicative subgroup coset has sum zero. At 17, A contains both 1 and the nonsquare 14 (x=2), so it is still not a coset. Therefore A·A has strictly more than half the group: equality of cardinality would make all translates aA for a∈A identical, forcing A to be a subgroup coset. It follows A·A·A is the entire unit group by intersection of sets of complementary sizes.

At p=5 the two d=3 letters with H product 1 give normalized unit values 1 and 4. A mixed pair supplies a nonsquare: use d=3 at x=2 and d=2 at x=3, so H product 1 and normalized Q product 3. These give all units by products. Whenever lifting such a correcting pair, prescribe its second H as the exact modular inverse of the first, not merely the inverse modulo p.

To lift a chosen Q residue, reserve one pair of main affine letters. When the target Q is a unit, take its H product b to be a nonsquare mod p and choose nonzero x with both f(x),f(b/x) nonzero. Such x exists (at most two excluded choices), and the first-digit derivative of f(x)f(b/x) is a nonzero constant times 1−b/x². When target Q is zero mod p, choose b=1 and x equal to the root of f (−4 or 3), which is neither 0 nor ±1. The pair then has Q zero mod p with the same derivative condition.

Supply the remaining required H after dividing by b. For a unit use two main letters with Q both units (avoid at most two forbidden first H choices). For a nonunit use one d=2 letter, whose Q is a unit whenever x=0 mod p. If Q is to be a unit, multiply in H=1 correcting pairs above to match its first digit. Fix these letters with their H values lifted to p^e. Finally vary x on its selected residue class and set the reserved second H to b/x exactly modulo p^e. Its Q derivative is a unit, and the other letters have unit Q. This supplies each remaining Q digit uniquely while retaining the exact prescribed H product. In the p=5 case, all d=3 letters used have unit H.

## Simultaneous prime-power packets

Process primes of M in increasing order. For each new letter at p, initially require n to be a multiple of p−1 and of r^(f+2)(r−1) for every already processed r^f. This progression has step coprime to p. Its letters have (H,Q(2))≡(1,2) at each earlier prime power: every displayed Φ vanishes there by the formula and the extra two powers, even when a denominator 1−r_base is divisible by the earlier prime. For d=2,3 the possible odd prime divisors of those nonzero numerators are only 3,5,7,19, with valuation at most one. Use the unit-H slope at p to choose the lengths and preserve their resulting profiles by adding multiples of the progression step times p^e.

Before committing them, ensure each letter has *both* coordinates units at every larger requested prime q. In increasing order, fix the current phase mod q−1 and then vary n modulo q without disturbing earlier conditions. In that fixed phase H and normalized Q are polynomials of degree at most two in n modulo q. They are nonzero polynomials: for d=2 their slopes are nonzero; for d=3 H is affine except for its nonconstant quadratic at 7, and Q is affine except for its nonconstant quadratics at 7 or 19 and its slope zero at 11. This last exceptional slope occurs only when d=3 is used in the p=5 packet. Directly in its even phases 0,2,4,6,8 modulo 10, the Q(2)/2^(3n) values modulo 11 are 2,0,10,8,10. The sole bad phase n≡2 mod 10 would give H=0 modulo 5, forbidden in all d=3 letters of that packet. Thus at most four residues of n are excluded for each q≥7, and a good new digit exists.

Previous packets have unit coordinates at the current prime; divide by their values when choosing the requested pair for the next packet. If that packet has K letters, choose singleton padding t so 2^(K+t)≡1 at the current prime power and all earlier ones. This realizes the normalized Q product used in the local proof and makes the packet identity at every earlier prime power. Positivity of the word lengths can always be enforced by adding their retained periods. This proves simultaneous joint surjectivity, with every step represented by explicit tournament factors.

## Why this still cannot select an arbitrary absolute pair

Allow every ordinal multiset of P₂(n),P₃(n) and singletons. The two coordinates only see the multiset, not the order. Since Q(2)≥2^order, pairs with Q(2)≤X have order O(log X). At each component size there are at most two pattern types plus a singleton letter. The elementary colored-partition bound gives at most exp(O(sqrt(log X))) distinct pairs below that Q budget. Thus even prescribing *both* responses simultaneously modulo arbitrary odd moduli prime to 3 cannot turn this sparse bank into an absolute selector.

Replay `replay_joint.py` uses division-free slot transfer, constructs the actual physical component lengths, and checks prime powers and mixed moduli including the nonunit q and p=5/q=11 avoidance boundary.

## Strong inflation upgrade away from 5

If M is coprime to 30, the same arbitrary joint residues (H,Q2) can be realized by a **strong** tournament. Construct an ordinal helper V with H(V) and Q_V(2) units modulo M, but Q_V(3)=0 modulo M. For each p^e dividing M use one P₂(n). For z=3 its normalized slot response is

    Q(3)/3^(2n+1) = 1 + [Φ_n(16/9)+Φ_n(4/9)]/9.

At p≥11 choose phase n=0 modulo p−1. Its slope modulo p is 2/35, so n=−35/2 modulo p gives a simple zero; at that root H=−32/3 and Q2=−8/3 modulo p, both units. Higher digits lift by the same affine-slope argument. For each smaller prime q of M, insist n be divisible by q and q−1, so its H,Q2 are units there; these constraints are coprime to the current p. For a larger prime, freeze its phase and choose a digit avoiding the two affine zeros as before.

At p=7 take phase n≡4 modulo 6 and first digit n≡3 modulo 7 (representative n=10). Direct evaluation gives Q3=0 and H=Q2=2 modulo 7. The response retains a simple root: the first ratio r=16/9 is 1 modulo 7, the other ratio is 2. For the first term use
```
Φ_n(1+s)=Σ_(j≥0) binom(n+1,j+2) s^j.
```
This follows by summing the geometric binomial expansion and hockey-stick twice. When shifting n by D·7^t with t≥1 and D coprime to 7, its terms at j≥1 change by multiples of 7^(t+1): use Vandermonde and v₇(binom(D·7^t,i))≥t−v₇(i) for 1≤i≤j+2, then j−v₇(i)≥1. The j=0 term thus has slope n+1/2 modulo 7, zero at n≡3. The other Φ term has slope 1/(1−2)=−1 along the retained phase, hence the total slope is a unit. This lifts all remaining digits. Subsequent primes are larger and their H,Q2 unit avoidance applies unchanged.

Ordinal multiplication combines these helpers. Given desired output residues A,B, set q=(B+2A)/6 and h=A−q modulo M. Independently prescribe a tournament T with residues (h/H(V),q/Q_V(2)) using the joint ordinal construction. The ordinal source U=V→T then has (H,Q2,Q3)=(h,q,0) modulo M. Its cyclic wrap C3[U,1,1] is strong and reads

    H_out=h+q,   Q2_out=4q+2Q_U(3)−2h,

which are exactly A,B modulo M. This proves the strong joint result for all such M. At 5 this particular helper cannot work: Q(3)≡(−1)^order Q(2) modulo 5 by the parity of Q. The earlier ordinal theorem still includes 5.

## The strong upgrade at 5 by a two-variable tangent, not annihilation

The obstruction at 5 to killing Q3 is specific to that helper. A different local lift gives the **strong theorem for every modulus coprime to 6, including powers of 5**, as follows.

Use the p=5 pair construction but reserve two parameters: the first H, x, of the derivative pair whose H product b is fixed, and the H coordinate u of one ordinary letter in the H-correction part. If the prescribed source H is a nonunit, this ordinary letter is the single d=2 factor with H≡0; if the H is a unit it is one of the two ordinary d=3 factors. Its Q is always a unit at the starting digit, and all other ordinary letters have unit H. The derivative matrix of the source pair (h,q) in parameters (u,x) is triangular and invertible: the h derivative in u is a product of units, the h derivative in x is zero, and the q derivative in x is the already established good-pair derivative times a product of units. This works even when both source coordinates vanish modulo 5. The u parameter stays in its chosen first-digit class; all d=3 factors still have unit H.

Make the total number of factors plus singleton padding a multiple of φ(5^e). All lengths have n divisible by four, so every pattern order dn+1 is odd and the source total order is even; the same padding normalizes the total 2 and 3 factors at 5^e. At the first digit one has Q3=q because 3≡−2 and the source order is even. Thus at the requested output first digit A,B choose source q=(B+2A)/8, h=A−q, and build the base pair as before.

At every higher 5-adic digit, the tangents of normalized Q(3)/3 and Q(2)/2 of a single pattern factor match. Here are the ingredients justifying this stronger statement (plain parity modulo 5 alone would not suffice). At z=2,3 one can write its normalized response as
```
Q(z)/z = z^(dn) [1 + Σ_(j=1)^(d−1) z^(−j−1)
                Σ_(ε=±1) (z+ε)^(j−1) Φ_n(((z+ε)/z)^d)].
```
On the retained phase n=0 modulo 4, shifting n by j·4·5^t at t≥1 fixes the outside powers at the next digit. The Φ increment slope modulo 5 is 1/(1−r) when r is not 1 modulo 5 and is n+1/2 when r is 1, using the same binomial proof as at 7 (the valuation bound is still valid for p=5). Replacing z by −z modulo 5 swaps the two signs, and the explicit signs in z^(−j−1)·(z+ε)^(j−1) cancel. Since the length exponent dn is even, the normalized tangent is unchanged. Reparameterizing by each factor's H divides by the same unit H slope. Consequently multiplication and the fixed padding give ΔQ3=Δq modulo 5^(t+1) for steps changing parameters by multiples of 5^t.

The wrapped pair therefore has tangent (Δh+Δq, 6Δq−2Δh), of determinant 8, composed with the invertible two-parameter source tangent. There is one pair of new digits at each level; solving these two linear equations lifts the requested wrapped pair to 5^e. No higher-digit assumption Q3=q is made.

Finally process the larger primes as in the simultaneous construction, but now make each new packet neutral in **all three** Q values at each earlier prime power. A letter at a later p is (1,2,3) there by the same Φ formula and progression n divisible by r^(f+2)(r−1); the extra normalized slot at z=3 has only the denominators 1−16/9 and 1−4/9 for these d=2 letters, each with at most one valuation at any earlier prime. Padding normalizes both 2 and 3 because its total exponent is a multiple of the earlier φ. At the current prime include its Q3-killer, with unit H,Q2, then correct those two coordinates to the determinant-6 solution. For p=7 its phase n=4 modulo 6 and n=3 modulo 7 is compatible with the initial earlier-5 condition n divisible by 4·5^(f+2); the step is their least common multiple with 6, coprime to 7. Later-prime unit avoidance is unchanged. At the end perform one common cyclic wrap of the whole ordinal source. It is strong and has the prescribed pair modulo M. The replay contains the independent actual-modulus checks of the two-parameter lift and the mixed-prime completion.

The exp(O(sqrt(log X))) bound also holds for the wrapped pattern bank at Q2≤X: the wrap increases order by two and reads only the ordinal component multiset. Thus the stronger joint congruence result does not change the absolute-capacity barrier.
