# A rank-four controller at arbitrary binary and decimal precision

**Proved.** For any positive integers s,d and any prescribed stream of odd classes modulo 2^s·5^d, actual literal span-3 return blocks can realize that stream online at chosen boundaries. Rank stays four at every letter and the **same** rational image returns after each eight-letter block. Precision may also increase between requests. No claim is made about a short uniform block budget as precision grows, strongness at intermediate zero letters, or prescribed absolute Hamiltonian counts.

Use the [168 certified return blocks and integral four-coordinate lattice](https://github.com/eliottcassidy2000/math/blob/3a2423f0bb11b8a63dc81b7b3616a881c101b3f6/05-knowledge/results/chanpak_eight_braid_20261116/THREE_LEASHES.md). Let A_i be their integer lattice actions (column convention) and f=(32,14,10,27) the H row. Each physical block is the three listed steering letters followed by the same five refresh letters. Each maps the rational image E isomorphically to E. In this document only 139 actions whose determinants are units mod 5 are used; the initial state z=(-43,-5,4,53) is nonzero mod 5, so this invariant is preserved.

## The two hands

One hand holds a one-dimensional dyadic state fixed while the other edits four five-adic coordinates. Its implementation matters: independent coverage at the two moduli would not alone imply a joint result.

### Dyadic projection without losing rational rank or five-adic primitivity

Let T=A_(4,3,5), where the subscript is the steering triple. Its exact characteristic polynomial is
```
x^4 − 905 x^3 + 1960 x² + 9424 x + 3456.
```
Modulo 2 this factors as x³(x−1); T is invertible modulo 5 and T^620=I modulo 25. Work at s≥3 (this proves smaller dyadic precisions too). Choose the smallest positive K≥3s with K≡16 modulo 620, so K<3s+620. Put R=T^K. Then R≡T^16 modulo 25, and K is a multiple of four.

Modulo 2^s, R sends every state onto a rank-one free line (on the line it may multiply by a unit). Here is a direct proof. The characteristic polynomial has a simple odd root modulo 2 and it lifts uniquely to an odd root λ modulo 2^s by adjusting one digit at a time. Factor p(x)=(x−λ)g(x); g(λ) is a unit and g≡x³ modulo 2. Cayley–Hamilton gives the idempotent Π=g(T)/g(λ), with T=λ on its image and g(T)=0 on its kernel. On the kernel T³ is twice an endomorphism, so T^(3s) vanishes there modulo 2^s. On the odd eigenspace R multiplies by the unit λ^K; thus R=λ^K Π sends every state onto the simple odd eigenspace. Modulo 8 its scalar is one because K is a multiple of four, so its reduction agrees there with T^16. For the free rank-one assertion, choose a lift v in its image generating the one-dimensional reduction modulo 2. For any y in the image, subtract its appropriate mod-2 multiple of v. The difference lies in 2·image(Π), since applying Π to a vector divisible by 2 leaves its image divisible in that image. Iterating proves that v generates the whole image. Here v can be normalized to f v=1: at mod 8 its normalized reduction is (4,0,0,3), verified directly from T^16.

Apply R to begin a request. The resulting physical state has H odd and is of the form H v modulo 2^s. For each usable A_i, employ the macro B_i=R A_i (first A_i, then the repeated physical reset). On the dyadic line it multiplies by a scalar. At mod 8 these scalars include 3 and 5. A lift of 5 generates the subgroup 1 mod 4 modulo every 2^s: v₂(β−1)=2 and successive squaring raises this valuation by one. A lift of 3 lies in the other coset. Thus two chosen macros can set the dyadic request to any odd value multiplicatively. They may change the five-adic coordinates; they do not destroy primitivity.

### Five-adic corrections that leave that dyadic choice alone

For all usable B_i take G_i=B_i^e for e=2^max(s−2,4). Their dyadic multipliers are all one because the exponent of the odd units modulo 2^s divides 2^(s−2). Calibration uses the fixed matrices
```
G_i^0 = (T^16 A_i)^16 mod 25.
```
It checks two exact facts: their directed orbit modulo 5 contains **all 624 nonzero vectors**, and with o_i=ord(G_i^0 mod 5), the 139 matrices
```
M_i = ((G_i^0)^o_i − I)/5 mod 5
```
span the **16-dimensional** space of all four-by-four matrices over F5. The verifier enumerates the finite orbit and performs rational-free elimination on 16 coordinates modulo 5. Since these actions are permutations on nonzero vectors, their directed orbit is the same group orbit from every one of the 624 points. Their mod-5 orders are odd: the exponent 96,720 used in the calibration has 2-part 16, and the sixteenth powers remove that part.

For any s, the actual G_i modulo 25 is the corresponding G_i^0 raised to a further power of two. Modulo 5 that does not change any generated cyclic group or o_i. Its first-order difference ((G_i)^o_i−I)/5 is the corresponding nonzero scalar multiple of M_i modulo 5, so the tangent matrices still span. This is why a single finite calibration suffices for every precision, not a finite coprimality experiment at a few sample levels.

Choose a desired primitive vector w modulo 5^d with f w equal to the odd-part target readout. Such a w exists for *every* target, including multiples of 5: fix its fourth coordinate to one and solve for the first, since f₁=32 is a unit. First steer the current nonzero vector to w modulo 5 with the finite orbit. When the current vector agrees with w modulo 5^j, use powers of G_i^o_i. The binomial identity gives
```
(G_i^o_i)^(5^(j−1)) ≡ I + 5^j M'_i  (mod 5^(j+1)),
```
where M'_i are the calibrated tangent matrices rescaled by units. Since the current vector mod 5 is nonzero, the vectors M'_i w span F5^4. Products of these near-identity maps therefore select any next correction (their cross terms vanish at this digit), preserving previous digits and the dyadic readout. Induct through j=d−1. This proves the joint exact residue steering.

All advertised operations are products of physical eight-letter returns, including exponentiation by positive repetition, so they preserve the same rational image at block boundaries and rank four at every intermediate letter by monotonicity of rank. A later request may use a new K and precision; the preserved primitivity and physical oddness let it start again. The reset is rank one only in the dyadic **quotient**. Keeping the scalar motion instead of forcing a true idempotent saves a long power: the reset length is O(s+1), rather than exponential in s. The constructive request budget is O((s+1)·2^s·5^d) letters since the finite mod-5 path and all exact mod-5 generator orders contribute only constants, and the lifting lengths form a geometric series.

## Replay and scale

Generate `return_automaton.json` with the earlier companion's `lattice.py`. Run `python mixed_controller.py` from this directory. Its finite calibration uses integer arithmetic and explicitly checks the mod-5 orbit and the tangent rank. It then constructs compressed physical instructions for three successive targets modulo 128,000, verifying their effects by modular multiplication. That modulus is lcm(2^10,10^3); it lets binary and decimal prescriptions coexist iff their common mod-8 constraint agrees. The useful fact about 1024=1000+24 here is not equality of scales: their low mod-8 alignment does not fix the missing powers of 5, which the second hand really controls. The witness plans may be long when expanded and do not land at the small absolute targets used as labels.

No Riemann hypothesis, prime-distribution conjecture, or floating-point spectral assertion is used.
