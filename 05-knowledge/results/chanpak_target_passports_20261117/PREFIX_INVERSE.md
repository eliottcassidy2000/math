# An exact binary-prefix decoder for a stackable sevenfold twin family

The [independent-digit transfer](ORDINAL_DIGITS.md) gives a genuinely stackable +14 switch and permits a free transitive prefix. This yields a physical two-parameter family that can now be inverted from an absolute target, not just approximated. It does not land at all targets.

## The same boundary, a different absolute correction

Take the deletion R=(50,104,19,21,98,76,13) and the two tag masks 31 and 47. The full boundary (H,P₂,weighted starts,weighted ends,full ends,full starts) is identical:

    (175,1988,1295,1295,125,125).

The intercepts K=H−U are −609 and −595; the difference is 14. Any common ordinal stack sees identical mixed contacts. Without a prefix the upper quotient is always divisible by 5.

Prepend a transitive tournament whose tags encode the binary integer v, the first vertex being the most significant bit. Its weighted start is v. Its full marked end count is v mod 2, so take **even v** (v=0 uses the empty prefix). If F(v) is the single-apex path count for this contact word, its intercept is F(v)−v, by the ledger or the explicit terminal count of 1. The ordinal transfer gives W=1988v+1295 and lower K=175(F(v)−v)−784−1295v. For odd v an extra term −125 appears in K and prevents divisibility by 7; it cannot be removed by the power of two in the output since W is already divisible by 7.

Append the outer transitive module TT_k, k≥1. The counts are 7Q(v,k) and 7(Q(v,k)+2), where

    Q(v,k) = C+A v+25F(v),
    C=185·2^k−112,   A=284·2^k−210,
    F(0)=1,  F(2t+1)=F(t),  F(2t)=F(t)+t+(t mod 2)  (t>0).

The recurrence is the [one-apex cut count](GAP_WORDS.md#exact-read-and-concatenation) read from the least significant bit. The graphs are strong for v>0 because q reaches the leading tagged vertex, which dominates every later block; all exterior vertices return via the outer module. The v=0 case is checked literally and uniform inflation preserves it. The known twins 1049,1051 occur at v=2,k=1, now as an instance of a freely selectable prefix with full boundary equality. Other verified twin quotients include 12539/12541 and 268529/268531. The +14 shift is exact, not an assertion of indefinitely many prime pairs.

## Choose a candidate address from the count

For a prescribed lower quotient p, only heights with 185·2^k−87≤p matter, since F≥1 and A>0. At each height gcd(A,25)=1 (its residue mod 5 is 4·2^k), so v is fixed modulo 25, and evenness fixes it modulo 50. This alone cannot force the exact target.

Let L=2^ℓ be the least power of 2 with AL≥2√p. Write v=Lq+r, 0≤r<L. For ℓ>0 let X_r be the ℓ-bit reverse complement of r and δ=1 if its top bit is zero. For ℓ=0 set both to zero. Iterating the low-bit recurrence gives

    F(Lq+r)=F(q)+X_r q+δ(q mod 2)+F(r)−1.

The contributions at a low zero are its remaining prefix plus the next bit. The prefix contributes q with the appropriate power of two; the next bit depends on the parity of q only at the top low position. Summing those contributions gives the identity and an exact carry between low and high blocks.

Set B=A L+25X_r and Y=p−C−Ar−25(F(r)−1). A solution has

    Bq ≡ Y (mod 25),
    max(0,ceil((Y−50)/(B+25))) ≤ q ≤ floor((Y−25)/B).

Indeed 1≤F(q)≤q+1 by binary induction, so the remaining bracket F(q)+δ(q mod2) lies between 1 and q+2. The real interval width is at most `25/B+25p/B²<25` for any eligible p (p≥283), since Y≤p and B≥2√p. Also B is invertible mod 25. Therefore at most one candidate in each low cell survives the label. Try it, check v even (skip odd r if ℓ>0), and substitute the **full** equation. An upper-rail target just calls the same inverse at p−2.

The total number of low cells over all eligible heights is O(√p+log p), because A grows geometrically. Each cell uses O(log p) arithmetic steps on O(log p)-bit integers. At heights where A≥2√p there is only one low cell, giving a polylogarithmic exact decision for that large-height portion. This is an inverse of this two-parameter family, not of the entire spectrum.

Replay: `python replay_prefix_inverse.py`. It verifies both physical witness counts by subset DP, the boundary profiles, the low-bit formula, and comparison to exhaustive arithmetic addresses at small targets, with additional larger controls. The proof is not inferred from the tests.

## Every congruence phase at one fixed height

There is an independent exact **modular** result for this same two-variable reader: for each fixed k≥1, every residue compatible with oddness modulo any M occurs as a sevenfold quotient. The decoder in the previous section is still needed for absolute equality.

To see it, put v=2r and `φ(r)=r+1−(r mod 2)−F(r)`. This is even (F is odd), has φ(0)=0, and satisfies

    φ(2r)=φ(r),   φ(2r+1)=φ(r)+r+(r mod 2).

Now restrict further to r=2u, and let g=φ(2u)/2. On this address space two literal bit operations keep r even:

    A: append 0 to r       (u,g)↦(2u,g),
    T: append 10 to r      (u,g)↦(4u+1,g+u).

In these coordinates the lower quotient is `C₀+2Aₖu−50g`, with C₀=185·2^k−87 and Aₖ=568·2^k−370. Thus it is always odd.

First take an odd modulus m. Both A and T are invertible modulo m; inverse moves can be implemented by positive repetitions. Composition below is right to left. Set D=T A⁻² and S=A D A⁻¹. Directly,

    D: (u,g)↦(u+1,g+u/4),
    S: (u,g)↦(u+2,g+u/8),
    S⁻¹D²: (u,g)↦(u,g+3u/8+1/4).

At u=0 the last operator adds the unit 1/4 to g. One can set g there, then apply D to reach any u while precompensating its known drift. Every pair (u,g) modulo m is physical. The quotient coefficients 2Aₖ and −50 are primitive modulo any odd m: their only possible shared prime there is 5, and 5∤Aₖ. So the resulting quotient can have any specified odd-modulus residue.

To preserve that phase while controlling the dyadic part, write the requested modulus as `2^(a+1)m` for a≥0 when the modulus is even. Choose j≥a with `2^(j+4)≡1 (mod m)` (any j≥a if m=1). The physical reset R=T A^j is invertible modulo m; choose the odd program using the preimage of the desired state under R. Modulo 2^a, R puts u at 1 while leaving g unchanged. Now use the physical loop L=T A^j T. On the dyadic state u=1 it returns u to 1 and increments g by one. Modulo m its u-coefficient is 1, so it acts by `(u,g)↦(u+b,g+c u+d)` for constants b,c,d. Because m is odd, m repetitions are the identity modulo m: both the terms linear in m and the drift proportional to the binomial m(m−1)/2 vanish. But L^m increments g by the odd unit m modulo 2^a. This lets the target set the dyadic part independently, since after subtracting the odd constant C₀ and dividing by two, the coefficient of g is −25.

This supplies a constructive phase program at fixed k without ruling out any modulus by a stock-prime coprimality condition. It does not bound the magnitude against the prescribed target. The replay `python replay_prefix_phase.py` checks the constructive word maps and the CRT combination for mixed moduli; its program only appends literal bits to a transitive prefix.


The dyadic correction need not take a number of iterations proportional to the dyadic modulus. At the dyadic state u=1, first apply A^i and then T, which adds 2^i to g. Apply A^d and then T to reset u to 1, taking d≥a and i+d+4 divisible by the order of 2 modulo m. The reset adds zero modulo 2^a. The combined loop is again a translation in u modulo m, with drift affine in u, so its m-fold repetition is odd-phase neutral and adds m·2^i to g modulo 2^a. Use the bits of the desired dyadic correction divided by the unit m. This yields O(m a(a+ord_m(2))) letters for the dyadic correction, only O(a²) for a pure dyadic modulus. It does not give a comparable magnitude bound for the entire program.
