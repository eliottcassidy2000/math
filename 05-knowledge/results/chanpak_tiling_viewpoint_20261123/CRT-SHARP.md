# Minimum transitive-precedence alphabets and exact CRT periods

We extend the six-block construction in [PRECEDENCE-FIBER.md](PRECEDENCE-FIBER.md) (commit `e97cb5c`). A block alphabet means a fixed finite set used with repetition under ordinal joins, so Hamilton-path counts are multiplied. As there, each block has a transitive precedence tournament, inherited by every join. The previous result supplies six verified seeds of heights 5,17,37,77,129,751 and orders at most 9.

## Both minima for the original prime set

Modulo 264=8·3·11 the unit group modulo its squares has rank 4. Only seeds coprime to 2·3·11 can appear in unit products. Thus at least four such seeds are needed. The valuation patterns (v3,v11)=(1,0) and (0,1), read modulo 9·121, each require a seed of that pattern, since valuations add. These are two different nonunit seeds. Hence six is optimal over all fixed ordinal alphabets, not merely within the inherited choices.

The maximum order of a seed is also optimally nine: the exhaustive order-eight result in the inherited note and source-padding show every smaller block has H≡1 mod 4; their products cannot hit 3 mod 4. The old six attain both bounds simultaneously.

## Exact period and deterministic normal form

For a≥3 and b,c≥1 put M=2^a3^b11^c and L=2^(a−3)3^(b−1)11^(c−1). Given an admissible residue r, strip the truncated valuations using A=129^min(b,v3(r))·77^min(c,v11(r)). Select a correcting unit u with Au≡r mod M: in a partially saturated prime channel divide the congruence by its prime power; for a saturated channel prescribe u=1; lift each determined unit digit string using least representatives; then apply CRT. This choice is deterministic, including r=0 in any channel.

The order of 5 modulo 264 is 10. The local valuations of 5^10−1 at 2,3,11 are exactly 3,1,1, respectively. Binomial lifting gives orders 2^(a−3), 3^(b−1), 11^(c−1) for 5^10 in the three principal kernels, so 5^10 generates their product, and the exact order of 5 mod M is 10L. The quotient U/<5> has order eight. Projection onto the four binary characters of the inherited note, followed by quotient by the character of 5, maps that quotient onto F₂³ and hence is an isomorphism. The characters of 17,37,751 give a basis. Every correcting unit thus has a unique normal form

    u = 5^t 17^d 37^e 751^f mod M,
    0≤t<10L,   d,e,f∈{0,1}.

Obtain the bits by binary elimination. For the residual 5-power, look up t0 mod10 using modulo264. Divide by 5^t0 and lift the exponent of 5^10 in each local principal kernel, one digit at a time (try 0,…,p−1 and require the next congruence). CRT of these exponents gives t=t0+10j with 0≤j<L. This avoids a search through M residues. The order bound for the actual join is 40L+17+9b+8c; a single source makes the nonsingleton join distinct from its converse if needed. Use a positive period if all factors vanish and a nonsingleton is required. Smaller a,b,c follow by projection. The period is exact, and for units the displayed normal-form exponent bound t≤10L−1 is sharp, irrespective of whether shorter alternative word representations exist.

A lower bound for the necessary common order of this particular bank follows from capacity rather than period. Unit words only use 5,17,37,751; permuting them does not change the height. Thus at order at most N their residue count is at most the number of nonnegative lattice vectors satisfying 4x+6y+6z+9w≤N. Unit-cube comparison with the enclosing weighted simplex gives N^4/31104+O(N^3). As the modulus above has 80L unit residues, a common selector order must be at least (2488320L)^(1/4)(1−o(1)) as L→∞. This does not show the normal-form upper bound is shortest. For individual nonunits one can instead run the unit correction at 2^a·3^(b−min(b,v3(r)))·11^(c−min(c,v11(r))), since saturated or stripped valuation digits impose no further congruence on the unit.

## Stable criterion and three new prime channels

For P a finite set of odd primes, include 2 as the odd-parity channel. Fix an alphabet. Unit seeds must be coprime to every p∈P. For each p∈P, require a seed which has v_p=1 and v_q=0 for all other q∈P. These valuation requirements are necessary by testing isolated valuation-one residues, and sufficient to separately prescribe all truncated valuations provided unit seeds generate.

For unit generation at *all* powers it suffices to test M0=8·∏_(p∈P)p². Indeed local cyclic coordinates at each odd p have orders (p−1)p^(e−1), and the first p² already sees the full elementary quotient at each prime divisor of that order for e≥2. Likewise the 2-adic factor for e≥3 is C₂×C_{2^(e−2)}, with the elementary quotient retained at 8. In a finite abelian group a subset generates iff its images span U/U^ℓ for every prime ℓ dividing its order: otherwise there is a nonzero quotient of prime order. Consequently at least |P|+2 unit seeds are needed from rank at ℓ=2, and at least |P| different nonunit seeds, for total 2|P|+2.

In each of the following *separate* extensions of P={3,11}, keep the old six seeds and insert the two displayed seeds. J_s is the transitive order s with its outer arc reversed; for s≥5 it has transitive precedence order and height 1+2^(s−2), as proved in the inherited note.

| new prime q | new unit height | isolated q-marker height |
|---:|---:|---:|
|13|1025 (J12)|65 (J8)|
|41|257 (J10)|1025 (J12)|
|257|65 (J8)|257 (J10)|

All original markers remain units at q. Each new marker has valuation one at q and is coprime to 3 and 11. The attached checker proves full elementary ranks for all three extensions by listing verified local generator logarithms and performing rank elimination. Thus eight blocks are optimal in each of these three families. No simultaneous three-prime extension is asserted.

Replay: run crt-sharp-check.py with Python 3.
