# The digit bank is a lattice tiling with a missing physical phase

This is an explanation of where the 99-offset construction succeeds and where it cannot, not a global spectrum theorem. Let a denote an odd digit-bank slope (a=99 in S325). For odd H with gcd(H,a)=1, use the odd digit c∈{1,3,…,2a−1}, and write Q=2j and c=2k+1. Then

    N = aQ+cH = H + 2(aj+Hk),       0≤k<a.

## 1. An exact arrangement of the odds into tiles

For all even Q of either sign this map is a bijection to the odd integers: reduction modulo 2a determines c uniquely, then Q is forced. At fixed residue of Q modulo 2H, increasing Q by 2H moves the whole digit progression by 2aH. Its a consecutive points in that residue of N modulo 2H abut the next progression, which starts exactly one progression step later. Thus this is a literal tiling of the odd line by H interleaved arithmetic progressions. No two different Q tiles at fixed H overlap, whether they are mathematically admissible or physical.

If Q≥0 is the only gate, the same arrangement gives the elementary two-generator integer-semigroup conductor. Put m=(N−H)/2. The largest m missed is (a−1)H−a=aH−a−H: in residue of the highest digit its exact j is −1, and every other residue has smaller last miss. Thus the largest odd N missed by the Q≥0 arithmetic shadow is 2aH−2a−H (for a,H>1). Requiring Q≥Q₀ for a fixed even Q₀ translates that front by aQ₀. This front is not a physical tournament threshold. In particular a tempting coincidence, H=151 on a S326 repair chart gives the idealized front 29549, cannot explain that repair's success: the 151 there is the left block H of a different scalar formula and the actual source's Q is not allowed to be zero.

## 2. The first physical phase is a rotation, not a prefix interval

Every nontrivial tournament has Q divisible by four. Temporarily forget the finer phase budget and all other obstructions, keeping only 4|Q. At fixed H the available odd classes modulo 4a are exactly cH, c=1,3,…,2a−1. That is **exactly half** the classes. The longest run of consecutive odd integers allowed by this arithmetic test is

    1 + floor((a−1)/s),   s=min(u,2a−u),   u=H⁻¹ mod 2a with 1≤u<2a.

In fact x=(N−H)/2, multiplied by H⁻¹ modulo 2a, must land in [0,a−1]. Increasing N by two makes x advance by u on that circle. Run the rotation in the equivalent direction of shortest absolute step s. A progression cannot wrap while it remains in the chosen half-circle, so the displayed length is both an upper bound and attained by starting at an end.

For a=99 this predicts run caps of 3,20,3 for H=23,79,91, respectively. H=79 is special here because H⁻¹=−5 modulo198. The often compared 99 consecutive odds belong to the easiest H=1 phase, not to a typical strong source. The exact profiles and the actual per-digit phase can only shorten these arithmetic runs, never lengthen them.

Two complementary H classes H and −H, or H and H+2a, modulo 4a partition these odd residues. In the second pair, the digit c for a prescribed target is unchanged modulo 2a; the two candidate Q values differ by 2c times an odd integer, so exactly one is divisible by four. This is a concrete way to use another source to repair a failed phase without treating phase agreement as full existence. [Equal-budget physical switches](EQUAL_CLOCK_ORBITS.md) supply both complementary pairs at a *literally common* Q.

## 3. Why the 30,030 neighborhood is a place to look

With a=99, a target near 30,000 forces Q below 304 before a digit is chosen. The elementary triangle term in the [small-profile atlas](../session326_exact_source_profiles_20261114/PRIMORIAL_SWITCH.md) excludes order-seven strong components below Q=320. This makes the physical local atlas small even though the absolute target looks large. The first-vertex floor is what makes the local result transferable: at any slope a, below 320a all source SCCs have order at most six. Looking past the floor, the atlas must change; there is no reason a 30,000-scale pattern persists solely because digit residues do.

There is also a purely arithmetic superposition: 30030=66·455 and 198=66·3, so translating by 30030 has order three on the digit ring modulo198. It changes a unit-source digit by a nonzero step of 66 or 132, cycling three digits. This fact does not conserve Q, nor promise that any of the three lifts is feasible. The primorial period governs exponential sieving and the atlas floor governs available source complexity. They are different mechanisms that happen to meet at a convenient scale.

## 4. Positive digits have a real phase price

In the one-return-star notation N=αβ Q+cH with c=h d−kλ, each marked start in I gives at least h legal two-covers by splitting every spanning path before that vertex. Hence α≥|S|h, where S is the nonempty marked set of starts. Also d≤β, as deleting the outgoing constraint only adds covers, and λ≥0. It follows that

    c ≤ hβ ≤ αβ/|S| ≤ αβ.

Consequently if αβ=a·2^r, a digit c>a always costs r≥1. In particular no clever choice of the star's arcs can give a phase-zero digit above 99. If c>αβ/2, the left cut has exactly one marked exit, not two. This is a genuine structural reason high digits appear only after another clock step.

This gives a two-tier arithmetic envelope at fixed H>1 prime to a. For c≤a the effective Q is at least divisible by four; for c>a it is divisible by eight. Modulo 8a the envelope contains (3a+1)/2 odd residue classes out of 4a: there are (a+1)/2 lower digits, each with two lifts, and (a−1)/2 upper digits with one. The classes remain disjoint for different digits by reduction modulo 2a. For 99 this is 149/396, stricter than the first-gate half. It is an upper envelope for the whole one-return-star slope family, not a claim of attainment.

For the literal 99 bank, the digit phase histogram at r=0,1,2,3,4 is (33,34,18,11,3). At fixed H>1 coprime99, the bank's phase gate admits 897 of the 3168 odd residue classes modulo 6336, fraction 299/1056. Because Q≡(−1)^order H mod3, each allowed class keeps exactly two of its three lifts modulo 19008 when imposing Q not divisible by three. The combined envelope has density 299/1584 among odd classes. Just as in the simpler rotation, absolute slot budgets can only reduce this envelope.

The equal-clock orbit provides the whole hierarchy of phase companions: at a fixed digit c, shifting H to H+2a·t leaves c unchanged modulo 2a and shifts its required effective Q by −2ct. As t ranges over 0,…,2^(r+1)−1, these are all even classes modulo 2^(r+2); exactly one pays the digit's first phase-divisibility test. At slope 99 all these H classes can occur at one exact common Q by the orbit theorem. That common Q is not the candidate budget just computed, so the theorem supplies complete phase choices, not exact N coverage.
