# Magic-square balance as a Chan–Pak repair and target test

Two archived pieces are used directly. [Session57](../session57_square_sum_magic_paths_20261025/README.md) separates the decorated additive frame (the labels and diagonals really matter) from simply having a path or a list of eligible values. Its congruent-curve bridge retains three arithmetic progressions at one fixed gap. [Session143](../session143_magic_square_cyclic_bridges_20261103/README.md) classifies real cubic rearrangements on nine consecutive indices: the nonlinear survivor has deviations proportional to k³−15k. Those are implemented as literal tournament readings below. See also the primary [Bremner article](https://pldml.icm.edu.pl/pldml/element/bwmeta1.element.bwnjournal-article-aav88i3p289bwm) and [geometric formulation](https://arxiv.org/abs/2609.09351) for the ordinary square-entry problem. We prove neither that square-entry problem nor Chan–Pak.

## A true three-choice repair gives a physical magic frame

The marked host with outgoing rows

    (52,365,176,69,266,280,437,315,13)

has three choices of tag mask, 113,43,39. Tags specify the vertices beaten by a new root. All three have the same profile (H,P₂,W,Z,E,S)=(1521,26412,10422,11997,750,631), but K=H−U takes −2387,−2385,−2383. Here W,Z are weighted marked starts/ends, E,S ordinary marked ends/starts, and U the joint marked two-cover count. A transitive return module of order k≥1 therefore reads 10422·2^k−2387, and two and four more, all from genuine strong tournaments. At k=1 the three values are 18457,18459,18461.

Under ordinal stacking an isolated change in local K scales by the product of other H values, with no mixed effect because all other profile coordinates are held fixed. With four copies of the new ternary tile, make one copy a free digit and the other three a coherent free digit. Each may be in any of its three states. The resulting nine exact counts are C+d(i+3j−4), 0≤i,j≤2, d=2·1521³. Arrange them in the ordinary Lo Shu order, or in the following balanced frame, where a=d,b=3d:

    C+a     C−a−b   C+b
    C−a+b   C       C+a−b
    C−b     C+a+b   C−a.

The eight lines are exactly 3C, every graph is strong, and the nine values are distinct. Strength follows because the root reaches the first marked block, which reaches each later block, and all vertices return through the terminal module. The underlying host is unchanged throughout, aside from literal root contacts on the independent tile groups. The progression has step d, **not** step two; the old scalar count cannot replace the marked interface.

For a literal implementation of Session143's nonlinear survivor, take eleven copies and use coherent groups of nine and two. Then the normalized deviations are 9x+2y, x,y∈{−1,0,1}. They are the exact multiset of (k³−15k)/2 over k=−4,…,4. In the balanced frame use a=9d,b=2d (this time d=2·1521^10). This realizes the archived cubic pattern as an exact magic square of strong tournament counts at every return height. Entrywise squaring would not preserve the balance and is not attempted.

## Two choices cannot carry the same additive frame in a coprime bank

In any 3×3 magic square over a ring with 3 invertible, every opposite pair sums to twice the center; the same identity also holds over the integers and thus modulo every m by reduction. Over any odd modulus, a magic square that takes at most two residue values must be monochromatic: the value opposite any noncentral residue must be its reflection about the center, neither of the two choices unless the residues coincide.

Use a fixed-height ordinal read with pairwise coprime odd tile counts h_i>1. Assume each tile offers only two K choices, distinct modulo h_i, and all other profile coordinates are fixed. Its effect on the result modulo h_i is detectable and all other tiles disappear in that residue (their difference weights have factor h_i). An integer magic grid from this bank must therefore use the same choice in all nine cells for each tile. This rules out any nonconstant grid from such a two-choice construction. More generally with a three-residue palette over a prime at least five, if there is any variation, the palette must be {C−d,C,C+d}. For the balanced frame this forces one of a,b to vanish modulo that prime, since otherwise one of a+b or a−b is ±2d. These are precise palette conditions, not just limitations on the number of witnesses scanned. They suggest splitting prime tags between the two physical axes rather than assuming a binary bank can span the desired exact window.

## Putting a prescribed number at the center: an inverse, not surjectivity

Take the full-boundary binary switches at counts 389 and 97 from the [Bézout repair](../chanpak_prime_dichotomy_20261117/BEZOUT.md). Stack two copies of each, taking the center with one of each switched; zero, one, or two switches of each type gives a balanced frame. The common deletion profile at the center is

    J (ordinary H)=1423779289,  P=205520896000000,
    W=46594073974899,           Z=142977375089376,
    K=−32411071075897.

The two axes step by 8J/389 and −2J/97. They are independent in the 3×3 range since 4·97 and 389 are coprime and exceed 2. Their sum is the primitive lattice step −2·37733 because 389=4·97+1, echoing the earlier Bézout move without claiming a unit two-count repair.

Prepend a transitive contact word encoding an even binary address v, with empty word at v=0. Its weighted start is v, ordinary H=1, and full marked end is zero. Its intercept is F(v)−v, with F(0)=1, F(2t+1)=F(t), F(2t)=F(t)+t+(t mod2) for t>0, as in the [earlier prefix method](../chanpak_target_passports_20261117/PREFIX_INVERSE.md). The common step sizes stay fixed (prefix ordinary H=1). At return height one the exact center is

    C(v)=C0+A v+J F(v),  C0=2W+K−J,  A=2P−J−Z.

For an absolute target N, choose the least power of two L with AL≥2√N. Write v=Lq+r, 0≤r<L, and let X_r be the reverse-complement of the ell-bit word r, where L=2^ell, and delta be one if its first bit is zero (both zero when ell=0). The recurrence gives F(v)=F(q)+X_r q+delta(q mod2)+F(r)−1. Put B=AL+JX_r and Y=N−C0−Ar−J(F(r)−1). A candidate must satisfy

    Bq≡Y mod J,
    max(0,ceil((Y−2J)/(B+J))) ≤ q ≤ floor((Y−J)/B).

This uses 1≤F(q)+delta(q mod2)≤q+2. Since gcd(A,J)=1, B is invertible modulo J. The real interval's width is at most JN/B²+J/B<J (for eligible N), so each suffix r leaves at most one q. Check v even and the exact equation. This is an exact O(√N) arithmetic-candidate test for the center, or for any designated cell after subtracting its fixed offset. Possible candidate existence for every N is not implied.

For the local/global contrast, [the separate residue construction](SQUARE_PHASE.md) produces nine distinct, literal magic counts all congruent to a square modulo any prescribed modulus in its stated class. For that family the nine-entry width stays fixed as the center grows, so at sufficiently large heights **at most one** entry can actually be a square. This models exactly why passing simultaneous phases does not choose a physical exact target.

Replays: `python replay.py` checks literal tile counts, the cubic and linear magic assemblies, and absolute-center inverse controls. `python replay_square_phase.py` checks the separate congruence construction.
