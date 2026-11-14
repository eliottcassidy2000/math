# The first lost trit is a nearest-neighbor interaction

**PROVED:** leading 3-adic readout for every one-apex marking of an ordinal triangle train, including arbitrary transitive gaps. Closest inherited divisibility is the [arbitrary-gap half-contact theorem](../berggren_target_return_20261111/arithmetic/TRIANGLE_ALL_K_PROOF.md), which gives the universal minimum over extensions. The formula below reads the leading digit for a specified marking, rather than reproving that gcd.

## All marked inputs

Start with an ordinal sum of r directed triangles separated by arbitrary transitive blocks, with optional transitive ends. Empty blocks are omitted. Attach one apex with an arbitrary mask of wins (painted 1). For each triangle in order write r_i in {0,1,2,3} for the number of 1s on it. Call i,i+1 abutting if their intervening gap is empty. Let N be the full Hamilton count. Then for r≥2,

    N / 3^(r−2) ≡ −sum_(abutting i,i+1) r_i r_(i+1)  (mod 3).

In particular the displayed quotient is integral. A nonzero right side characterizes precisely the marked inputs with valuation r−2. If no triangles abut, the next formula is stronger. For each triangle put L_i equal to the last mark of the preceding transitive gap (0 if absent), and E_i equal to one minus the first mark of the following transitive gap (0 if absent). Then for r≥1 and no abutments,

    N / 3^(r−1) ≡ sum_i [2[0<r_i<3] − r_i L_i + r_i E_i]  (mod 3).

These are sums in ordinary integers reduced modulo 3. A single interior transitive gap has two seam marks, so its two contributions cannot in general be replaced by a single arbitrary trit.

## Derivation from physical cut states

Use the [four-coordinate transfer](../chanpak_selective_target_20261114/CUT_STATE.md). Write u=S/h, y=d/h, v=N/h for the normalized old state, where h is the Hamilton count without apex. An appended block of ordinary count a sends

    u'=(q/a)u+s/a,  y'=e/a,
    v'=v+(n_R−a)/a+(t/a)u+(f/a)y.

For a triangle with r_i selected, a=3, q=12, (s,e)=(3r_i,r_i), (t,f)=(3(3−r_i),3−r_i), and n_R−a=2[0<r_i<3]. For example each designated starting vertex has one path of each of lengths 1,2,3 in vertices, and every complement has count one. The directions are read from the marked cut, not independently assigned. For a nonempty transitive block a=1, q=2^(block size), e is its last bit and f is one minus its first bit, since a free nonempty set of remaining vertices cancels in the signed sum except at the extremal endpoint. The other coordinates are integers.

The first physical block has u=s/a, y=e/a and v=n_R/a, so no empty-block transport identity is needed. Thereafter u is integral and y has denominator dividing 3, which can occur only right after a triangle. Multiplying each change of v by 9, only the f*y/a term when two triangles abut can contribute mod 3; it contributes (3−r_(i+1))r_i. The final h is 3^r. In the no-abutment case multiply by 3 instead: the internal triangle change, transitive-to-triangle bias, and triangle-to-transitive bias are the only terms that can survive. This gives the second formula, with the initial block supplying the internal term if it is a triangle.

## The abutting train exhibits asymptotic trit balance

For r≥2 all triangles abutting, choose every vertex mark by independent fair bits. Under this physical sampling the r_i (mod 3) are independent with probabilities (1/4,3/8,3/8). The sum of nearest-neighbor products is asymptotically uniform mod 3, with error O((485/512)^floor(r/4)) in each residue. Indeed divide into disjoint blocks of four and condition on all block endpoints and on all interiors of blocks whose endpoints are not both zero. In a block whose endpoints are both zero, the remaining contribution is the product of two independent trits, insulated from both neighboring blocks. The product has probabilities (7/16,9/32,9/32), so at either nontrivial cubic character its mean is 5/32. Averaging the absolute bound (5/32)^(number of insulated blocks), each eligible with independent probability 1/16, gives the stated decay by character orthogonality. Hence about 2/3 of *markings*, not necessarily isomorphism classes or distinct counts, hit the lowest possible valuation r−2.

## Chan–Pak implication and precise scope

A literal Collatz-like F(H)=3H+2 sends H+1 to 3(H+1); every positive iterate is 2 mod 3. Once r≥3, no marking of this one-apex triangle train can realize such an iterate. The carry that drops the first two factors of three is exactly the nearest-neighbor product sum, not the product of separate local count residues. This says which selective contacts spend the available two factors, but no choice of those contacts can spend a third. The additive induction barrier in [TRIANGLE_TARIFF](TRIANGLE_TARIFF.md) applies to a different, broader induced-extension question and does not require the one-apex assumption.

[Replay](replay_triangle_clock.py) checks the symbolic trits against direct subset path counting. No claim that either obstruction holds for unrestricted replacement of the base.
