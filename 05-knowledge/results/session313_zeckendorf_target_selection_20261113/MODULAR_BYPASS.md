# A uniform all-depth bypass for the first modular guard

This complements the [physical marks and contact lower bound](PHYSICAL_MARKS.md).

## Bypassing selection for the first modular guard at every depth

There is also an actual all-depth construction of a 3-coprime extension with a bounded multiplicative overhead over the necessary contact budget. This does not choose a specific exact count.

First, the cyclic join has a simple modulo-three law. For a block i write a_i=H_i and b_i for its ordered two-path cover count with both slots nonempty, reduced modulo 3. A block's ordered nonempty r-path cover count is divisible by r!, so only r=1,2 can contribute. A cyclic run word visiting all three blocks has length between 3 and 6. Lengths 3 and 6 each give 3 starts for the same product, hence zero. Length 4 contributes the sum of products with one b, length 5 the sum with two b. The odd-cycle expansion of the full polynomial Q_i(t) has powers of parity |T_i|, so

    b_i=Q_i(2)−2H_i ≡ Q_i(−1)+H_i = ((−1)^|T_i|+1)a_i mod 3.

Thus b_i=0 for odd blocks and b_i=−a_i for even blocks. For e even blocks the resulting join H modulo 3 is [−e+binom(e,2)] times the product of the three H_i. This coefficient is 0 at e=0,3 and −1 at e=1,2. In particular, three 3-coprime child counts and mixed size parities give a 3-coprime parent, and no target-dependent selection of child counts is needed.

At depth two, use the 72-vertex tower, choosing one marked vertex per top child in its first B-leaf, with label 1 in B. Color these three marked leaves red-blue-blue in cyclic order. The physical two-probe extension gives H≡1 mod 3 by quotient walk enumeration; its order is even. The red-red-blue variant gives H≡2 mod 3. One more universal sink preserves H and gives an odd-order seed. The modulo-three computation enumerates all orders of the probes and all cyclic runs in each intervening core segment, allowing empty segments. Its local tables count unconstrained ordered runs and runs with the designated mark as exact start or end. A forbidden endpoint is treated by inclusion–exclusion. At most two unrooted ordered slots or three rooted slots need be retained, because permuting the unrestricted nonempty slots gives a factor 3! when there are at least three. Inside each 24-vertex module, repeat the same quotient decomposition and compute the local tables by subset DP on B. The [exact modulo-three replay](two_probe_mod3.py) implements these finite checks.

For d≥2 build above these seeds by cyclic joins. An even parent uses one even and two odd children; an odd parent uses two even and one odd. Both are 3-coprime by the law above. Each child contains the corresponding monochromatic core as an induced tournament, so the final extension contains the full original depth-d core. Starting at depth two with costs (E,O)=(2,3), the two constructed costs satisfy (E,O)↦(E+2O,2E+O). Their sum triples and their difference changes sign. Choosing the cheaper final parity gives

    at most (5·3^(d−2)−1)/2 exterior vertices for a 3-coprime extension.

Compare the necessary lower bound ceil(3^(d−1)/2). They agree at d=2, and the ratio of upper to lower tends to 5/3. Both grow linearly with the number of bottom triples. The construction solves a prescribed divisibility guard at every depth, not integer count coverage or any arbitrary residue.

Both nonzero residues are available at the same budget: at depth two use one of the two base colorings, and at higher depth change a single base sign, which multiplies the final residue by $-1$. Within the grammar that keeps the cyclic joins above depth-two tiles uniform and only attaches new vertices within individual tiles, the upper cost is exact. At depth two the packet bound forces at least two added vertices, and parity forces costs at least two and three for the even and odd choices. Above that level mixed child parities and three nonzero residues are necessary by the cyclic-join law. The stated recurrence is therefore a matching lower bound in that grammar. Improving the global lower-bound gap requires a nonuniform attachment across tile boundaries or a changed core.
