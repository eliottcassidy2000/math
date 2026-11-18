# Midpoint-out certificates, the fifteen-prime chambers, and all mixed local levels

The goal is to make “prove 7 and 21 absent last” precise without smuggling in Chan–Pak. Three results are separated below: a literal shared host for the midpoint interval, a symbolic infinite family in which only two leaves need individual checks, and an explicit local selector at every mixed (3,5) precision. None asserts full absolute coverage.

## 1. Middle out on [7,21]: one host, not one reflection

Take the strong four-vertex core with outgoing masks (6,12,8,1), H=5, and append two ordinal singleton sinks s→t. Its end counts E_x are (1,1,2,1). For each marked core vertex x, let R_x be the number of ordered covers by an optional prefix and a nonempty suffix starting at x; the four reads are (7,9,6,5). The replay obtains both vectors directly from subset DP.

Reverse x→s. The old paths not using x→s contribute 5−E_x. A path using the new s→x reads one of the rooted covers, with t still last. Thus the four counts are 5−E_x+R_x=(11,13,9,9). Alternatively reverse x→t. The old paths survive (the join of the core first touches s), and the first sink s can be allocated either before t or after the second core run. The result is 5+2R_x=(19,23,17,15). With k≥1 earlier sinks, the same partition of the sink list gives 5+2^k R_x. Neither 7 nor 21 has that form: 7−5=2 and 21−5=16, while none of 7,9,6,5 is a power of two.

For the requested midpoint 14, the interior is certified by the common host, in these pairs:

| Distance from 14 | Counts | Least orders | Reversed contacts, from the common host |
|---:|---|---|---|
| 1 | 13,15 | 5,5 | 1→s, 3→t |
| 3 | 11,17 | 5,6 | 0→s, 2→t |
| 5 | 9,19 | 5,6 | 2→s, 0→t |
| 7 | 7,21 | absent, absent | separate global exclusion |

The literal order-six host need not be the least-order witness for each value. The replay enumerates every labeled tournament through six to certify the displayed minimum orders. To make the last row global, use [Busch's strong minimum](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3): a strong witness for 7 would have order at most 4, and one for 21 at most 6. For a nonstrong 21 witness, stripping singleton SCCs leaves either that single 21-core or a factorization 3·7; the 7 case is already excluded. Thus the two endpoint proofs really are logically closed after the finite census.

**Reflection obstruction (all orders).** The pair 11,17 cannot be obtained from one tournament by reversing just one arc. In general when both counts are distinct positive multiplicative atoms of the physical spectrum (in particular, primes), a one-arc adjacency is witnessed after dropping common outer singleton SCCs at an order at most the larger of their two strong-order cutoffs. To see this, each side has one nonsingleton SCC. If the reversed arc is internal to an old component, no outside component can merge with it, and all the outside singleton components are irrelevant on both sides. If the old arc goes between components i<j, reversing it either just swaps two consecutive singleton components, leaving the count unchanged, or merges precisely the whole component interval i…j. The latter follows by finding a forward return route through an intermediate block or within a nonsingleton endpoint block and using the backward arc; all connections outside the interval still point the same way. The merged interval must contain the original nontrivial component, since otherwise its old count survives as a factor of the new distinct atom. Strip the unchanged outer singletons. For 11,17 the larger cutoff is six, and the exhaustive six-vertex edit check gives distance exactly two. The other mirror pairs 13,15 and 9,19 admit a single flip at order six. Thus counting symmetry about 14 does not lift to a uniform physical reflection. In the table, the middle pair really needs two separate branches.

## 2. What “the two holes last” means

There is no ordinary serial schedule in which finitely many operations occur at each time and an actual **finite** final proof is emitted after infinitely many distinct witness outputs. Nor can there be a fixed finite worst-case time for the two holes that exceeds witness emission time for every unbounded input, in a standard finite-symbol model. The [prior selector](../chanpak_prime_frontier_selector_20261118/README.md) terminates pointwise, but it supplies no such bounded guarantee.

There are two useful interpretations that avoid this trap.

1. **Finite windows.** For each radius r, decide all positive odd H in the window |H−14|≤r, using the finite least-order test, sort witnesses from the center, and emit 7,21 (if in the window) only after all other tickets have been checked. A third hole, if found, is reported as a counterexample, not buried by the schedule. The first nontrivial such window is the table above.
2. **Structural lastness.** Order positive odds by their least witness order, using distance from 14 to break ties and placing absence at infinite order. The witness part has order type omega: finitely many values are born at any finite order, and there are infinitely many distinct realizable counts by ordinal products of triangles. The full order has type omega+2 exactly when Chan–Pak is true. This equivalence describes a completed record, not a finite-time test for when one has reached stage omega.

More useful still is a finite symbolic certificate for **an infinite chamber**, where the endpoint checks really can be postponed until after a finite factor proof certifies every other member.

## 3. Ogg-prime chambers and two residual leaves

[Ogg's fifteen](https://celebratio.org/Ogg_AP/article/1157/) are primes p for which every supersingular j-invariant in characteristic p is rational over F_p, equivalently X₀(p)/w_p has genus zero. They are 2,3,5,7,11,13,17,19,23,29,31,41,47,59,71. This is not the list of supersingular primes of a single elliptic curve.

Let P consist of the odd Ogg primes other than 7. The [earlier full first-birth atlas](../chanpak_birth_atlas_metagraph_20261117/metagraph/PROOF.md) supplies literal witnesses for all p in P and all 7p with p in P other than 3, plus 49,63,343. These 29 graphs including the identity seed are included; every one has been independently recounted, largest order eight. The argument works for any finite bank meeting the same conditions, so Ogg membership is an organizing choice, not a prerequisite.

Write uniquely

    N = m · 7^b · product_(p in P) p^a_p ,   gcd(m, product of Ogg primes)=1,

for positive odd N. Fix any m with a witness. **Every point in this infinite chamber is covered by ordinal multiplication except possibly 7m and 21m.** If b is even build the 7 factors with 49; if b≥3 is odd use one 343 and then 49s. If b=1 and a_p>0 for some p≠3, trade one p for the seed 7p. If only 3s occur, use the seed 63 when a_3≥2. The exceptions are exactly a_3=0 or 1, giving the endpoint pair on the scaled segment [7m,21m] with midpoint 14m. For m=1 both endpoints really are absent. For general attained m the possible witness-status states for (7m,21m) are (0,0),(0,1),(1,1); (1,0) is impossible because multiplying by 3 preserves realization. The distinction at 1 and 3 is therefore transparent: they are the first two levels of a single 3-adic ray inside the 7-once sheet. This is an exact three-state implication structure, not a forced tournament isomorphism.

With x_p marking exponent a_p, and y marking b, the exact count-set series on a witnessed chamber is

    1 / [(1−y) product_(p in P)(1−x_p)]
       − [7m absent] y − [21m absent] x_3 y.

The exception tickets can be emitted last **after a finite symbolic proof of the rest**, even though infinitely many individual counts were not printed. For arbitrary m not yet witnessed, split b into the two parities. Each half is an upward-closed subset of an exponent lattice under multiplication by p∈P and 49, so Dickson's lemma gives finitely many Pareto-minimal corners per half and a rational multivariate series by inclusion–exclusion of orthants. This is a structural finiteness statement; merely enumerating witnesses does not announce when one has found all Pareto corners. The existing pointwise decision algorithm handles each ticket but cannot yet prove that every new chamber has the required base and 7m witnesses.

The exact compiler route(n) returns these factors or the outstanding exceptional leaf, and can materialize a witness with ordinal joins if the base is supplied. The allowed stripped-base residue density among odd inputs is the exact product prod_(odd Ogg p)(1−1/p)=297271296000/1050750592307≈0.282913. This is a density of **input cores**, not a density of missing counts or a proof that the remaining bases are possible. The earlier prime/7p atom reduction remains sharper if the only goal is an all-order proof.

## 4. Widen every mixed modular level using a target-dependent address

The [three-probe wheel](../chanpak_fifteen_probe_wheel_20261118/README.md) was independently finite-verified through modulus 810. There is now an all-level **one-shot** construction at these primes; it is not a claim about all-level strongly connected components. Let Q be the same four-vertex strong probe at mark 0. For a seed with h=H and B=f(2), add an ordinal transitive block of order t, then apply Q. The output is strong and has

    H_out = h + 2 B · 2^t.

Indeed ordinal padding keeps h and multiplies B by 2^t, and Q reads h+2B. The uniform substitution is strong because one can travel through the strong quotient to move between modules, or return to the module to reach a different vertex in it.

Four six-vertex seeds with (h,B)=(17,224),(23,296),(37,376),(43,424), and t=j+4a for j=1,3, give respectively mod-30 pairs (13,1),(7,19),(11,23),(29,17). The coefficient 2B·2^j is a unit at 3 and 5 for each. Given a target unit modulo 2·3^k·5^ell, select its base chart, subtract h and divide by that coefficient in the odd modulus. The result is 1 modulo 15. Since 16=1+15, binomial lifting gives v_p(16^u−1)=1+v_p(u) at p=3,5. Thus 16 generates the whole principal-unit subgroup at each local modulus. Their exact orders, 3^(k−1) and 5^(ell−1), are coprime, so CRT chooses a single a satisfying both. mixed.py implements digit-by-digit lifting and prints the constructive address.

This proves that every odd unit at **every mixed precision** has a physical representative by just four seeds. But for fixed seeds the eight resulting rays are still exponentially sparse (only O(log X) distinct values ≤X). A strong all-level congruence theorem is therefore not the missing exact placement theorem. The template also contains no absolute 7 or 21, even though it realizes their applicable local shadows.

## 5. Second-hop archive tests of the analogy

* The [Ogg / triangular / Chowla packet](../ogg_triangular_chowla_20260927.md) distinguishes the true genus-defect cutoff from a residue filter. In prime level p>3, Riemann–Hurwitz gives g⁺=(2g+2−f)/4 with f the Fricke fixed points. Equivalently g⁺ counts half the nonrational supersingular j-values. The genus grows linearly while the CM fixed-point contribution is O(sqrt(p) log p), so a finite cutoff is a theorem. The exact shared tournament object is the **marked Paley half-row**: for p≡3 mod4, p>3, its residue/miss imbalance is (2−(2/p))h(−p). This is not H(Paley). For Chan–Pak a corresponding universal height/rebate inequality that *forces* the existence of every remaining witness is still missing.
* Following the [Monster-level square-sum study](../monster_levels_square_sum_7_21_20261022.md) leads to Pell arrival clocks. The newest vertex sees floor(sqrt(2n−1))−floor(sqrt(n)) old neighbors. On the special Pell rays an entering square and a departing square leave the degree unchanged but exchange actual ports; the first nondegenerate n=9 is still disconnected. This is a concrete warning that equal numeric “middle” signals do not mean the path interface survived. The n=15 square-sum onset instead has a compulsory ear.
* The [Lucas orbit study](../reciprocal_lucas_carriers_20261028/THEORY.md) provides an exact positive-control theorem for coordinating local clocks. Its half-turn criterion says that for odd p≠5, divisibility p|L_n requires n≡d_p/2 mod d_p with even projective order d_p. Multiple primes share an index iff all their orders are even and have the same 2-adic valuation: necessity is v₂(n)=v₂(d_p)−1, and sufficiency takes half the least common multiple (all quotients odd). The solutions are precisely that half-period class modulo the least common multiple. For instance p=3 and p=11 have different valuations, so no Lucas number is divisible by 33 despite separate realizations. Our 3/5 selector succeeds locally for an explicit reason (coprime lift orders after choosing the chart), not because independent local works can generally be combined.
* The [fifteen-projective atlas](../relation_gluing_7_15_20261028/ATLAS.md) has a positive control: Hamming codewords really do become the 15 points of PG(3,2), and a chosen symplectic form picks 15 of the 35 lines. The expression “fifteen” alone does not identify that incidence structure with Ogg primes or chamber coordinates. The constructive content is the chosen observer and the reconstruction certificate.

The finite window search, the chamber factor proof, and the modular selector overlap but answer different questions. A last two-hole conclusion over *all* positive odds still needs an absolute theorem for bases/primes outside the seed bank.

Reproduce with `python replay.py | diff -u replay.out -`, `python mixed.py | diff -u mixed.out -`, and `g++ -O3 -std=c++17 edits.cpp -o /tmp/mid-edits; /tmp/mid-edits | diff -u edits.out -`. The midpoint host and the all-level four-seed lift came from a background research lane and were independently recounted here.
