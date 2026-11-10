# Value-refined cycle-index engines and a sharp strong seven-symmetry gap

Session239, 2026-11-10. **PROVED from the explicitly inherited theorems below, plus reproducible FINITE-EXACT inputs.** The recurrence and sharp sevenfold gap were independently audited in `../audit/README.md` and `../audit/strong_gap_audit.md`. No literature-priority claim. Chan–Pak completeness remains OPEN.

## 0. Inheritance, types, and research board

H(T) counts ordered directed spanning Hamiltonian paths; H(empty)=1. A000568 counts ordinary tournament isomorphism classes, not H values, H maxima, or converse-merged classes.

Closest mechanisms:

- THM-4803, `01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md`: the scalar SCC identity A=1/(1-S) is inherited, not new.
- THM-3380, `01-canon/theorems/THM-3380-hamiltonian-deletion-layer-monoid-semiring-and-small-order-boundaries.md`: the multiplicative monoid algebra [a][b]=[ab] and H multiplicativity under ordinal sum are inherited.
- THM-5094, `01-canon/theorems/THM-5094-cyclic-orbit-profile-cost-and-sharp-sevenfold-path-floor.md`: exact cyclic profile floor, its ordinal-sum equality characterization, and mu_3=3, mu_5=15, mu_7=175. The referenced Session238 source directory was absent in this checkout; the canonical theorem statement was present and used explicitly as an inherited dependency.
- LEM-003, `01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md`: the automorphism group acts freely on Hamiltonian paths. The directed-cycle analogue is false, already witnessed by C3.
- THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`: the full path-cover substitution kernel. Its two-path specialization below is inherited; the new step is using it after a prime-action/module reduction.
- THM-4987, `01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md`: inherited module-interaction costs 27 for R5 and105 for P7 after division by their full groups. These equal raw costs135 and2205 below. Its explicit warning says those are not bounds for arbitrary seven-symmetric tournaments. The orbit-profile split here supplies the missing action-specific reduction, also covering Q7.
- THM-4744, `01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`, and `05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`: fixed-H rational series and eventual polynomiality were already proved, including stronger arithmetic degree information. They are not rediscovered as new results here.

Canonical hostiles: C3 followed by C3 has H=9, not6; C3 followed by a singleton and the reverse SCC word are distinct classes; a transposition fixes no tournament; Q7 has H/7=25 whereas P7 has H/21=9, so a subgroup denominator cannot be replaced by the full-group denominator. Corrected near miss: THM-4803's old cyclic-triple/automorphism invariance is false and is not used.

Board: (i) scalar Burnside fast counts; (ii) H-valued SCC algebra; (iii) cyclic orbit-profile floors; (iv) two-path-cover inflation cost; (v) fixed-target divisor quotient. Anchor is a spectrum-preserving bridge from enumeration; niche is strong class generation; wildcard is the action-gap theorem in section4. The least-used relevant sidecar is the ordered two-path-cover count, which turns out to be essential.

## 1. Fixed-permutation SCC recurrence retaining H

Let m=(m_1,m_2,...) be a finite cycle multiplicity vector, |m|=sum d m_d, and let sigma_m be any labelled permutation of that type. Work in the integral multiplicative monoid algebra with basis [h], h>=1, and [a][b]=[ab]. Define

    F_m = sum_{sigma_m T=T} [H(T)],
    G_m = sum_{sigma_m T=T, T strong} [H(T)],
    F_0=[1], G_0=0.

Strong includes a singleton. Put binom(m,k)=product_d binom(m_d,k_d).

**Theorem1 (PROVED).** For every nonzero m,

    F_m = sum_{0<k<=m} binom(m,k) G_k F_(m-k),
    G_m = F_m - sum_{0<k<m} binom(m,k) G_k F_(m-k).       (1)

Here k<m means k<=m and k!=m, not coordinatewise strictness.

Proof. The condensation of a tournament is a uniquely ordered list of strong components. An automorphism fixes each position of this list, hence its first SCC is a union of whole sigma_m cycles. Choosing k_d of the m_d labelled d-cycles has the displayed binomial multiplicity. The first block and the remainder are independently invariant, and every cross arc is forced from the first block to the remainder. Conversely this construction uniquely recovers the fixed tournament. A Hamiltonian path crosses the block boundary once and H multiplies. This proves the identity coefficient by coefficient. A proper nonempty complement lowers vertex order, so the subtraction is triangular and uniquely defines G.

Write z_m=product d^(m_d) m_d!. Weighted Burnside gives

    U_n = sum_{|m|=n} F_m/z_m,
    V_n = sum_{|m|=n} G_m/z_m.                          (2)

These are integral H histograms of ordinary and strong ordinary classes. This follows by applying the ordinary orbit-counting proof separately to every invariant H fibre. Define their order series U(x), V(x). Then

    U(x)=1/(1-V(x))                                     (3)

in the formal power-series ring over the monoid algebra. All coefficients at a fixed order have finite support. No analytic convergence is asserted. More precisely, the cycle-variable series sum_m F_m y^m/z_m and sum_m G_m y^m/z_m satisfy the same reciprocal identity. Equation(1) is its integral unnormalized implementation.

Map / retained predicate / loss: fixed labelled tournaments -> class histogram by weighted Burnside preserves existence and multiplicity of every H value; evaluation [h]->1 loses the value and recovers only A000568. Restoring the H monoid coefficients is the necessary sidecar. No computation of the scalar coefficients alone recovers (2).

## 2. A finite exact target algebra and genuine pruning

Fix h>0, h=product p^e_p. The span of [a] with a not dividing h is an ideal: ab dividing h would imply a divides h. Thus projection onto divisors is a ring quotient R_h with

    [a][b]=[ab] if ab divides h, and 0 otherwise.

It has rank tau(h). Equivalently

    R_h = Z[y_p : p divides h]/(y_p^(e_p+1): p divides h),
    [a] <-> product y_p^v_p(a).

All other prime variables are set to zero. For h=1 the ring is Z. Equations(1)–(3) hold exactly in R_h, so discarding nondivisor H values at every step cannot damage the target answer. A dense multiplication needs precisely at most

    B(h)=product_p binom(e_p+2,2)

ordered divisor-pair products (a,b) with ab|h, not h^2 arithmetic. This does not bound the cost of producing the input fixed-tournament histograms. For a prime target the quotient is a dual-number ring; this isolates the one nontrivial strong factor, an inherited SCC consequence rather than a coverage proof.

For a permutation of odd cycle lengths d_1,...,d_t, put M=lcm(d_i), P=product_i mu_(d_i), with mu_1=1. LEM-003 and THM-5094 give

    coefficient [h] F_m = 0 unless M|h and P<=h.          (4)

An even cycle fixes no tournament, because its half-power reverses an antipodal pair. For fixed target h only cycle lengths dividing h can contribute, and each nontrivial orbit consumes its mu_d factor from the budget. The same conditions apply to G_m.

**Strong refinement.** If t>=2 and G_m has coefficient at h, then

    h >= P+2M.                                          (5)

Indeed equality in THM-5094 forces the orbit blocks to form an ordinal sum, which is not strong for two or more blocks. Both H and P are odd multiples of M: mu_d is divisible by d, and lcm(d_i) divides their product. The next admissible multiple after P is therefore P+2M. This refinement is a direct credited corollary, not a new proof of the inherited cyclic floor.

Cost boundary: the identity profile has binom(n,2) independent pair orbits and naive production still requires 2^binom(n,2) tournaments. Neither the divisor quotient nor weighted Burnside makes unrestricted spectrum extraction polynomial-time. The script is a correctness verifier, not a claimed large-n speed record.

## 3. Related scalar sequences, efficiently once A000568 is known

Let a_0=1 and a_n=A000568(n). Let s_n count nonempty strong ordinary classes, with s_0=0 and s_1=1. For n>=1 these are [A051337](https://oeis.org/A051337); OEIS assigns the extra convention A051337(0)=1, which is not included in the nonempty strong-component series. The inherited identity gives

    s_n=a_n-sum_(k=1)^(n-1) s_k a_(n-k).

Thus a prefix of strong class counts needs O(N^2) extra scalar arithmetic, or formal series inversion using faster polynomial multiplication when that is useful. No faster bit-complexity than the output cost is claimed. The prefix for n=1..10 is

    1,0,1,1,6,35,353,6008,178133,9355949.

For labelled tournaments L_n=2^binom(n,2), the analogous EGF recurrence is

    ell_n=L_n-sum_(k=1)^(n-1) binom(n,k) ell_k L_(n-k),

with prefix 1,0,2,24,544,22320,1677488,... ([A054946](https://oeis.org/A054946), offset1, already cited in THM-4746). Ordinary and labelled recurrences must not be interchanged. Marking SCC number gives U(x,u)=1/(1-u V(x)); the coefficient at u^r counts exactly r ordered SCCs, not multisets. The verifier prints both strong sequences through20 and independently verifies their values through6 by direct tournaments/classes.

Fixed-H rationality and multiplicative support closure follow by the same coefficient algebra, but were already established in the credited Session24 packet. In particular, these transformations do not prove that the only missing odd H values are7 and21.

## 4. A sharp all-order gap for strong seven-symmetry

Let p be an odd prime. Define mu_p as in THM-5094 and

    gamma_p = min_R H(C3[R,1,1]),

where R ranges over the 2^((p-1)/2) labelled p-circulants for a specified p-cycle. For a tournament R put F2(R)=sum_(empty!=X!=V) H(R[X])H(R[V\X]), the number of ordered two-path covers.

The inherited THM-3121 specialization, also proved directly below, is

    H(C3[R,1,1])=3H(R)+F2(R).                            (6)

**Theorem2 (PROVED).** If T is strong and has an automorphism of order p, then either |T|=p and T is a p-circulant, or

    H(T) >= min(gamma_p, mu_p^2+2p).                     (7)

Proof. The action has only singleton and p-element orbits. It has at least one p-orbit because its order is p. If it has at least two p-orbits, (5) gives H>=mu_p^s+2p>=mu_p^2+2p.

Otherwise it has a single p-orbit R and at least one fixed vertex. Every fixed vertex sees all of R in the same direction, so R is a module. Contract R to a vertex r. The quotient Q is strong: project paths between representatives; a directed walk yields reachability in the quotient. Every vertex of a strong tournament on at least three vertices lies in a triangle. For r, its nonempty in- and out-neighborhoods would otherwise have all cross edges pointing from the in-neighborhood to the out-neighborhood, preventing r from reaching its in-neighborhood. Thus Q has r->b->c->r, and T contains C3[R,1,1] induced.

Hamiltonian-path count cannot decrease when vertices are added to an induced tournament. To see this without appealing to a named theorem, insert one new vertex into each old path: prepend if it beats the first vertex, append if the last vertex beats it, and otherwise some consecutive pair crosses from beating it to being beaten by it. Choose a deterministic available insertion. Deleting the new vertex recovers the old path, so the map is injective. Repeated insertion yields H(T)>=H(C3[R,1,1])>=gamma_p. This proves (7).

For (6), let R->b->c->R be the quotient. A spanning path either uses R as one run, in one of the three cyclic orders, or as exactly two runs in the unique pattern R_1,b,c,R_2. More than two R-runs are impossible. The first case contributes3H(R); the second is in bijection with all ordered spanning two-path covers of R. This proves (6), with all label/ordering factors retained.

**Corollary3 (PROVED using exact small inputs).** A strong tournament with an automorphism of order7 either has order7 and H in{175,189}, or has H>=2163. The latter threshold is sharp, attained on9 vertices by C3[Q7,1,1]. Consequently every strong tournament with H<2163 other than the order7 cyclic cases has automorphism group order coprime to7 (Cauchy's theorem). In fact every strong realization of203 or217 is asymmetric. The inherited Busch floor in THM-4744 is225 at order11 and nondecreasing, so such a realization has at most10 vertices. Its automorphism order divides H by LEM-003: the new gap excludes its possible factor7, and an embedding into S_n for n<=10 excludes29 or31. This is a constraint on realizers, not a nonexistence claim. The supplementary audit checks this deduction. In fact those two strong fibres are entirely asymmetric: the inherited Busch floor f(11)=225 forces their order to be at most10, excluding automorphisms of order29 or31, while LEM-003 forces every automorphism-group prime to divide H. No existence assertion is needed for this implication.

The complete eight-circulant computation gives F2(R)=1638 in all cases. Six have H175, hence inflation2163; two have H189, hence inflation2205. Meanwhile mu_7^2+14=30639>2163. This verifies all inputs of(7), with the inherited Q7/P7 identification and without any upper bound on ambient order. At p=3 the bound is15, sharp at C3[C3,1,1]; at p=5 it is135, sharp at C3[R5,1,1]. The raw135 and2205 inflation costs recover the specifically scoped module facts from THM-4987. The new sevenfold statement includes the cheaper Q7 case and arbitrary orbit profiles.

Boundary / hostile: dropping strongness destroys the corollary immediately. Any number of singleton SCCs padded around Q7 retains H175 and sevenfold symmetry at arbitrarily large order. Replacing H by H/|Aut(T)| also changes the problem. The theorem bounds raw paths under an actual order7 subgroup.

This creates a concrete bridge to Chan–Pak work: weighted Burnside's strong target coefficients for any nontrivial7-cycle profile vanish below2163 unless the profile is exactly(7). It removes an infinite strong symmetry sector, but it does not construct asymmetric realizers or prove unrestricted coverage.

## 5. Reproduction and exact universes

Run from the repository root:

    python3 05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/verify.py
    python3 -O 05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/verify.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/prime_gap.py
    python3 -O 05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/prime_gap.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/c7_low_fibres.py

- `verify.py`: all33,867 nonempty labelled tournaments at orders1–6, plus the auxiliary empty object; every permutation cycle type (including reversing even-cycle hostiles); independent direct fixedness vs edge-orbit generation; fixed-permutation SCC recurrence; weighted Burnside vs an independent isomorphism-orbit BFS using adjacent transpositions; all319 specified target/profile truncation checks. All pass. Normal/optimized output agrees bytewise. Explicit runtime checks are not removed by `-O`.
- `prime_gap.py`: all2,4,8 circulant orientations at p=3,5,7. Exact subset-DP covers versus whole inflated-tournament Hamiltonian DP; permutation-enumerated core paths; independently enumerated all9! permutations of a sharp9-vertex witness. All pass. Outputs record masks and both cover/H values; normal/optimized output agrees bytewise.
- `c7_low_fibres.py`: an additional bounded hostile probe, all specified C7-fixed tournaments at orders7–10:8,16,64,512 orientations, with strong minima175,none,2163,4151. Its H<225 completeness uses the separately inherited Busch floor f(11)=225 from THM-4744. Corollary3 does not depend on this bounded-order argument.

Masks index pairs (i,j), i<j, in lexicographic order; a set bit means i->j. Cyclic profiles use consecutive labelled cycles. Source/output hashes are in `SHA256SUMS`.

## 6. Claims and remaining frontier

Direct new deliverables are the value-retaining fixed-permutation SCC recurrence, its finite divisor algebra implementation, and the action-gap argument yielding the sharp2163 threshold. Their elementary ingredients and inherited numerical facts are credited above; no external novelty search was completed. Related scalar sequences and fixed-H rationality are recovered, not promoted as new. The minimum gamma_p for larger primes, better multiple-orbit strong floors, and a genuinely faster identity-profile spectrum engine remain open work. Chan–Pak and LRC(14) remain OPEN.
