# Deficit-six provenance, refined channels, and the almost-regular strongness gate

Session30, 2026-10-22. Status: PROVED algebra and structural lemma; inherited classifications cited with exact scope; claimed six-order closure NOT RECOVERED in the current repository. No new census was performed.

## 1. Recovered mathematical type and evidence boundary

For a strong self-converse tournament T let h=H(T)/|Aut(T)| and let b be its number of blue Hamilton-word orbits. A word is blue precisely when positional reversal is an anti-automorphism. The deficit is Δ=h−b, always even. R6 counts strong self-converse isomorphism types with Δ=6 by number of vertices. It is not a score variance, feedback deficit, labelled path count, or adjacency polynomial.

Sources read:
- `05-knowledge/results/deficit_six_centers_20261021b/README.md`, §§1,5,6: exact definitions; three small centers with (order,h,b)=(5,9,3),(6,15,9),(7,9,3); original cutoff and trivalent SCC formula.
- `05-knowledge/results/deficit_six_centers_20261021b/COMBINED_CUTOFF.md`: combined necessary support, complete order-ten input, minimal recurrence order18.
- `05-knowledge/results/linked_tail_geometry_and_six_center_orders_20261021c.md`, §3: exact remaining support15,18,20,21,24,28, degree bound28, safe onset41.
- `05-knowledge/results/center_transitivity_20261021c/README.md`, `ORDER35.md`, `VT15.md`: full-orbit versus regularity distinction; exclusions25,35; vertex-transitive order15 exclusion only.
- `05-knowledge/results/deficit_twelve_20261021c/README.md` and `deficit_fourteen_20261021c/README.md`: complete exclusions12,14 and audit scope.
- `05-knowledge/results/deficiency_three_20261020/centers.md`, §§2–3: exact SCC divisor restriction.
- `05-knowledge/results/fixed_degree_metagraph_rationality_20261021/README.md`: inherited fixed-degree framework; polynomial unknowns not wiring.

The four preceding complete universes contain262144,2048,8192,131072 frames at orders12,14,25,35 respectively. Frames can repeat isomorphism types. Order12 records independent proof/API review but not an independent full bank replay; its normal/optimized bank runs are producer runs. Order14 records independent root review and optimized replay. Order35 records independent proof/API review and optimized complete-bank replay. None of these audits excludes the six residual orders.

The user reports that all six residual orders have now passed independent audits. The publisher's targeted current-tree search did not locate that subsequent packet. Therefore this integration does not represent the six closures as locally source-verified. It retains

R6 = z5+z6+z7+C(z), C=c15 z15+c18 z18+c20 z20+c21 z21+c24 z24+c28 z28,

with nonnegative unknown integer coefficients in the recovered sources. Claims below that set C=0 are explicitly conditional on the missing closure. This is a provenance gap, not a mathematical refutation of the reported result.

## 2. Exact known-family refinement

Put p(z)=1−z−z3 and u(z)=p(z2). The unique self-converse SCC decomposition is A→C→A^op, with empty center allowed. If a=h(A), c=h(C), d=b(C), then h=a²c, b=ad, and Δ=a(ac−d). At Δ6 the only cases are a1 with a deficit-six center, or a3 with c=d1. The proved three-presentation half series is H3(z)=z5/p(z)².

It follows that the families built from the three known centers and the nonrigid half have exactly two signatures:

F_(9,3)(z)=(z5+z7)/u + z10(1+z+z3)/u²,
F_(15,9)(z)=z6/u.

These are unconditional exact series for those specified families. They exhaust all self-converse deficit-six types if and only if C=0. Each additional unknown strong center contributes its own signature with rigid palindromic padding; no unverified signature restriction is imposed on it.

Within the known-family series, (15,9) occurs only at even orders. Setting t=z² and q=1−t−t³, the even-order channels are t3/q and t5/q². The odd-order (9,3) channel after extracting a factor z is t2(1+t)/q+t5(1+t)/q². These are exact formal identities, not asymptotic approximations. If C=0 they are complete parity-refined global counts.

## 3. Onset19 is sharp, but recurrence order is18

The degree-three merged-node series is

G=z5/2(1/p²−1/u)+R6/u+z10(1+z+z3)/u².

Here degree-three means fundamental black degree3 of the counted vertices. It does NOT mean a recurrence of order3.

For C=0 the exactly reduced denominator is D=p²u², degree18. Its numerator is

N=z5−2z7−2z8+z9+3z10−z14+2z15−z16+z17−z18.

The checker independently obtains gcd(N,D)=1, degree(N)=degree(D)=18. The homogeneous recurrence from D therefore holds from index19 and fails at18, where its residual is −1. Thus19 is the sharp onset for this canonical denominator recurrence, not just an upper bound. Minimal eventual recurrence order remains18.

The generating-function subtraction isolates precisely the outstanding source:

u(z)[G(z)−G0(z)] = C(z),

where G0 uses the three known centers. Thus the six residual coefficients are the six potentially nonzero coefficients of a finite innovation polynomial. This gives a source-sensitive audit identity: independent true merged-node counts can recover c_j by g_j−g_(j−2)−g_(j−6) minus the same expression for G0. It is NOT a substitute for an actual count or exclusion certificate. Missing center contributions cannot be hidden by recurrence fitting.

## 4. A proved bridge to almost-regular tournaments

Lemma. Every almost-regular tournament of even order n≥4, with every score in {n/2−1,n/2}, is strong.

Proof. If not, choose the source SCC of order a and put b=n−a≥1. Every edge from it to the complement points outward. Its average full outdegree is b+(a−1)/2=(n+b−1)/2. Since every score is at most n/2, b≤1. Consequently b=1, and the sole complementary vertex is a sink of score0, contradicting the lower score bound n/2−1≥1. QED.

This is an actual map-compatible restriction: passing from a tournament to its ordered SCC word, the almost-regular condition forces word length1. It removes every nontrivial rigid or nonrigid padding branch above. It does not identify score variance with black deficit or the forest coefficients.

Among the three known deficit-six centers, only the order6 source C3[K2] has even order. Its literal rows (36,17,10,19,37,14) have scores (2,2,2,3,3,3) in sorted order, and exact subset DP gives H45. Inherited Aut3 and b9 give h15 and Δ6. Therefore:

- Unconditionally, the intersection of the known deficit-six families with almost-regular even tournaments is exactly C3[K2].
- Conditional on C=0, it is the unique self-converse almost-regular deficit-six tournament, at order6; no such source exists at any even n≥8.
- Without the missing closure, the possible additional intersections remain confined to strong centers at18,20,24,28. This excludes15,21 by parity, not by a purported new census.

The no-padding lemma is relevant to interpreting almost-regular loop-expansion experiments: SCC growth mechanisms responsible for the rational metagraph count series are absent in that sampling family. It supplies no asymptotic formula for kappa/rho and no complete wiring classification.

## 5. Reproduction

Python3 with SymPy:

python -B check.py
python -B -O check.py

Both runs agree byte-for-byte with receipt.json. Tests use explicit exceptions and remain active under optimization. They check rational identities, numerator/denominator gcd, residual isolation, known coefficients5..14, nonnegative coefficients through80 for G0, and the literal n6 score/H control. They do not inspect or certify absent six-order proof banks.
