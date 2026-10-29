# Session91 independent audit: operations, metrics, and retained information

Audit date: 2026-10-29. Source pin: `90e704a8f35614d2085778d1b518f03ff0d750ca`.

**Verdict:** the reviewed all-order arguments are sound within the explicit scopes below. Ten independent audit programs passed **3,253,649 active checks**, with byte-identical JSON results under ordinary Python and `python3 -O`. No external-priority claim, physical packing optimum, major conjecture closure, or new unrestricted tournament enumeration follows from these checks. Finite tests support, but do not replace, the proofs.

## Reproducible evidence

Run each `independent_*.py` with ordinary Python and `python3 -O`. Checks raise exceptions explicitly, rather than relying on removable `assert` statements. Recorded results are the corresponding `*.normal.json` and `*.optimized.json` files.

| Program | Active checks | Independent scope |
|---|---:|---|
| `independent_golden_quotients.py` | 46,310 | Every rooted cyclic legal word through period20; exact integer quotient keys; 70 period lifts |
| `independent_residue_automata.py` | 2,621 | All2,025 observer/target cases modulo1..9; independent complete-DFA minimization; four11/121 controls |
| `independent_guard_leakage.py` | 107,688 | Every graph through5 vertices,120 deterministic random graphs on6..9 vertices, cycles3..10; exact integer kernels and modular lower bounds |
| `independent_metric_carriers.py` | 6,817 | Golden and hexagonal matrices; toroidalK7 distances; affine group21;14 triangular faces and two Hamming codes |
| `independent_golden_dynamics.py` | 117 | Every fixed named toggle through period18; even zero-fibre transition probabilities |
| `independent_geometric_padding.py` | 3,041,519 | 4,656 padding trees across beta1..96; independent AHU rigidity/nonisomorphism and exact distances |
| `independent_lrc_geometry.py` | 7,904 | 162 positive-speed sets; full rational sawtooth arrangement; good-set intervals, reflection images and AP13 hostile |
| `independent_lucas_operator.py` | 81 | LiteralW5 forbidden-neighbor matrix; exact Smith minors, modular ranks, depth3, weighted gauge |
| `independent_late_repairs.py` | 40,583 | Primitive prime-power norm controls; isolated alternating rows for even cycles4..14 and exact source recovery |

The dynamics program imports quotient primitives from the audit's own quotient program; the Lucas program imports the audit's own modular elimination. The padding audit deliberately imports only the producer's point-set constructor, then independently reconstructs all distances, graph edges and rooted-tree types. Other mathematical implementations do not import producer algorithms. Standard-library exact arithmetic is used; no floating eigensolver supplies proof evidence.

The additional `independent_integral_frame.py` supplies9 exact checks of the integral five-cycle frame, its determinants, and the two characteristic-five rank profiles. Both execution modes agree.

## 1. Correct inherited golden quotient

The reviewed Session74 correction of THM4592 is essential. For O=Z[phi], phi²=phi+1, and R_n=O/(phi^n−1), the cyclic legal-word reader is bijective for oddn, but for evenn has exactly one three-word fibre: empty and the two alternating words. Its order is L_n−1−(−1)^n, not uniformly L_n. The inherited seven-carry proof supplies the all-order claim; our integer-pair quotient implementation exhausts n<=20.

Period lifting m to rm multiplies by 1+phi^m+...+phi^((r−1)m). It is injective, but reduction back multiplies the source by r. Source preservation is not the same as a section of reduction. Rational2-adic coding does not prove positive-integer Collatz behavior.

For evenn>=4 any *one fixed named* totalized toggle, together with the old reader output, separates the zero triple. Empty becomes a singleton; one alternating word stays put; the other loses an occupied bit and has a distinct nonzero image. All other fibres were already singletons. Thus the minimal future carrier for these specified observations is the full word set. Uniform attempted toggles are not strongly lumpable: probabilities into the zero fibre are0 and1/2. At n=2 quotient-only observations remain constant; legality as a separate output changes the problem.

## 2. Observer norm and compatible joint information

The two-register high-to-low Fibonacci update z -> Az+d(1,2), A=[[0,1],[1,1]], with one previous-digit bit, is correct. For observer ell=(a,b), the future-invisible subgroup is ker[[a,b],[b,a+b]] modulo m. Cayley–Hamilton makes it invariant under all future shifts. Legal fixed-length blocks produce translations(1,1) and(0,1), proving reachability and distinguishability, including legality-bit and dead-state distinctions.

For a reachable target, the minimum complete DFA has2J+1 states, where J is the index of that kernel. For nonzero(a,b), g=gcd(a,b), Delta=a²+ab−b², Smith form gives

    J = m² / (gcd(m,g) gcd(m,abs(Delta)/g)).

Unreachable targets give the one-state empty language. The zero observer with target0 gives three states. The classical ordinary-value baseline is credited in the producer note, not claimed new.

The exact four control counts are243,23,2663,243 for (m,a,b)=(11,1,0),(11,1,4),(121,1,4),(121,1,37). A root must be lifted with the modulus: keeping4 loses precision at121. For p-primitive observers, inert primes force J=p^(2k), split roots permit J=p^k, and ramified5 has minimum5^(2k−1). The last statement follows from a=2b+5c and Delta=5(b²+5bc+5c²), with b a unit.

Joint ordinary moduli m,n share both register compatibility modulo gcd(m,n) and the same legality bit. The joint live carrier is2 lcm(m,n)², not the independent product. General simultaneous linear targets require membership in the *joint* image, not merely separate reachability. The intersection-kernel formula follows by the same translation argument. This supplies genuine information recovery from explicit compatibility equations, not a license to reconstruct arbitrary tournament ports from marginals.

## 3. Exact operator loss, including hostile characteristics

For independent-set feasibility of any graph H, L=PB(1−P) and D=LL^t=PB²P−(PBP)². Forbidden-neighbor columns have either two legal predecessors (one conflicting edge) or one (a conflicting star with at least two edges). Thus D is a signless exchange Laplacian plus pins.

In a connected graph every independent set with at least two tokens reaches a pin by sliding along a closest-pair geodesic; any premature collision would contradict minimal distance. Only the empty state and, for bipartite H, the alternating singleton vector survive. Tensor-product decomposition yields nullity2^b for b bipartite connected components, including isolates. The theorem is over Q/R. Our literal incidence matrices, explicit integer kernel basis, and modular rank lower bounds independently certify rational ranks in the tested universe.

For W5, D=[0] direct-sum (2I+C5) direct-sum (3I+C5). Exact minor gcds give Smith diagonals(1,1,1,1,4) and(1,1,1,5,25). Modulo5, L still has rank10 but D has rank8: Gram rank need not equal factor rank in finite characteristic. The pair block has nonzero square and zero cube. This is not loss of invertibility of the odd-characteristic Walsh change of coordinates.

The15 legal insertion arcs embed in Paley11; they are not its55 arcs, and the full plus/minus Cayley graph is K11 rather than the guarded graph. Nonzero real weights permit a diagonal gauge preserving rank, not the integral lattice or arbitrary zero-weight behavior.

For evenn>=4 the reader pushforward and D have equal rank but different kernels with zero intersection. The alternating-word rows of D are isolated diagonals n/2. Reader output plus those two scalar observations recovers all linear amplitudes; two extra linear measurements are minimal. This statement is distinct from the discrete three-way tag, and division by n/2 cannot be silently used modulo a prime.

## 4. Arithmetic clocks do not supply physical isometries

The global11-square optimum is **source-reported**, not reproduced by these tests. The primary `11SquaresFormalized` repository describes Lean-kernel reasoning plus native numerical certificates. The older golden contact-family dictionary was inherited as an accepted supplied premise; its original PDF was unavailable for independent rereading. It concerns a restricted marked family, not every physical packing. Octic/S8 claims are also inherited, not recomputed here.

The golden quotient's order-five cell action does not extend to a metric rotation. Q^5−I=Q³(4I−Q), while Q^5 is not I. No nonzero symmetric form satisfies Q^tGQ=G; the actual quadratic norm is split and anti-invariant. Quotient labels erase integer displacement.

The affine Paley group Cp:C_((p−1)/2), p=3 mod4, embeds in flat2-torus isometries exactly for p=3 or7. A finite odd linear image has order at most3 by the integral-trace restriction; an abelian normal translation kernel in the nonabelian affine group lies in Cp. Hence (p−1)/2<=3. The explicit hexagonal construction realizes21; translation realizes3.

The seven points x(3,1)/7 in metric[[2,−1],[−1,2]] have pairwise squared geodesic distance2/7. Their numerator norms are positive multiples of7 and explicit representatives attain7. This is toroidalK7, not planarK7 and not a Hadwiger–Nelson bound. Chosen shortest displacement lifts have nonzero cycle sums: the missing lift data prevents integration into one planar configuration. The14 oriented faces triangulate a torus and split into two Fano systems; their binary spans are Hamming codes intersecting in{0,all-one}. These classical structures are not claimed newly discovered.

## 5. Euler repair preserves the right symmetry, not every symmetry

Attaching beta distinct rigid trees to a separate connected core of cycle rank beta changes Euler data to e=v−1 while preserving abstract graph automorphisms. The spine/tooth construction has a unique diameter and forbids reflection; the audit independently confirms4,656 instances with exact AHU types. Separated components preserve faithful unit-distance realization and, when present, minimum separation1.

The displayed point-set symmetry does **not** remain the core's ambient Euclidean symmetry. Every abstract automorphism fixes the distinct rigid padding components pointwise; one has three noncollinear points, forcing any Euclidean isometry of the whole displayed set to be identity. A C4 core may retain abstract D8 but lose all ambient isometries. No packing density or enclosing-region guarantee survives this padding operation.

For square tiles distinguish positive-length side contact from any boundary contact: the14-vertex hostile has13 side contacts but15 boundary contacts. The asymptotic padding lower bound applies only to the stated disjoint core-preserving odd-automorphism contract.

## 6. Reflection persistence and rational separation

For positive integer speeds including an even one, positive superlevel good sets avoid0 and1/2; reflection swaps the two halves. The entire H0 persistence module, not merely each dimension, is a regular C2 factor times the half-module. For all-odd speeds the fixed contribution is different. Unit period and any gcd normalization must remain explicit.

The classical opposite-slope/pair-sum maximum-clock lemma gives denominator at most2M, M=max speed. A rational optimum lambda cannot equal2−phi. The nonzero integer numerator of lambda²−3lambda+1 gives the strict gap >1/(12M²). Our independent full sawtooth arrangement verifies this without relying on the producer's candidate-clock search. The AP13 and mixed-parity singleton/empty controls pass. There is no height-independent gap and no lonely-runner conjecture proof.

## Synthesis supported by the audit

The common mechanism is **operation-relative information loss**. A quotient can preserve a clock while losing metric displacement; preserve outputs while losing future legal moves; preserve rank while having the wrong kernel; or preserve abstract automorphisms while destroying ambient geometry. Positive recovery results require explicit compatibility equations, observable kernels, or additional measured rows. These are structural explanations rather than numerical analogies between11,7,15 or21.

No mathematical correction remains outstanding in the audited lane notes. Original inaccessible packing premises and unreproduced external formalization remain clearly marked evidence boundaries, not silently promoted to independently verified results.

## Late integral-frame bridge

The displayed P in `zeckendorf/GOLDEN_FIVE_CYCLE_FRAME.md` has augmentation minor−1 and intertwines the five-cycle block with two copies of A² over the integers. The ambient constant-plus-augmentation frame has determinant−5. Thus the index5 obstruction belongs to splitting off the constant mode, not to splitting the two golden blocks. Modulo5 the constant vector lies inside augmentation; the ambient shifted operator has Jordan blocks3,2 while its augmentation restriction has blocks2,2. Independent matrix multiplication, determinant expansion, and modular elimination verify all displayed identities. This is an operator-module isomorphism, not an isometry or an identification of legal binary configurations.
