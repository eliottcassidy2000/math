# Target placement by signed port contacts

**Exact placement:** a changed tile puts the seam center at **1090**, with endpoints **1089 and 1091**. The upper endpoint is the previously prescribed test value (2^5+1)^2+2. **New composition:** multiple nontransitive tiles can retain the exact ±1 switch by a signed product rule. **Limit:** varying both exterior tapes on a fixed single tile does not give all targets. Chan–Pak remains open.

## 1. The signed product rule

Put arbitrary tournament tiles D_i in ordinal order, so every arc between distinct tiles points forward. Add two vertices x,a. Inside a tile, let I_x,I_a be the subsets which beat x,a, respectively. The two switch variants differ only on x→a (minus) versus a→x (plus). Unlike the old twin-contact case, a vertex inside a tile may beat one port and lose to the other.

For a tile, retain four endpoint vectors: the start and end vectors of its Hamiltonian paths (W1), and the end of the first path and start of the second path in its ordered spanning two-path covers (W2). Require that each of these four vectors have equal sums on I_x and I_a. Call this an **endpoint match**. A singleton with identical port contacts trivially matches.

Let h_i=H(D_i), and form M_i(u,v) by counting ordered spanning two-path covers with the first ending at u and the second starting at v. Treat the masks as 0/1 column vectors. Define

    2 s_i = I_x^T (M_i − M_i^T) I_a.

Then the literal count difference is

    H(plus) − H(minus) = 2 Σ_i s_i ∏_{j ≠ i} h_j.             (1)

**Proof.** Paths avoiding the switched edge already agree. Contract the used ordered port pair to q. With q deleted the remaining path segments visit the ordinal tiles in order, so a tile can be visited at most twice. Fix the tile labels along both segments, including the two labels at the neighbors of q, and let the local paths vary independently. Any tile adjacent to q on only one side contributes a sum of an appropriate W1 or W2 endpoint vector over an arrival/departure set. The equality of sums on the two masks also gives equality on their complements, so these patterns cancel. This includes patterns where q touches two different tiles, and patterns at an endpoint. When both q-neighbors belong to tile i, all earlier tiles are entirely to its left and later tiles entirely to its right. The two local counts are I_a^T M_i (1−I_x) and I_x^T M_i (1−I_a), with the first term for plus. Their difference is 2s_i by the matched W2 end sums. Other tiles each contribute one factor H(D_j).

Equation (1) is the ordinary product derivative with signed local inputs: [t] ∏(h_i+s_i t). It is not a claim about higher-order coefficients or independently switchable global signs. Swapping the two masks of a tile reverses s_i, without changing h_i or destroying the match. A transitive tape of singleton tiles can have arbitrary shared contacts and contributes factors 1, signed input 0.

## 2. Two charges larger than one cancel to a physical unit

Outgoing row masks use bit j for an arc to vertex j. The following literal endpoint-matched tiles have charge pairs (h,s):

| (h,s) | outgoing rows | I_x, I_a |
|---|---|---|
| (3,−1) | 2,4,1 | 4,2 |
| (161,+54) | 2,12,17,37,75,87,143,63 | 37,67 |
| (155,+1) | 14,28,88,96,73,23,35 | 5,18 |
| (89,+1) | 30,60,120,48,64,81,11 | 68,33 |

The first two combine to −161+3·54=1. In that order the literal counts are **224143,224145**. With their order reversed they are **178695,178697**. The old triangle-only identity with quotient count one is no longer the source of the unit: the nonport tournament here has H=483 and the unit comes from cancellation. Reordering changes the absolute center while leaving the unit invariant.

There are also four unbounded center families with the same signed charge. Between the two tiles insert m≥1 transitive vertices that either all lose to both ports (0) or all beat both ports (1). Write the minus count as A·3^m+B·2^m+C; add 1 for the center, 2 for plus:

| Tile order | Middle contact | A | B | C |
|---|---:|---:|---:|---:|
| 3,161 | 0 | 203826 | 26490 | −6227 |
| 3,161 | 1 | 96384 | 126288 | 4641 |
| 161,3 | 0 | 60544 | 98010 | 21391 |
| 161,3 | 1 | 237578 | −51588 | −8555 |

Removing x and a splits each Hamiltonian path into at most three pieces; within each, the middle block is a single run when present. A compressed trace using j middle runs lifts in j! S(m,j) ways. Sizes 1,2,3 therefore determine each row exactly; direct subset DP gives the three initial reads. This supplies exact computable movable centers, not interval coverage.

For the prescribed target 1091 use a tile with rows **(126,76,88,112,98,70,0)** and masks **(19,7)**. Its profile has h=15, s=1, and both masks read (15,0,99,309) on the four vectors. Append a single vertex beaten by every tile vertex but beating **both** ports. The two literal counts are (1089,1091), so its even center is 1090 and the positive digit hits the target. This is a new way to place that target, not a claim that 1091 previously had no realization. The vertex masks and counts are independently replayed in replay.py. The bare (155,+1) tile gives (2071,2073) as a further example not relying on a directed triangle module.

## 3. Which signed arithmetic transfers

* **Bézout / determinant.** For two tiles, call the half-gap s_1 h_2+s_2 h_1. A negative sign turns this into a determinant of two count/charge columns. The example is the exact 3×54−161=1 identity. In the archived [unimodular Bézout transport](../../../01-canon/theorems/THM-829-unimodular-bezout-owner-transport.md), arithmetic ownership likewise has to survive a transport; here the four endpoint equalities are the physical ownership condition.
* **Signed CRT and carries.** Write P=∏h_i. A unit requires all nontrivial h_i to be pairwise coprime: a shared prime divides every summand of (1). With pairwise coprime h_i, the signed residues must satisfy s_i·(P/h_i)≡1 (mod h_i) for every i. The residues alone give half-gap 1+kP, not 1. The actual bound |Σ s_i/h_i|≤1/2 would force k=0 for P>1. In the literal 3×54 example, cancellation pays the exact carry directly.
* **Balanced ternary / signed Collatz.** In the coordinate x=(H+1)/2, the three scalar operations H↦3H,3H+2,3H+4 become 3x−1,3x,3x+1. Sign selection from the residue gives the inverse digit, but an actual tournament requires a reusable endpoint state, not just a scalar equation. The archived [signed sibling ladders](../collatz_complement_ladders_20261007c.md) distinguish the analogous affine identity from its source/terminal guard. The four menu equalities allow the sign here; they do not establish ready children for all three ternary digits.
* **Pell / signed norm.** The square-triangular form solves a unit equation (2u+1)^2−8v²=1. The [central sign computation](../../../01-canon/theorems/THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle.md) retains the sign after quotienting. An analogous proposed deformation of the (161,54) tile changes (h,s) by (3t,t), keeping 3s−h=1. That is an exact invariant to search for in a physical tile family, not an existing graph operation.
* **Different roles of ±1.** Switching the global port arc reads C±s_total, where C is the midpoint. At s_total=1 the physical digit is solved. Placing C=N±1 is a separate absolute equation. Neither endpoint matching, signed scaling, CRT, nor a growing family of charges alone pays it.

## 4. Why both free tapes still do not place all centers

For the original triangle take a transitive quotient and let arbitrary exterior vertices precede or follow the triangle, with shared binary contacts to x,a. This is **every** unit case of the original shared-contact triangle construction, because its gap is 2H(Q) and a tournament has H=1 exactly when transitive. The [two-sided trace proof](TWO_SIDED.md) gives O(X^log_3(2) log(2X)) distinct counts through X, so infinitely many primes are missed even after both tapes vary.

The same upper bound holds for any fixed single endpoint-matched tile of nonzero response, with arbitrary shared contact words before and after it. To adapt the traces, split a Hamiltonian path of D into three nonempty runs when both sides have at least two gates. For a high right side alone, use a path through D minus two vertices followed by two singletons, choosing distinct successors for x,a; this is possible because each departure mask is nonempty and the two masks differ. The high left argument uses distinct predecessors. All this remains true if the other tape is placed whole at its outside. For a one-gate right tape use D' R0 x a d R1 in the minus orientation, with d a successor of a and D' a Hamiltonian path of the rest; reverse the port labels in the other orientation. For a one-gate left tape use L0 d x a L1 D', with d a predecessor of x. Nonempty proper masks follow from equality of the strictly positive W2 endpoint vectors, together with nonzero skew. A nonzero tile has at least three vertices. These substitutions give the same address-count argument. The bound therefore survives a fixed finite bank of single tiles. This does not bound arbitrary growing tiles or arbitrary compositions of many tiles.

The useful next question is now a two-coordinate physical one: obtain growing tiles whose matched masks retain 3s−h=1 (or another exact Bézout relation), **and** whose attainable even centers cover the necessary absolute intervals. The target hit at 1091 and the two-charge calculation address different halves of that question. Neither is extrapolated to universal coverage.

## Reproduction

Run python replay.py for direct permutation evaluation of the tile tables and independent literal subset DP for the full pairs. The two-sided proof and short path checks are in TWO_SIDED.md and verify_two_sided.py. Earlier work: [triangle switch](../three_visit_triangle_switch_20261117/README.md), [fixed-core suffix obstruction](../three_visit_triangle_switch_20261117/BINARY_DENSITY.md). For the published conjecture see [Chan and Pak, §6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).
