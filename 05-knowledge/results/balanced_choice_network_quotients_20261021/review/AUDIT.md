# Independent audit of the choice-network quotient

Date: 2026-10-21 UTC. Verdict: **PASS after two scope qualifications**.

The corrected producer proof, its unchanged checker and both saved receipts pass this audit. All network formulas and the critical-group presentation are valid with their stated domains. The associated root synthesis passes for its choice-network paragraphs after adding the positive-depth qualification. This is not an audit of the synthesis's other mathematical lanes, an identification with the tournament metagraph, or a literature-priority claim.

## Frozen target and corrections

The final audited producer proof is `../session20-choice-networks/PROOF.md`, SHA256 `09959ef9f50e716e0df6965d3c86ff564a89ce9a3f4a7f629667a7e02120f7d0`.
The unchanged producer checker has SHA256 `117ff01b044c4dba9cd99928278ccce027726ed0ad6880aa2ef8cd3bcee5de0b`.
Both producer receipt files have SHA256 `97835a5189e0b0734f4baf7e70d948707ec07b4928aa4453823f3795a920f1fa`.
The choice-network statements in the final root synthesis are bound to SHA256 `efd074c6e9b29c7526957128bf265b3e730090c64f671bac129d8ce92317cf13`.

Two review findings were corrected by the producer, not silently changed by this reviewer:

1. The initial proof said loop mass could not be recovered from the electrical Laplacian unless separately retained. This was too strong for the specified pointed complete-source family: the known exceptional weight one gives `w_x = -L_(x,infinity)`, from which the source assumption recovers `binom(w_x,2)`. The repaired paragraph correctly distinguishes a Laplacian's insensitivity to arbitrary loops from reconstruction using known fibre/source provenance. The earlier proof hash `bfd048cbee73f300e241bfddb7dc09bedf02bc053b7bbd0e49f4a77598ae1180` is retained in the manifest as superseded, not as approved current content.
2. The synthesis's two-primary formula needed `d >= 1`. At `s=1` the general Smith formula still holds, but `Q` need not be odd. For example, `q=4,s=1` gives `K_4`, whose critical group has order 16, not a trivial two-primary subgroup. The main proof already had the correct positive-depth hypothesis.

The pin checker rejected the changed proof before its new hash was approved. `PIN_REJECTION.json` records that deliberate version reconciliation. No producer file, repository file, or Git state was modified by this reviewer. The producer's original “awaiting independent audit” heading can be superseded by a publication wrapper pointing to this frozen PASS; changing the target proof again would require a new pin/review.

## Inheritance and incoming overlap

`SOURCE_PINS.json` binds the complete source bytes and identifies roles:

- `zero_heavy_ports_20261020g/microcosm.md`, sections 2 and 4: the actual pointed choice square, reflection/blue diagonal, and central-pair deletion.
- `balanced_fibre_endpoint_extensions_20261020/BALANCED_DELETION_EXTENSIONS.md`, sections 1, 1.1 and 2: arbitrary reflected-depth deletion, its product fibres, and collapsed internal black paths.
- `strong_blue_games_20261020f/double_return_counts.md`: the upstream path decoder.
- Incoming `metagraph_fibre_network_20261020k/README.md`, section 4, reported by the parent as commit `118eec39549`: one-step weighted contraction, resistance, held Markov law and conditional entropy were already proved there. Its generic wiring identity already distinguishes shorting from Schur complement. The present packet continues that work with general depth/weights, a full integral critical-group presentation and an explicit complete-graph Kron formula; it must not present the overlapping one-step formulas as new discoveries. The corrected root synthesis explicitly attributes this overlap.

The concurrent balanced port-defect and conductor-tower documents are additionally pinned only as boundary/comparison sources. Their general results are not independently certified by this network audit. In particular, the arithmetic comparator `G_16 = C4 x C2` is supplied by that separate lane; no group identification is inferred from equal orders.

Actual balanced source parameters satisfy `q=2^(k-1)+1 >= 3`. The standalone positive-integer network theorem legitimately permits all `q >= 2, s >= 1`; its synthetic `q=2` tests are not claims that this parameter occurs in the original tournament family.

## Proof audit

### Contraction, energy, and lost loops

Every fine edge between distinct fibres contributes one unit, giving `c_xy=w_x w_y`. Each internal pair becomes a loop. Expanding the edge energy yields `L=Q diag(w)-w w^T`; all loop vectors are zero in this electrical form. The diagonal-parent count is `w_x^2`, split into `w_x` blue lifts and `w_x(w_x-1)` ordered black lifts. Passing modulo swap divides the black contribution by two. The exceptional `1` is a diagonal state, not a third nonloop edge-weight type. The repaired provenance qualification is necessary and sufficient.

### Resistance and spanning trees

For `b=e_a-e_b`, the displayed potential `u_i=b_i/(Q w_i)` has weighted mean zero and satisfies `Lu=b`. Its potential difference proves the resistance formula. The root cofactor is a positive diagonal matrix minus one rank-one matrix; the determinant factor is exactly `w_r/Q`, giving `Q^(q-2) product(w)`. This works for arbitrary positive real weights as a weighted tree sum, and for integer fibres as a multigraph tree count. Loops never contribute. All denominators are nonzero in the stated domain.

### Critical group and an explicit integral certificate

With `m=q-1`, the integer reduced matrix is `A=s(QI-sJ)` and `Q-sm=1`. Its column relations satisfy `sum_i r_i=sS`. Therefore replacing the relation set by `sS` and all `sQ e_i` is an equality of integral lattices, not merely a rational row reduction.

There is also a concrete matrix certificate. Set

- `U=(e_1,...,e_(m-1),S)`;
- `V=(e_1+s*1,...,e_(m-1)+s*1,1)`;
- `D=diag(sQ,...,sQ,s)`.

Then `det U=det V=1` and `AV=UD`. For `V`, subtracting `s` times its last column from every earlier column gives `U`. Thus `A` is integrally equivalent to `D`. Reorder diagonal entries to obtain the Smith chain `s,sQ,...,sQ`; the divisibility condition is automatic. This includes `m=1` and `s=1`. For `s=2^d,d>=1`, `Q` is odd and the two-primary subgroup is `(C_(2^d))^m`. The odd critical group of the fine `K_Q` follows from the same formula at unit weights. Preservation of restricted real energies does not assert preservation of the integral firing lattice.

### Kron reduction

For `h=Q-q>0`, the interior block has eigenvalue `q` on the all-ones vector, giving `1^T L_II^-1 1=h/q`. The Schur complement is therefore `Q I_q-(Q/q)J_q`. For `h=0`, it reduces to the original complete-graph matrix. The response and pair resistance do not depend on which representatives were selected. This procedure frees interior potentials; it is not the constant-fibre restriction. Clearing denominators changes the firing lattice and physical conductance scale unless a separate lattice convention is declared.

### Strong lumpability and loop normalization

Every microstate in fibre `x` has `w_y-delta_xy` available neighbours in fibre `y`; this proves strong lumpability directly. The proposed stationary law satisfies detailed balance and sums to one. To realize the same walk on the contracted multigraph, each loop contributes two directed incidences. Consequently the degree is `w_x(Q-1)` and the hold is `(w_x-1)/(Q-1)`. Dropping loops changes both the transition law and generally its stationary distribution. No corresponding claim is made for physical tournament flips.

### Conditional entropy and blue correction

Under a uniform child ordered pair, the parent mass is `w_x w_y/Q^2`. Its conditional fibre is uniform. When `s=2^d`, a finite coordinate contributes `d` bits with probability `ms/Q`, so the two-coordinate conditional entropy is `2dms/Q`. A uniform-parent sampling law instead gives `2dm/q`; conflating them changes the result. Enumerating the child diagonal and the inverse image of the parent diagonal gives probabilities `1/Q` and `(ms^2+1)/Q^2`, with difference `ms(s-1)/Q^2`. These formulas concern ordered path states, not unordered choice edges.

## Reproduction and independent finite scope

From this workspace:

    python3 session20-choice-networks-review/check_review.py > review.out
    python3 -O session20-choice-networks-review/check_review.py > review.optimized.out
    cmp review.out review.optimized.out

The checker accepts `--workspace`, `--producer`, and `--repo` to identify relocated source roots. The manifest's workspace-relative companion paths must remain available for full source-pin verification. The checker imports no producer or repository Python code; producer normal/optimized runs are isolated subprocess replays. Independent mathematical routines use literal edge/path enumeration, a fraction-free Bareiss determinant, a separate exact linear-system solver, explicit unimodular transforms, and Prüfer tree enumeration. All validation failures are explicit exceptions and remain active with `-O`.

Both review modes produce the same **29,865-check PASS** receipt. Both producer modes reproduce the frozen **1,328-check** receipts byte for byte.

The independent scope is:

- 360 arbitrary positive-integer weight profiles: `q=2,...,5`, each weight in `{1,2,3}`. Literal contraction, energy, every root cofactor, Prüfer weighted-tree sum, all pair resistances, every microstate's quotient transitions, loop-dart normalization, reversibility and stationarity.
- 84 integral matrix certificates: `m=1,...,12`, `s in {1,2,3,4,5,8,16}`. Explicit `AV=UD` and both unimodular determinants, not just matching group order.
- 30 Kron networks: `q=2,...,6`, `s=1,...,6`, retaining one actual vertex per fibre. Boundary responses are obtained by solving each Dirichlet basis problem; all pair resistances are checked.
- 36 entropy/blue cases: `m=1,...,6`, `d=0,...,5`, including direct enumeration of every ordered child pair and comparison with the different uniform-parent law.
- Inherited-carrier corroboration by literal Hamiltonian DFS for `k=2,...,5`, giving `9,25,81,289` paths. These match the product decoder, the reflection swap and fixed-point diagonal. All source/target pairs in this range with positive deletion depth give **1,054 actual path deletions**, with exact product fibres.

Retained hostile witnesses distinguish shorting from Kron; holding from loopless jumps; child-uniform from parent-uniform entropy; ordered from unordered collapse; the full critical group from its two-primary subgroup; the minimal V4 coincidence from larger rank/depth; known-source reconstruction from arbitrary loop blindness; the `d=0` boundary; and denominator-cleared Kron from the original integer contraction. At `q=3,s=2`, clearing Kron's denominator produces tree count 75 rather than the contraction's 20.

The producer checker was sound but had weaker coverage than its prose might suggest: its entropy check was only the literal `8/5` calculation, and its blue correction used aggregate formulas rather than explicit ordered-pair enumeration. Those are testing limitations, not faults in the general proofs; the independent checks above close that finite-control gap.

## Limits and final disposition

The all-parameter conclusions rest on the displayed algebraic proofs, not extrapolation from these finite ranges. No unrestricted tournament census, metagraph homomorphism, full conductor-theory audit, complete classification, or novelty claim is approved. The physical port/flip transport obligation and the balanced two-seam sidecar remain separate. Within these boundaries, the corrected choice-network packet and the scoped synthesis statements are approved.
