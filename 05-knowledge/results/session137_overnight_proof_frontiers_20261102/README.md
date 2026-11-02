# Session137: overnight proof frontiers

**The global Hamiltonian determinant bound is proved.** This session also proves polynomial weighted stability, all-modulus strong width-three local universality, exact modular context signatures, and structural limits on several exact-spectrum constructions. Finite classifications and independent checks are recorded separately from all-order proofs.

[Results and scope](SYNTHESIS.md) · [Verification and replay](VERIFICATION.md) · [First checkpoint](CHECKPOINT1.md)

## 1. Global determinant closure and counting certificates

For every nonempty tournament with adjacency matrix A,

```text
H(T) >= det(I+A) + |det A| + C_proper(T)
     >= max(1,s(T)) disc(T) >= disc(T),
```

where `C_proper` counts directed cycles not spanning T,
`disc(T)=det(I+A-A^T)/2^(n-1)`, and
`s(T)=1^T(I+A-A^T)^(-1)1`.
Equality in the first inequality holds exactly when no two directed triangles are vertex-disjoint. In particular, `H=disc` holds exactly for transitive tournaments. This closes the formerly open target in THM-1950 and THM-3729.

The weaker bound `H>=det(I+A)+|det A|` holds for every nonempty real weighted tournament matrix. The proper-cycle refinement is unweighted; the empty tournament and weighted order-two equality boundary are explicitly retained.

Four specified determinant deficit functions have canonical polynomial-time certificate verifiers and therefore belong to **#P**. This is not a general closure theorem for nonnegative GapP functions or an efficient algorithm for counting Hamiltonian paths.

[Global proof](discriminant/GLOBAL_PROOF.md) · [Independent proof audit](independent_global/GLOBAL_DETERMINANT_AUDIT.md) · [Counting construction](global_referee/COUNTING_GAPS.md)

## 2. Polynomial weighted stability

For `V=sum_(i<j) A_ij A_ji` and `Delta=H-det(I+A)-|det A|`,

```text
Delta >= beta_n V,                         n>=3,
beta_n = [24-2^(4-n)(n^2-n+3)]/[n(n-1)]
       >= 6/[n(n-1)],     beta_n ~ 24/n^2.
```

The coefficient is sharp at orders three and four: `beta_3=1`, `beta_4=3/4`. No higher-order optimality is asserted here. The bound gives an explicit unordered-edge rounding distance to a deterministic tournament, not a claim that the rounded tournament is transitive.

[Proof](discriminant/WEIGHTED_STABILITY.md) · [Independent audit](independent/WEIGHTED_POLYNOMIAL_AUDIT.md)

## 3. Exact-spectrum frontiers

| Result | Status and boundary | Read |
|---|---|---|
| Every parity-compatible residue modulo every integer has strong width-three witnesses with arbitrarily large distinct exact H values | **PROVED.** Local universality is not exact integer realization; inherited `w(13)=w(19)=5`, `w(31)=6` remain intact | [Width-three theorem](width_phase2/PROOF.md) |
| Exact homogeneous-context signatures use `p_T(c) mod [M/gcd(M,c!(c-1)!)]`; universal class counts at moduli 3,4,8,12 are 5,3,14,11 | **PROVED**, with finite witness banks; branching separates constructor grammars, not global spectral width | [Context theorem](width/PROOF.md) |
| Actual-source cone and lattice, original-port sealing, and finite-base normal forms | **PROVED.** Fixing both height h and interior L yields only `O(log X)` outputs; positivity is not realizability | [Constructors](constructors/README.md) |
| Sharp natural cycle-carrier packing cost; 29 connected carriers through cycle-budget 15 and 300 H values at total cost at most 15 | **PROVED / FINITE-EXACT.** The rigorous source-order bound makes the budget classification complete; H=31 separates minimum source order from minimum cycle cost | [Carrier classification](cycle_carriers/THEORY.md) |
| Common-child connectivity, edit divisors, H147 deletion-deck obstruction, and canonical deleted-triangle card reconstruction | **PROVED / FINITE-EXACT.** Scalar strong-quotient zero-rigidity through quotient order seven is producer **FINITE-EXACT**, not an independent full census | [Coverage](coverage/README.md) |
| First-missing-speed harmonic chambers, primitive boundary grids, all tied owners, and exact maximum certificates | **PROVED** for the specified second-band family, not arbitrary speed sets | [Harmonic theorem](harmonic_owners/README.md) |

## What remains open

- Full positive-odd H-spectrum completeness, apart from the established omissions 7 and 21.
- Whether exact spectral substitution width is globally bounded or unbounded.
- Scalar zero-rigidity for strong quotients of arbitrary order.
- LRC(14), arbitrary-speed harmonic closure, and a general below-chamber classification.

No finite census is extrapolated to any of these claims. No formal proof-assistant certification or external literature-priority claim is made.

## Reading the evidence

[Verification](VERIFICATION.md) distinguishes independent proofs, independent implementations, same-source replays, and finite scope. The [independent audit packet](independent/AUDIT.md) retains reviewed-source pins and exact controls. `width_phase2/probe.json` is **exploratory**, not a theorem dependency. First-checkpoint proof and replay files remain available. The [current concept board](BOARD.md), this README, and the [synthesis](SYNTHESIS.md) give the completed session status.
