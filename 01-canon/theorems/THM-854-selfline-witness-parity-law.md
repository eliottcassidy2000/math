---
id: THM-854
title: The self-line witness parity and path-holonomy laws
status: PROVED ALL n (flip-vector iteration and witness-order parity) + FINITE-EXACT n=5..8 (cycle type, square obstruction, and two-path holonomy)
source: opus-2026-07-15-S312; n=8 path-holonomy extension codex-2026-07-15-S15
depends_on:
  - THM-852 (Klein-four skeleton; the n=8 refutation)
related: [THM-790 (leg law), THM-849, THM-855 (F1-F5), MISTAKE-150, codex-S14 H-drift]
verification:
  - 05-knowledge/results/selfk_type_ledger_fixed_opus_S312.out
  - 05-knowledge/results/selfk_debug_witness_opus_S312.out
  - 04-computation/path_relative_kappa_holonomy_codex_S15.py
  - 05-knowledge/results/path_relative_kappa_holonomy_codex_S15.out
---

# THM-854 — the self-line witness parity and path-holonomy laws

**Setup.** Encode a tournament X on [n] as a GF(2) vector over pair-positions
(bit per pair with a fixed orientation convention). Reversal is X ↦ X ⊕ 𝟙
(flip every pair); the complement tiling is X ↦ X ⊕ 𝟙 ⊕ p₀ = X^{p₀} (flip
every pair OFF the base path p₀). A permutation π acts affinely on tournament
bit-vectors; write `L_pi` for its linear relabelling action on pair masks. A
**witness** of a self-line at tiling t is π with π·X = X^{p₀}, X = T(t).

## (1) The iteration identity

π^k·X = X ⊕ [k odd]·𝟙 ⊕ (p₀ Δ L_pi p₀ Δ ... Δ L_pi^{k-1}p₀)
for all k ≥ 0.

*Proof.* Induction: `pi·(X xor F)=(pi·X) xor L_pi(F)` for every flip-set
`F`, while `L_pi(1)=1`. ∎

In particular **σ² realizes the pair-flip p₀ Δ L_sigma p₀** for any witness
σ. So `sigma^2 in Aut(X)` iff `L_sigma p_0=p_0`, meaning that σ preserves the
base path's undirected pair-set. Its only setwise automorphisms are the
identity and path reflection. For `n>=3`, the identity cannot witness because
the off-path flip-set is nonempty, so `X != X^{p_0}`. Path reflection reverses
the standard orientation on the path arcs, while both `X` and `X^{p_0}`
contain the same standard directed path, so it cannot witness either. Thus
**no black self-line witness at `n>=3` has `sigma^2 in Aut(X)`**: the premise
of THM-852(ii)'s involution corollary never holds. (Verified directly through
`n=8`, where every direct black realizer also has order greater than two.)

## (2) The order equation

Setting k = ord(π) in (1): **[k odd]·𝟙 = Δ_{j<k} L_pi^j p₀.**

- k odd ⟹ the π-orbit of the base path XORs to the COMPLETE GRAPH: a
  **mod-2 path-decomposition of K_n** by ≤ k translates of one path.
- k even ⟹ the translates XOR to ∅ (each pair covered evenly).

There is a compact cyclic-module form. Let `C=<L_pi>` act linearly on the
pair-mask module `V=F_2^(E(K_n))`, put
`N_(L_pi)=1+L_pi+...+L_pi^(k-1)`, and keep the affine displacement
`delta_pi(X):=pi·X+X` separate. The witness, square-holonomy, and order
equations are respectively

```text
delta_pi(X)=1+p_0,          D_pi=(1+L_pi)p_0,
N_(L_pi) p_0=[k odd]1.                                  (2.1)
```

Thus existence is an affine-displacement problem, the failed square hypothesis
is the one-step linear path coboundary, and witness parity is its cyclic norm
obstruction.  This connects the self-line object directly to the defect
algebra without identifying its class quotient with the underlying module.

## (3) The parity law

Counting cardinalities mod 2 in (2): k odd requires
C(n,2) ≡ k(n−1) ≡ n−1 (mod 2), i.e. **n(n−1)/2 ≡ n−1 (mod 2)**, which holds
iff **n ≡ 1 or 2 (mod 4)**.

> **Odd-order witnesses can exist only at n ≡ 1, 2 (mod 4); at
> n ≡ 0, 3 (mod 4) every self-line witness has even order.**

## (4) The computed direct-realizer spectrum

| n | witness cycle types (mass) | orders | n mod 4 | odd allowed? |
|---|---------------------------|--------|---------|--------------|
| 5 | (1,1,3): 4, (5,): 4       | 3, 5   | 1       | yes — and ALL odd |
| 6 | (1,5): 4, (2,4): 8        | 5, 4   | 2       | yes — mixed |
| 7 | (1,6): 32, (2,5): 48, (3,4): 8 | 6, 10, 12 | 3 | no — all even ✓ |
| 8 | (2,2,4): 128, (1,2,5): 188, (1,3,4): 32, (8,): 56 | 4, 10, 12, 8 | 0 | no — all even ✓ |

At n=7 every witness has exactly TWO cycles (all two-part partitions of 7
appear).  No direct realizer in the audited range is an involution.  All black
endpoint tournaments in this range are asymmetric, so each endpoint has one
direct realizer and the displayed masses equal the black endpoint counts.
In general the witness-pair mass is instead
`sum_(t in Q_n^K)|Aut(T_t)|` and need not equal `2 selfK(n)`.

## (5) The path-holonomy stalk

The operation is path-relative, with exact equivariance

```text
sigma kappa_P = kappa_(sigma P) sigma.                     (5.1)
```

For a direct realizer `sigma T=kappa_P T`, its square therefore satisfies

```text
sigma^2 T=Flip_(D_sigma)T,
D_sigma=E(P) symmetric_difference E(sigma P).             (5.2)
```

The tournament-class quotient forgets `D_sigma`.  That is precisely why the
conditional coset lemma in THM-852 could not be applied: `sigma^2` is an
automorphism iff `D_sigma` is empty.  An independent complete witness search
finds this never occurs for a black direct realizer at `n=5,6,7,8`.

At `n=8`, the exact holonomy-size census is

```text
|D_sigma|             6    8   10   12   14
endpoint realizers  136    8   84  128   48
Klein-four orbits    34    2   21   32   12.              (5.3)
```

The joint spectrum is stricter than cycle type.  For example, 8-cycles occur
at every holonomy size in (5.3), while type `(5,2,1)` occurs at sizes 10 and
12.  The complete joint table is stored in the verifier output.  Thus the
minimal witness-level object is not merely a permutation type or tournament
class; it must retain the marked path pair `(P,sigma P)`, or at least its
symmetric-difference flip holonomy.

## Reading

The black endpoint counts `2 selfK=8,12,88,404` at `n=5..8` are not `SC(n)`;
the equality fails at `n=8`.  In this finite range every endpoint is
asymmetric, so these numbers also equal the direct-realizer pair mass.  The
all-size pair mass is the separately weighted quantity
`sum_(t in Q_n^K)|Aut(T_t)|`.

The parity law decomposes that mass by witness order with a mod-four dichotomy
in `n`. Any odd-order witness in the permitted classes `n=1,2 mod 4` yields a
mod-two path decomposition of `K_n`; the parity law does not assert that such
a witness exists for every permitted `n`. At `n=0,3 mod 4`, every witness is
forced into an even-order channel. Equation
(5.3) adds the path-relative transport coordinate destroyed by the class and
Klein quotients.  The next counting problem is consequently joint: enumerate
witness cycle type together with path-orbit XOR, rather than seek a class-only
decoder for the self-line counts `4,6,44,202` or the Klein-orbit counts
`2,3,22,101`.  ∎
