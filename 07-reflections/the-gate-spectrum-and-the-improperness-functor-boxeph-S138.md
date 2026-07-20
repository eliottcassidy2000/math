# The gate spectrum: one functor under seven sessions of laws (boxeph-2026-07-19-S138)

Owner brief: hypothesize new constructs and combine ideas. Three definitions, one
identity, and two conjectures — each grounded in a measured fact from S130–S138.

## 1. The improperness functor I(k, Q, θ)

For modulus Q and threshold θ, let **I(k, Q, θ)** = the set of k-tuples of residues
mod Q admitting no witness scale (no b with all dist_Q(v·b) > θQ). Everything this
thread has touched is one instance:

- the **certificate rungs** (S130–S133): the spread lemmas at p ∈ {13,…,25} are
  emptiness/structure certificates for I(12 or 13, p, ·) — the small-Q end where
  per-scale conditions are decisive;
- the **level-1 seas** (kind-pasteur c88): I(13, p, 1/14), generic and cheap;
- the **tower gates** (S134–S138): I(13, 2^a·7^b·p, 1/14) with CRT restriction maps;
- the **n=12 gap program**: the same functor at threshold family 1/13.

The restriction maps I(k, Q, θ) → I(k, Q′, θ) for Q′ | Q (scale multiplication) are
what the ladder measurements probe: S137/S138's "refinement admit" = fiber
non-emptiness of I(13, 49p) → I(13, 7p).

## 2. The identity: loneliness is the phase boundary of I

M(V) ≥ θ ⟺ some modulus gives V's residues a witness ⟺ V ∉ I(13, Q, θ) for some Q.
So **M is exactly the phase boundary of the I-functor**, and the Pinch rung ladder
(M = D/s, S121–S123) is the discrete set of boundary values — the rung-quantization
observations (11-for-11, HYP-7945) say the boundary is attained only on the D/s grid.
The rung-certificate identity of S130 and the S-T tower are the same object seen at
small Q and large Q.

## 3. The banded CRT correspondence (the general slaving law)

For Q = c·p, every mixed scale imposes the relation Γ_c,p ⊂ Z/p × Z/c: the covering
runner's (mod-p, mod-c) pair lies in the band's CRT image. Structure theorem (proved
by inspection at c = 2, 7, 14, 49; general c the same argument): Γ has **density
exactly 2θ = 1/7**, and its fibers over Z/p are **APs with common difference p mod c**
of cardinality ~c/7 (1-of-7 at c=7, 2-of-14, 7-of-49 = the seven lifts of the forced
mod-7 class). The 1/7 density is conserved through every gate — the threshold never
dilutes and never concentrates; only the fiber GEOMETRY changes. The ×7 fiber being a
singleton is the apex-7 wall as pure CRT rigidity (HYP-7975).

## 4. The tower profile and the first-rung collapse conjecture

Define **Z_p(b) = |I(13, 7^b·p, 1/14)|** (death-star's |S| at tower height b).
Measured accessibility (one-sided, min-conflicts): at p=43 every rung is soft
(40/40, 20/20, 20/20); at p=191: ×14 3/15, ×7 2/15, **×49 refinement 2/2**.

**Conjecture (first-rung collapse):** Z_p(1)/(Z_p(0)·7^13) → 0 as p → ∞, while the
conditional ratios Z_p(b+1)/(Z_p(b)·7^13) for b ≥ 1 stay bounded below — the sea
thins at the sea→tower boundary, and refinement inside the tower stays cheap because
the slaving fibers are APs aligned with the previous rung's forced classes.
Testable: |S|-per-height should collapse at height 1 and plateau. If true, the
S-T pipeline's expensive object is ONLY the first 7-rung per prime, and wall (b)'s
7^13 cost concentrates where pruning by the fiber-AP structure is available.

## 5. The service-capacity program (the quantitative cage)

Combine the blocking budget (S130–S131: covering duties + blockers on 12–13 slots)
with the slaving law: each runner has a **service capacity** — the set of scale
constraints its residue pair can cover across ALL gates simultaneously. The cage,
if provable, is a counting theorem: 13 capacities cannot jointly service every gate
of every prime unless the tuple is ansatz-structured. The two-prime datum (S136:
coupling weak at the 14-level, 5.8%) says capacity binds only through the shared
2^a·7^b coordinates — so the theorem to attempt is a capacity bound at the FIRST
7-rung as a function of p (matching §4), not a cross-prime counting bound.

## Cross-links

HYP-7945/7950 (rung censuses) · HYP-7955/7975/7990/8010/8025 (gate laws + ladder) ·
kind-pasteur FINITE-CHECK ledger + c88 · death-star |S| experiment (the §4 test) ·
LRCMod19/23Spread.lean (the sieve template that formalizes §3 instances).
