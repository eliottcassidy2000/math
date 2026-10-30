# Pinned attachment counts and a smaller rigidity test

2026-10-30. **PROVED scoped / FINITE-EXACT / CONDITIONAL external application.**
Same-agent proof review; no independent formal certification or novelty claim.

The universal ramified tournament action now has a full cycle-index formula
after ANY nonempty set of vertices is named. Its geometry compresses to
occupied residue fibres and, within a single fibre, the LCA depth. This
extends the preceding stabilizer-order theorem to exact attachment counts
for every named palette, every joint binary-port arity, and every depth.

For one occupied fibre, increasing the allowed clock depth from ell-1
to ell obeys the short recurrence

    p B_ell(Q)=B_(ell-1)(Q)+(p-1)Q^(2p^(2k-ell)-p^(2k-2ell)).

Here B counts colorings under the specified pointwise stabilizer, on
N=p^(2k) vertices. Prescribing the colors at s pins divides B by Q^s.

The Boolean rigidity question needs LESS information. Every remaining
symmetry has a nonidentity order-p power among the independent shifts on
unoccupied fibres. Therefore, for uniform IID Q colors,

    Pr(rigid)=(1-Q^[-N(p-1)/p^2])^(p-j),

where j is the number of occupied residue fibres. LCA depth disappears
from this predicate but remains essential to enumeration. With pinned
vertices white, the minimum binary black support is exactly p-j.
At the rooted625 action this means four black vertices, or random-rigidity
probability (1-2^-100)^4. This kills the universal group after pinning;
it does not assert unrooted rigidity for every bank tournament.

Incoming Session112 also supplies the exact surviving unsigned reflection
phase. Combining its theorem with our coloring test gives a new full
phase correction: after a random minimal transversal and a random binary
attachment at625 vertices, unsigned symmetry survives with probability
exactly **2^-300/125**. The common phase and joint colours are retained;
this is not a claim about every individual bank's full group.

The requested [factoring preprint](https://hexagonmath.org/2610.00600v1)
supplies inspiration, not a formal dependency of those theorems. An
elementary lemma removes a power-of-two normalization exponent from
the modulus range needed to test a bounded count for zero. Session112
published this same lemma during our work and receives explicit credit;
our tests are a replay. Its use inside that engine remains an implementation
question. Conditional on the preprint's claimed
factoring machine, the existing mixed-multiplier aliquot inverse algorithm
has a subpolynomial validation bound for fixed multiplier.

The 1001 ->343 ->57 lane inherits Session78/81 reconstruction, Session82
uniqueness, Session108 squareclass boundaries, and Session111's578 merge.
It checks the actual prime7 query table and all19 sources of343; it does
not claim these were newly discovered or identify distinct edge types.

- [Proofs P1--P5](PROOFS.md): cycle indices, recurrence, exact directed
  and unsigned rigidity, and the same-occupancy/different-LCA hostile.
- [Arithmetic and primary-source audit](ARITHMETIC_AND_SOURCE.md): scale
  guards, prime-dilation sources, conditional complexity, connection ledger.
- [Audit](AUDIT.md), [checker](check.py), [receipt](receipt.json),
  [integrity manifest](manifest.json).

## Inheritance and concept board

Anchor: pinned ramified attachment enumeration. Niche: source-preserving
aliquot response validation. Wildcard: remove normalization height from
exact modular count queries. Live concepts: **cycle index; LCA; order-p
witness; joint palette; valuation/realization; unit modulus and height**.

Closest proved mechanism: U1--U2 in the
[preceding closure packet](../universal_ramified_closure_20261030/PROOFS.md).
Hostile: the actual81/27 pair, together with the earlier k=1 cycle-law
failure. Corrected near miss: a group base need not give injective direct
signatures, and rigidity is weaker than counting. Least-used sidecars:
the complete cycle partition after pinning, and the scale-prime support
of a modular count. Session78/81 already contains the arithmetic inverse;
recovering it prevented another pass through that saturated route.

Incoming `dce33fb05c`/`b51ef5cd73` supply unmarked affine-flat unions,
overlap correction and phase-aware complement transport. They remain the
right interface after forgetting the parent. Our fixed-parent formula
does not replace that interface. The native complement closure, full bank
automorphism classification, and unrestricted A000568 continuation grammar
remain **OPEN**. No new A000568 term is claimed.

Later incoming `fea3b07542`/`ecd0185d34` supplies the modular cutoff,
canonical saturated-cut ranking, signed norm tables and random phase axes.
The cutoff is credited; the axes enter P5 through their precise fixed-set
and coset data. Canonical ranking still needs an efficient saturated-cut
oracle, which the present fixed-parent formula does not automatically give.

The concept-board comparison has a sharp result: the same ell coordinate
is dispensable for rigidity and indispensable for orbit counts. Likewise
E is dispensable for the scaled zero predicate, while the unscaled height
bound is indispensable. Neither observation erases source or joint tensor
information for arbitrary later operations. This applies the existing
META-PATTERNS cards on controlled forgetting and recomputing certificates
after fibre-changing operations; no new general method card is promoted.
