---
id: HYP-9762
title: "No linear tape: in a one-apex host whose tape spells the binary digits of X with fixed digit gadgets, the Hamiltonian-path count is never affine in X (nor affine plus a digit-sum term); exact absolute-count hosts of arithmetic-progression type therefore need a unipotent (eigenvalue-1 Jordan) element in some port algebra, which positive path counts do not supply"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-15-hosts, wave 31 of chanpak_waves_20261113.
  FINITE-EXACT: no affine host among 442,598 one-apex hosts with digit gadgets of at most 3 vertices (L = 3), and no
  affine-plus-digit-sum host among 563,500 (w31_linear_tape_search.c, w31_affine_digitsum_tape_search.c).
  ARGUMENT for the mechanism: the cut formula H(T + x) = H(T) + 2W(X), with W the cut function of the nonnegative
  path-weight matrix Pi, has quadratic part -sum_{u,v in X} Pi(u,v), and intermediate gadgets scale long-range weights
  multiplicatively.
session: mac-mini-2026-11-15-hosts
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 31)
related:
  - THM-5192 (A) (carrier forms, surcharge p), HYP-9752 (row collisions), THM-5212 (Hankel laws with negative coefficients),
    ribbon THM-4880 (H) (no AP-inserting blocks in span 2)
---

# HYP-9762 - why a prescribed absolute count cannot come from a linear tape

**The idea that fails.** Suppose two digit gadgets B0, B1 made the host count exactly affine, H(host(X)) = D·X + E.
Then a handful of such progressions covering every odd residue would give every large N directly. That would be an
exact absolute-count host and a proof of Chan–Pak.

**What happens.**
- *Exhaustive check.* No gadget pair with at most 3 vertices, under any core marks or guard marks, makes the count
  affine. Allowing an extra digit-sum term β·s₂(X) does not help either.
- *Mechanism.* Write the apex as x with out-set X. Then H(T + x) = H(T) + 2W(X). Here W is the cut function of the
  nonnegative matrix Π(u, v) = Σ over even-order paths u ⇝ v of H(T − V(P)).
- *Why the quadratic part survives.* The quadratic part of a cut function is −Σ_{u,v∈X} Π(u,v), which positivity
  forbids from vanishing. Restricted to gadget patterns it could cancel only through exact cross-gadget identities.
  Yet each intermediate gadget multiplies long-range weights through H(T − V(P)), so all digit orders interact.

**What a circumvention needs.** A counter v ↦ 2v + b is a pair of port matrices with eigenvalues {2, 1} and a
Jordan-type coupling: a unipotent element. Two facts point the other way:
- every floor unit has a single nonzero eigenvalue (THM-5212 (D));
- span-2 blocks of length ≤ 4 have eigenvalue 1 at most simple on the reachable space (THM-4880 (H)).

The Hankel laws do have negative coefficients (H(WAB) = 4H(WA) − 3H(W)). So cancellation is possible in principle
in port algebras with several ports, and that is where a counter would have to be found.

**Conjecture.** No one-apex host with fixed digit gadgets has a count affine in the tape value. More strongly, no
finite family of fixed-gadget hosts contains an exact arithmetic progression of counts of unbounded length.

## Update (wave 32): the counter search in the full span-3 representation, and the sparse linear tapes

**No counter among short span-3 blocks (FINITE-EXACT, 62-dimensional representation, exact modulo two primes).**
- *No unipotent element.* Of the 4,680 span-3 blocks of length ≤ 4, 233 have eigenvalue 1 and none has a Jordan
  block there. This extends THM-4880 (H) from span 2 to span 3.
- *Eigenvalue 1 in strong blocks.* 141 of those blocks are bulk-strong. For them eigenvalue 1 is a semisimple
  sub-dominant eigenvalue: a constant offset in H(p u^k s), not a counter.
- *Shared fixed vectors.* A two-digit counter needs two strong blocks with a common fixed vector. Only 8 pairs of
  distinct primitive strong blocks share one, and each is an H-twin pair: in every tested context the digit weights
  are 0, because the blocks differ by one arc that never changes the count.

**Where the surcharge is linear: sparse tapes (PROVED).**
- *Padding law.* p(2^k u + 1) = p(u) + (2^(k−1) − 1)u + (u mod 2) for k ≥ 2, by p(2w + 1) = p(w) and
  p(2^j u) = p(u) + (2^j − 1)u + (u mod 2).
- *The triple.* For u = 1 this gives p(2^k + 1) = 2^(k−1) = (X − 1)/2. The binary addresses 1, 1001, 1000001,
  1000000001 (X = 1, 9, 65, 513) are exactly where the carrier count is affine in X.
- *Padded payloads.* An address 2^k u + 1 carries a surcharge only on the short payload u, with slope
  2^k(qY + t) + a(2^(k−1) − 1).

So linear tapes exist only on geometrically sparse address families. That is the precise form of the obstruction.
