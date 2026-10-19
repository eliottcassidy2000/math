# A parity obstruction to uniform protected edges in a common-pair affine lift

Date: 2026-10-19. Proved for the architecture stated below. This is a constraint on strengthening an existing construction, not a challenge to Gerbicz's established unrestricted interval theorem.

## 1. Architecture and theorem

Let n be large enough for two fixed vertices u,v to lie below 49n+24. Let S₀ and S₁ be square-sum paths on [1,n] and [1,n+1] respectively. Assume they are parity-compatible: every common integer occupies positions of the same parity in the two paths. For c in [−24,24], their affine blocks are

    T(c,i)_j = 49 S_i[j] + (−1)^j c,       i=0,1.

Changing the global sign convention for j merely relabels c. A construction for each N=49n+d, 24<=d<=72, selects exactly one of T(c,0),T(c,1) for each c, uses each selected block intact (possibly reversed), adds singleton vertices 1,...,24, and concatenates these blocks and singletons into a square-sum path or cycle on [1,N]. The same pair S₀,S₁ is used throughout this 49-value residue range. These are the relevant structural features of the affine-lift method described in [Gerbicz's original discussion](https://www.mersenneforum.org/node/17331/page2); the result below is conditional on this explicit model and does not require access to his particular glue tables.

**Theorem.** Write the unique centered residue decompositions

    u=49r+δ,       v=49t+ε,       −24<=δ,ε<=24.

Suppose u,v>24, r!=t, |δ|<24, |ε|<24, and u+v is a square whose positive root is not divisible by 7. Then no construction of the stated form can contain the fixed edge u–v at both N=49n+24 and N=49n+72. In particular it cannot preserve this edge for every d in [24,72].

The theorem allows arbitrary endpoint values and arbitrary residue-dependent glue choices. It only requires the common pair and intact-block architecture.

## 2. Proof

Every edge internal to a block has sum 49 times a square. The prescribed edge therefore cannot be internal. In every output containing it, both its vertices must be endpoints of their respective blocks, because u,v>24 excludes singleton connectors.

At d=24, all blocks must use S₀: each S₁ block has one extra vertex corresponding to n+1, and that vertex is at least 49(n+1)−24=49n+25, outside [1,N]. Thus r and t must both be endpoints of S₀.

At d=72, the 49 possible extra vertices are the distinct integers

    49(n+1)−24, ..., 49(n+1)+24.

Their block labels are permuted by the parity sign of n+1. Exactly the largest extra vertex lies outside [1,N]; it belongs to one of the two extreme blocks c=±24. Every block with |c|<24 must therefore use S₁. Since a fixed vertex's block label has absolute value equal to its centered residue, both u and v belong to these non-extreme blocks. Hence r and t must both be endpoints of S₁.

But two parity-compatible paths of consecutive lengths cannot have the same two endpoint vertices. The even-length path has one endpoint in an even position and one in an odd position. The odd-length path has both endpoints in odd positions. The even-position endpoint therefore cannot be an endpoint of the odd-length path while preserving its positional parity. This contradicts r!=t and proves the theorem.

The strict centered-residue hypotheses matter: one extreme block still uses S₀ at d=72, so the argument does not automatically cover edges with |δ|=24 or |ε|=24.

## 3. Exact seed obstruction

The AP seed edge 482–3362 satisfies

    482+3362=62²,
    482=49·10−8,
    3362=49·69−19.

Both centered residues are non-extreme, the quotients differ, and 7 does not divide 62. Therefore the theorem rules out preserving this edge throughout a full residue range using one common parity-compatible input pair, even if one enlarges the permitted endpoint constants.

The second seed edge has the different descent

    3362+6242=98²,
    3362=49·69−19,
    6242=49·127+19,
    69+127=14²,
    69=49·1+20,
    127=49·3−20.

It can be internal through two levels, descending finally to 1–3. This does not remove the first-edge obstruction.

## 4. What remains possible

This obstruction does not apply to a construction which chooses different input pairs for different output residues, splits or edits internal blocks, adds extra exceptional blocks, uses a different dilation architecture, or chooses scale-dependent absorber hosts rather than requiring a fixed edge at every level. Those changes require their own exact coverage and parity proofs. It also does not obstruct a path or cycle containing a prescribed seed for any individual N.

A finite base list containing the desired seed is insufficient by itself to establish an all-large-N prescribed-seed theorem in the common-pair architecture. A valid induction must explicitly verify protected-edge state transitions as well as ordinary vertex coverage.

## 5. Reproduction

`verify_affine_lift_protected_edge_obstruction_20261019.py` uses only the Python standard library. It checks the forced block selections at the two residue endpoints for many n, the centered seed decompositions, and the endpoint-parity contradiction. The infinite statement follows from the proof, rather than these finite controls.

From the repository root run `python 05-knowledge/experiments/verify_affine_lift_protected_edge_obstruction_20261019.py`. The receipt is [affine_lift_protected_edge_controls.json](affine_lift_protected_edge_controls.json).
