# Whole matching flips are exactly optimizable, but not enough

**PROVED; finite calibration included.** Fix d branches and n=d−1 sheet labels per branch. The block matrix B has diagonal zero and reciprocal permutation off-diagonal blocks P_ij. Define its exact defect by

    F = B²+B−n I−(J_d−I_d)⊗J_n.

The diagonal blocks are zero. Appending the root and prescribed parent vertices makes the full degree-d graph's defect equal to F on leaves and zero elsewhere. Thus the zero layer of this metagraph is precisely the Moore completions.

Replace one matching P=P_ij by any permutation Q, with transpose on the reverse block, and put D=Q−P. The only changes of F, up to transpose, are

    (i,j):  D,
    (i,k):  D P_jk,
    (j,k):  D^T P_ik       (k ≠ i,j).

No squared term survives in a diagonal block because PP^T=QQ^T=I. Orthogonality gives

    ||F_new−F||² = 2(2d−3)||D||² = 4(2d−3)s,

where s is the number of labels with changed image. This norm formula is true away from an exact completion as well. From an exact completion F=0 it is the energy cost. There is one new triangle per inserted leaf edge: its pair was formerly a nonedge with one alternate two-step path, the middle branch is not the edited pair, and no triangle uses two added matching edges. Combining with the [cycle-energy law](moore57-multiedit-energy.md) gives exactly s triangles and (d−3)s four-cycles. A single transposition at degree 57 therefore costs 888 and gives 2 triangles and 108 four-cycles.

## Exact cost of the best replacement

With Frobenius inner products on single matrices define

    W = F_ij + Σ_(k≠i,j) [F_ik P_jk^T + P_ik F_jk^T] − (2d−3)P.

The full gain identity, with E=||F||², is

    E_new − E = 4⟨W,Q−P⟩.

Indeed, take the inner product of F with the three changed-block expressions, double for transposes, and use ||Q−P||²=−2⟨P,Q−P⟩. The best replacement over all permutations is one linear assignment problem. The anchor-gauge derangement constraint just forbids the diagonal, so it is also a linear assignment. If the admissible permutations instead form the even Gewirtz catalogue, minimize the same linear cost by scanning the catalogue; unrestricted assignment need not respect that group.

An exact completion is a strict local minimum. The converse fails: on the explicit Hoffman–Singleton degree-seven chart, [seam-descent.py](seam-descent.py) checks the formula against all 720 permutations in a fixed test position, then starts 32 deterministic random anchor-gauged states and repeatedly makes the best strictly improving whole-matching replacement. All 32 terminate at nonzero energy (smallest 488). This certifies local traps, not their abundance or their global energy frontier.

For the even Gewirtz catalogue at d=57, every nonidentity relative move within the group displaces at least 48 labels, by [even-seed.py](even-seed.py). Thus the minimum gap from a zero to a different admissible single-seam state is at least 21,312 in this score. This is about one move near a hypothetical exact state, not about a trade between two exact states.
