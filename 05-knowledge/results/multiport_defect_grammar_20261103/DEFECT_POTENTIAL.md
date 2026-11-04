# The inner defect is a positive quadratic potential on the root's incoming set

2026-11-03. **PROVED reformulation of the inherited weighted-cut/rooted-defect
law**, with a new exact evaluation at the failed collar. This is not a global
decreasing rank for arbitrary repairs.

Source mechanisms:
[ADAPTIVE, rooted defect](../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md)
and [READY, Section1](../two_entry_funnels_20261103/READY.md). Their positive
pair term is retained explicitly here rather than hidden in one integer.

Fix a tournament D and adjoin a root r. Let I be the old vertices that beat r;
r beats every vertex outside I. For i,j in D write

    q=H(D), s_i=start_i(D), e_i=end_i(D), F_i=E_i(D,2),
    J_ij=ordered spanning two-path covers ending at i, then starting at j,
    kappa_ij=J_ij+J_ji (i != j).

Both paths counted by J are nonempty, and J_ii=0. The first slot in F is
nonempty, the second may be empty. All names refer to the fixed deletion D,
not to the tournament after r is adjoined. Then

    F_r(I) = q + sum_(i in I) F_i,
    b_r(I) = F_r(I)-H(I)
           = sum_(i in I) s_i + sum_({i,j} subset I) kappa_ij.       (Q1)

Proof. If D is empty, the resulting singleton has H=F_r=1 and b_r=0,
so (Q1) is immediate. Assume D is nonempty in the path partition below.
In a two-slot cover whose first path ends at r, either that path is
the singleton r, giving q, or deleting r leaves a first path ending at some
i in I, giving F_i. For Hamiltonian paths, separate r at the beginning,
end, or interior:

    H(I)=sum_(j notin I) s_j + sum_(i in I) e_i
         +sum_(i in I,j notin I) J_ij.

Use q=sum_i s_i and F_i=e_i+sum_j J_ij and subtract. This proves (Q1),
including empty I and full I; no strongness assumption is needed.

Thus b is monotone and supermodular in the incoming set, F is modular,
and H is submodular. Adding one new inport v has the exact nonnegative cost

    b(I union {v})-b(I)=s_v+sum_(i in I) kappa_iv,                  (Q2)

and the mixed second difference in distinct v,w is exactly kappa_vw>=0.
All third mixed differences of the root-star bits vanish. These statements
hold while D is fixed. Flipping edges *inside D* changes the coefficients;
the path-cover transport in [transport](transport/PROOF.md) supplies those
changes. It would be incorrect to infer monotonicity for arbitrary repairs.

## The actual 2960 is a pure interaction between two named ports

Take the failed grandchild suffix B14 of H1537279, with H(B)=80319,
root r=32 and incoming set I={22,33}. Its deletion D=B-32 has H(D)=13043.
Exact deletion queries give

    start_22(D)=start_33(D)=0,
    J_D(22,33)=890, J_D(33,22)=2070.

Therefore (Q1) reads

    b_32(B)=0+0+890+2070=2960.

The difference is neither an unexplained rounding residue nor a scalar that
can be silently reset. It is the joint interaction of the two incoming ports.
For the same fixed D, either singleton incoming set has zero defect. That
does not preserve H, the outer attachment, or a legal return word. The
[neutral-repair audit](neutral/README.md) shows that no allowed three-suffix-edge
repair restores any inner zero-ready mark while preserving the full1537279
count, even with arbitrary re-marking and nonshrinking transitive prefixes.

The incoming Session181
[multiport attachment identity](../session181_multiport_defect_20261104/multiport/PROOF.md)
is a useful specialization. If D=A->B, I=W lies in B, and A is nonempty,
all s_w(D)=0 and J_D(u,v)=H(A)J_B(u,v). Hence

    b_r = H(A) sum_(u,v in W) J_B(u,v).

This recovers its transported defect, including its distinct actual +4
example with outer defect3690. It does not identify that graph with the
zero-ready +4 graph from the previous packet.

## Independent finite controls

`python 05-knowledge/results/multiport_defect_grammar_20261103/check_defect.py`
checks (Q1) against fresh endpoint DP for every labelled deletion of order
0..4 and every incoming set. It also verifies (Q2), all mixed quadratic
coefficients, the actual deleted collar values, and the distinction between
deletion-level J_D(22,33) and full-tournament J_B(32,33). Assertions use
explicit exceptions and remain enabled under Python -O.
