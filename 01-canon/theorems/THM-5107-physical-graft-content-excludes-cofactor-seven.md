---
id: THM-5107
title: Physical graft content excludes cofactor seven
status: PROVED all proper rooted grafts; independently audited with separate finite exact controls. Chan–Pak completeness remains open.
source: 05-knowledge/results/session242_congruent_descent_20261110/exception/PROOF.md
session: Session242, 2026-11-10
depends_on:
  - THM-001-redei
  - THM-1370-h-spectrum-omits-7-21-all-n
  - THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak
---

# THM-5107: the seven-cofactor obstruction retains physical normalization

Let (A,a) be a finite rooted tournament with a a nonsink, B a nonempty
finite tournament, H=H(B), and

    Q = sum over X⊆V(B) of H(B[X]) H(B[V(B)\X]),

with H(empty)=1. Form the literal graft by adding B→a and A−a→B, retaining
all internal arcs, and write N for its Hamiltonian-path count and g=gcd(H,Q).
Then N/g is odd and N/g≠7. If |B|≥2, also N/g≥5.

Consequently a proper graft with count 7p, p an odd prime, has g∈{1,7}.
The case p=3 has no output; for p=7 the same conclusion follows directly
from the size bound. No existence of a graft at every admissible target is
asserted.

## Mechanism

Use the inherited physical source state (D,E,U,d), and f=U−D−d+E. The
Session229 proper-graft row and inequalities give

    N = Hd + Q(D−E) + (Q−H)f ≥ H+Q,
    Q≥2H, D,d positive odd, f≥0 even, E≤(D−1)/2.

For |B|≥2, complementary subsets pair the summands defining Q. Every product
is odd and there are an even number of pairs, proving 4|Q. Since g is odd,
h=H/g is odd and q=Q/g is divisible by 4. If N/g=7, the bounds h+q≤7 and
q≥2h force (h,q)=(1,4), the actual pair of the transitive two-vertex gadget.
Replacing B with that gadget on the unchanged source would realize 7.
Alternatively, its equation 7=d+4(D−E)+3f is impossible: D≥3 makes the
right side at least 9; D=1 gives 5+3f with f even. The singleton-gadget
boundary has g=1 and uses the established exclusion of count 7.

This is an instance of the physical normalization principle: if the divided
pair (H/g,Q/g) has an actual gadget witness, replacing the gadget realizes
N/g. A mere integer pair does not authorize that replacement.

Both allowed contents occur in actual representations. The strong gadget with rows
[22,60,8,81,4,29,55] has (H,Q)=(49,812), g=7. Grafted at a vertex of C3,
it gives the strong count 1673=7·239. Its normalized pair (7,116) is
nonphysical, since 7 is forbidden. A separate strong rooted seed followed
by a sink and singleton bridge gives 791=7·113 with gadget content 1.
The source's rooted content 7 is a different invariant from its gadget
content, which is 1.

The nonsink condition is essential. A singleton root and the H=49,Q=812
gadget give N=49 and N/g=7, a literal counterexample without properness.
For an n-vertex gadget, cutting each Hamilton path at every position gives
Q≥(n+1)H and N≥(n+2)H, an additional elementary inverse budget.

[Complete proof and source pins](../../05-knowledge/results/session242_congruent_descent_20261110/exception/PROOF.md),
[hostile independent audit](../../05-knowledge/results/session242_congruent_descent_20261110/graft_audit/REPORT.md).
The exact tests include every labelled gadget through order five and full
permutation counts of the two ten-vertex outputs. The all-order theorem
comes from the proof, not that finite census. Cofactor 21 is not excluded
by this theorem. The remaining prime/7p predecessor-existence problem is open.
