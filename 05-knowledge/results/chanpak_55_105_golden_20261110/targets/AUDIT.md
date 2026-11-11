# Independent audit of the common-source obstruction

**Verdict: PASS.** P1's all-order padding reduction and the P3 obstruction
at H105 are valid with their stated operation and inheritance scope. No
internal-padding loophole was found. The unrestricted spectrum conjecture
and the weaker induction allowing different actual sources remain open.

## Proof audit

For a child T=S+v with gcd(H(S),H(T))=1, let C be the SCC containing v.
The other SCCs contribute a common factor to both path counts. That factor
is1, so every other SCC is a singleton. Since H(T)>H(S), C is not a
singleton. Tournament SCCs have a unique total order; hence

    T = TT_l -> C -> TT_r,
    S = TT_l -> (C-v) -> TT_r.

This gives a particularly direct proof that stripping these vertices
cannot remove an internal singleton gap of K(S): every stripped vertex
is before **all** of C-v or after **all** of C-v. Every nontrivial SCC
of S lies inside C-v. A singleton of S between the first and last such
SCCs therefore lies inside C-v as well. Normalizing C-v removes at most
more exterior singletons and leaves exactly K(S).

Thus every branch from one S must appear in the marked-child inventory
with the same middle core, including every internal singleton. Conversely,
the coordinatewise maximum of finitely many exterior padding pairs works:
add universal sources before each entire child and universal sinks after
it. The old sources become identical after the recorded core isomorphisms,
and each child count is unchanged. The cuts to newly added exterior
vertices are determined, so no unproved module-contiguity assumption is
used. The theorem concerns alternative one-vertex branches; it does not
assert simultaneous insertion of all their new vertices.

The hypotheses matter: h>1 ensures a middle core exists, N_i>h ensures
the distinguished component is nontrivial, and coprimality removes all
nontrivial exterior factors. The core intersection alone proves no
claim about preservation of a previously specified joint-port observer,
F49/F343 phase or automorphism-order data. PROOF.md states this boundary.

For P3, every target is below375. The inherited strong inventory is
complete through order11, and the canonical bound in
[THM-4744, interval chains and strong inventories](../../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md)
gives minimum375 at order12 and a nondecreasing minimum thereafter.
Consequently every reduced strong child of count317 or319 is in the
supplied inventory. This is an all-order conclusion conditional only on
the explicitly inherited audited census and cited bound, not an inference
from a fixed source-order experiment.

The SCC generating functions also pass the type audit. A strong55 core
contributes four order7 choices and two exterior singleton gaps. Reducible55
has the two orders of the unique strong5 and strong11 components and three
singleton gaps. For105, the only reducible nontrivial factorization allowed
by the inherited omissions is3·35; its four ordered skeletons have order10
and three gaps. Internal padding is included by the third denominator
factor and is never silently stripped. Converse classes remain distinct.

## Independent exact computation

The auditor imports none of targets/check.py. It uses recursive predecessor
path counting, independent forward reachability for SCCs, and a different
complete isomorphism canonicalizer: sort vertices by out-degree and minimize
over **all** permutations inside the degree blocks. Since isomorphisms
preserve out-degree, this is a complete canonical form, not a degree-sequence
heuristic or a Hamilton-word canonicalization.

The audit rereads the same explicitly inherited bank, whose SHA256 is
`2c66c872b7af50346624a367150a6fb094af916fa1397ef707f342a897a73bd3`.
It checks all300 selected strong child classes and all2584 vertex deletions,
then compares every source-core key and **every recorded padding pair**.
The resulting branch-core counts are:

| Branch | Middle cores |
|---|---:|
| 55 -> 111 | 2 |
| 55 -> 167 | 2 |
| 55 -> 169 | 2 |
| 105 -> 317 | 3 |
| 105 -> 319 | 1 |
| 35 -> 107 | 0 |
| 35 -> 109 | 0 |

The three H55 sets coincide. The two H105 sets are disjoint. Both displayed
simultaneous H55 source witnesses are independently replayed on their
literal labels. A separate hostile constructs C3 -> singleton -> Q4:
normalization retains all eight vertices and its internal singleton;
adding two exterior sources and one exterior sink strips exactly those
three new vertices and leaves the same middle core.

Reproduction from the repository root:

    python3 -I -B 05-knowledge/results/chanpak_55_105_golden_20261110/golden/audit_targets.py
    python3 -I -B -O 05-knowledge/results/chanpak_55_105_golden_20261110/golden/audit_targets.py --output TARGET-AUDIT-optimized.json

The normal and optimized JSON receipts are byte-identical:
[normal](../golden/TARGET-AUDIT.json),
[optimized](../golden/TARGET-AUDIT-optimized.json).
The original bank enumeration is inherited, not rerun by this audit.

No correction to P1–P4 is required. In particular, the result refutes a
common-source branching rule at105 while leaving the different-source
branch rule logically intact. Any future reset theorem must genuinely
change the middle core between the displayed branch families, and must
separately certify the port and phase data it advertises.
