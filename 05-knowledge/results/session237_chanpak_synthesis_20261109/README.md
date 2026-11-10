# Session237: four bridges between the Chan–Pak approaches

**PROVED, independently audited partial results. Chan–Pak completeness remains OPEN.**
This session recovered older graph, tiling, orbit, endpoint, rooted-tail and
modular approaches and made their connections explicit. It establishes no new
unrestricted finite coverage cutoff and no proof of all odd counts except7,21.

## 1. One arithmetic statistic joins the endpoint and tail pictures

[THM-5093](../../../01-canon/theorems/THM-5093-endpoint-half-bound-equality-and-binary-pair-energy.md)
identifies the inherited binary pair-energy Q with exact endpoint slack in
an explicit tournament family:

    E=m-Q(m),   H=1+2m-Q(m),   H-2E=1+Q(m).

The same Q already drives fixed-port tail counts and modular pumping. This
is a physical construction with an actual endpoint, not a matching formula.
The theorem also classifies equality in the general endpoint half bound:
for a positive nonsink endpoint, H=2E+1 occurs only in the reversed-extreme-arc
family, with a possible transitive prefix. Thus the smallest possible full
count for endpoint N equals2N+1 exactly when N is a power of2.

This combines THM-4876's bound and finite equality observation, Session231's
join-extremal face, and THM-5001's endpoint universality. It explains precisely
why universal endpoint counts do not automatically produce universal full H.

[Proof](endpoint_slack/PROOF.md) · [Independent audit](endpoint_slack_audit/AUDIT.md).

## 2. Exact modular memory, not integer coverage

[THM-5091](../../../01-canon/theorems/THM-5091-minimal-modular-tail-memory-and-joint-seed-images.md)
proves that the normalized full-count observer (H(m)-1)/2 modulo M has exactly
2M² minimal deterministic MSB-first states for every M>=2. Every state is
reachable by inherited block pumping; every pair is distinguished by at most
one further input bit.

It also gives the exact simultaneous residues of any finite seed list sharing
one tail: two translates of a two-generator module. One explicit pair of
physical seeds lies on two affine lines; another fills every residue pair
at each fixed parameter parity. These are parameter-synchronized statements,
not independence claims about banks using unrelated parameters.

[Proof](arithmetic_states/PROOF.md) · [Independent audit](arithmetic_states_audit/AUDIT.md).
The [endpoint-observer supplement](arithmetic_states/ENDPOINT_OBSERVER.md)
proves that E modulo M needs2M² states for odd M>=3 but only M²/2 for even M
(one state at M=1). This is fourfold smaller at even M than the normalized
full-count observer above. Actual H modulo M is a different observer; its
state count must not be confused with that of (H-1)/2 modulo M.
[Independent extension audit](endpoint_slack_audit/AUDIT.md).

## 3. Root changes escape traps that the old interface cannot see

[THM-5092](../../../01-canon/theorems/THM-5092-four-test-tail-content-and-invisible-reroot-escape.md)
reduces the gcd of infinitely many fixed-port tail counts to four words:
empty, B, SB, BB. This arithmetic content can exceed one even in a strong,
fully asymmetric tournament. In the explicit H=9 example, changing its root
then adding one bridge gives19; taking its converse then a bridge gives29.

More strongly, two **six-vertex** rooted tournaments have identical
(h,R,e,d)=(9,9,3,3), and therefore identical counts after every fixed-port
word. Only one can escape to content one by rerooting; their converse states
also differ. Even adding the vertex count to the interface does not repair
the loss. Root-resolved endpoint/deletion/response data are necessary before
using a compressed state in a grammar that changes roots.

[Proof](symmetry_grafts/PROOF.md) · [Independent audit](symmetry_grafts_audit/AUDIT.md).

## 4. Arbitrary graph carriers exist, but extra cycles carry positive cost

[THM-5090](../../../01-canon/theorems/THM-5090-universal-induced-odd-cycle-carriers-and-positive-completion-tax.md)
embeds every finite graph G induced in a tournament's odd-cycle conflict graph.
The construction uses at most |E(G)|+6|V(G)| vertices and also installs a
private triangle at each selected cycle. Every completion satisfies

    I(Omega(T),x) >=_coeff (1+x)^m I(G,x/(1+x)) + m*x.

Thus no nontrivial hereditary property can characterize the full cycle
carriers. The obstacle to importing graph/poset/tiling ideas is not the
ability to install a selected intersection pattern. It is control of all
extra cycles and their positive packing contributions.

The proof connects to the exact tiling dictionary and audited even-graph
domination theorem, while preserving their distinct statistics. In particular,
even graphs attain21 and tournaments do not, so an unrestricted count-preserving
map from that graph class is impossible.

[Proof and connection ledger](cycle_tilings/REPORT.md) · [Independent audit](cycle_tilings_audit/AUDIT.md).

## What these bridges suggest next

- Combine root changes with a root-resolved response state. Content-one is a
  necessary large-prime filter for a finite fixed-port bank, not sufficient
  target reachability.
- Use graph/tiling carriers only when their extra-cycle contributions close
  in a controlled transfer state. Induced universality alone gives no selector.
- Treat residues and endpoint values as exact tools with known loss, then
  attack integer occupancy and baseline placement directly. Full modular
  projections, long consecutive windows and primitive content do not close it.

The novel feature of this session is the linkage: the same pair-energy is
now an endpoint defect and an arithmetic memory coordinate, while the
root-change and completion examples identify two concrete kinds of missing
state. No external literature-priority claim is made.

## Verification and publication

[Reproduction and exact scopes](REPRODUCE.md) · [Synthesis and boundaries](SYNTHESIS.md).
The source packet, independent proofs and independent algorithms are retained.
The first verified nonforce checkpoint is
`ebedba7d7d8207b922e7df4ee81f21cc976b20ee`.
Final compact replay passed all ten producer/independent programs in normal and
optimized modes, with both deliberate-failure controls for each. The separately
replayed full202013-state symmetry census has identical preserved outputs.
See REPLAY.json and SHA256SUMS.

**Session end:** research, independent audits and verification are complete.
The published results are partial theorems; Chan–Pak completeness remains OPEN.
