# Target-preserving returns, growing defect, and triangle metagraphs

**Completed research packet, 11 November 2026. PROVED IN NOTES / FINITE-EXACT.
Chan–Pak remains OPEN.** Started at `6512e16e0729`; first checkpoint
`8939375fbf`. Proof statuses describe the internal reviews recorded below,
not a literature-priority claim or an external referee certification.

[Synthesis](SYNTHESIS.md) · [Inheritance and connection ledger](INHERITANCE.md) ·
[Concept board](BOARD.md)

The inherited (+4,+2) path now has two successive **uniform, absolute-count
preserving returns**. Its base sequence reaches23→73→221→665→1997.
The failure of readiness at1997 is overcome by an actual indegree-one,
defect3 source, with finite named returns to both5993 and5995 and further
continuations. A uniform law permitting unlimited repetitions is still missing.

A separate complete theorem classifies every extension-gcd layer of every
ordered chain of triangles and singletons. This gives exact all-size
enumerators, lifetime layers, and every birth/deletion or arc flip that
stays in that family. It exposes information lost by the usual rooted
seven-coordinate observer.

| Question | Proof and evidence |
|---|---|
| Uniform ready returns at the actual (+4,+2) target | [TARGET_RETURN](reset/TARGET_RETURN.md), [independent audit](berggren/RESET_AUDIT.md) |
| Indegree-one defect3 escape and the5995 repair | [DEFECT_THREE_ESCAPE](reset/DEFECT_THREE_ESCAPE.md), [audit](berggren/DEFECT_ESCAPE_AUDIT.md) |
| The5993 branch, next ready return, and larger-core inverse | [H15_MENU](reset/H15_MENU.md), [audit](berggren/H15_AUDIT.md) |
| Alternative alternating path and failure of the growing cap | [ALTERNATING_AND_CAP](reset/ALTERNATING_AND_CAP.md) |
| Berggren decoration and actual wide-prefix returns | [Compiler](berggren/PROOF.md), [WIDE_RETURN](berggren/WIDE_RETURN.md), [size boundary](fibres/WIDE_REWRITE_SIZE.md) |
| All-k extension content for every triangle chain | [Half-contact proof](arithmetic/TRIANGLE_ALL_K_PROOF.md), [independent audit](berggren/ALL_K_AUDIT.md) |
| Complete internal family metagraph and enumeration | [TRIANGLE_METAGRAPH](fibres/TRIANGLE_METAGRAPH.md), [gap transport](fibres/GAP_TRANSPORT.md) |
| Same rooted observations, different extension obstructions | [ROOT_OBSERVER_LOSS](fibres/ROOT_OBSERVER_LOSS.md) |
| Complete decorated H9 fibres and branch counts | [H9_FIBRES](fibres/H9_FIBRES.md), [extension content](fibres/EXTENSION_CONTENT.md) |
| First25 primes,30000 translation, and square-prime pairs | [Arithmetic proof](arithmetic/PROOF.md), [compensation boundary](arithmetic/COMPENSATION.md) |

From the repository root:

```text
python -B 05-knowledge/results/berggren_target_return_20261111/reproduce_all.py
python -B 05-knowledge/results/berggren_target_return_20261111/reproduce_all.py --optimized
python -B 05-knowledge/results/berggren_target_return_20261111/seal_packet.py --check
```

Python and a C++ compiler (`g++`) are required. The replay receipts state the
finite universes; infinite claims use the linked proofs. Source pins, link
checks and a normalized-text publication manifest accompany the packet.
Discovery-only `reset/probe_ready.py` is retained as provenance, outside the
proved-claim replay. A mistaken contiguous-triangle argument in the first
checkpoint was repaired and logged in [MISTAKES](../../../01-canon/MISTAKES.md).
