# Independent review: complete low-cycle two-port compiler

**PASS. No correction requested.** This review binds the exact source manifest reproduced below. The original packet was not edited; all replays took place in a separate copy.

## Mathematical audit

The start-one compiler is a valid induction on the initial strong component after deleting a selected root. With a nonempty tail, the old-root count law excludes every other old start-one vertex except the entry. With an empty tail, equality in s_x−p_(x,a)+D_(x,a)=1 forces every old x-start Hamilton path to close to the unique old cycle, hence s_x=1. Therefore the requested second old vertex is already legitimately exposed by the recursive call, never re-exposed after forgetting. Tail births need only the retained new root; a requested tail vertex can stay as the second port.

The start-two compiler uses exactly the two entry cases in the complete structural decomposition. Single-entry recursion decreases the order; the two-entry case uses the proved all-pair start-one compiler. The regular-five exception is explicitly constructed. Normalizing only the initial two-vertex edge preserves exposed sets, with every later mask/drop re-indexed, so it is not an unauthorized port swap. The theorem promises final exposed sets, not an independently chosen port order.

The cycle readout is a bijection obtained by deleting the newest vertex from a Hamilton cycle. Its start must be selected and its end unselected. This proves the stated four nonnegative readouts and independence from the forgotten port. A prefix with more than two cycles is strong, and strong-induced cycle monotonicity compares it directly with any final strong target. Intermediate weak graphs do not invalidate that argument. The one-vertex insertion and two-vertex ear cases give the required monotonicity; no false cycle monotonicity for arbitrary weak intermediates is assumed.

Coordinate clipping commutes with every nonnegative integral transfer and readout. Thus the raw cap 3 closure and quotient recognize arbitrary-depth cycle admissibility, not a finite word-length sample. The bounded H search records accepted incoming edges before deduplication, so histories with different current cycle counts are not silently lost. The H cap is safe by induced Hamilton-path monotonicity. At retained states all coordinates contributing to H and its cycle subreadouts are exact. The finite closure is therefore an exact target-bounded certificate.

The all-order graph compiler, paired with cycle readout and safe pruning, proves the full low-cycle inverse-language criterion. The implication to the {7,21} spectrum remains conditional on unproved all-target low-cycle coverage. The packet states this distinction correctly.

## Computation and independent controls

- Copied source normal/-O compiler replays both match 3bfcb8b2b0c4eb5deeb4a1837b47504488e9b5a576509e10ec1ae9e39b7c37fe.
- Copied full support normal/-O replays, including the 3,455,583-state cap 501 closure, both match 9f5d89f5caf32fe6ccfe06de3c295ad1479b7cbc4625101a7cfe053dbe0f0d3f.
- I independently computed root-start counts and cycle counts by a forward subset DP on the complete uncapped unlabelled bank through 7. On all 396 strong graphs it checks 90 single start-one marks, 78 prescribed start-one pairs, 84 start-two roots, and the regular-five graph. A separately written literal birth replay checks all graph edges and retained original vertex identities in every compiled word.
- A separate pair-distinguishability calculation checks all 741 distinct pairs of the 39 quotient classes, verifies every raw-state merge, and finds exactly 38 coaccessible classes and one rejecting class. This is independent of the source's partition-refinement algorithm.
- The unchanged 27-coordinate transfer and its saturation justification were previously independently reviewed. I checked the serialized transition regeneration, source support algorithm, complete-processing condition, and uncapped witness replay in this packet.

Independent checks are in `audits/lowcycle-two-port/independent.py` and `independent.json`. These computational controls supplement, rather than replace, the structural induction and clipping proof.

## Exact frozen source manifest

- `PROOF.md`: `c4a16105ebc9f9aa9cb3cd61887ec4c437db58c361b2469af630899d536c5c73`.
- `compiler.py`: `71d99b1e86f8af3503c2461f0563b774ef699c0c78e95286f03196befb994c8f`.
- `verify.py`: `808284da62ded941abd0c920292c827bd25309286c22e6444f0f710d63c72cc3`.
- `VALIDATION.json`: `3bfcb8b2b0c4eb5deeb4a1837b47504488e9b5a576509e10ec1ae9e39b7c37fe`.
- `VALIDATION_OPTIMIZED.json`: `3bfcb8b2b0c4eb5deeb4a1837b47504488e9b5a576509e10ec1ae9e39b7c37fe`.
- `twoport_core.py`: `733155121ac7a5936a38a0ec1f3a602759f868c1e4e547fd4b6f89bced25f9cb`.
- `automaton.inc`: `de75a17dfb95fc1cbdfd789a969ea0ca03b1f573b9600e51a6bd2c9c4d914b10`.
- `support.cpp`: `0f97d243c87a7fefe4cc6ffd00f6ce86a06a755ec7cf74093d54014aea33f45c`.
- `support173.json`: `4700d5056b3d70dbefe1eb8e8a59fe9c240b55eb1c282b447516d5fc04583e57`.
- `support501.json`: `0fd9a15f93f1f5f6964155dc90bb68d763306e6191b72f46bfac76bc3b848705`.
- `control.py`: `d984bc6fce48cf011abf6d011250833b5a4b26b32f89790bca96119a72bf2498`.
- `control.json`: `0bdac3b9c0db42abf3c6c3c13ca547da2e6b5241c7e86a90cab07b870c699833`.
- `verify_support.py`: `20e79ff1b4520e6bfd89c879fb512073eaa49ab30049cea34e30aba7b45ebbba`.
- `SUPPORT_VALIDATION.json`: `9f5d89f5caf32fe6ccfe06de3c295ad1479b7cbc4625101a7cfe053dbe0f0d3f`.
- `SUPPORT_VALIDATION_OPTIMIZED.json`: `9f5d89f5caf32fe6ccfe06de3c295ad1479b7cbc4625101a7cfe053dbe0f0d3f`.
