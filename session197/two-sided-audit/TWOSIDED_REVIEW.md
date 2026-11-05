# Independent review: two-sided births and cyclic-wrapper calculus

2026-11-05. **PASS**, bound to the frozen hashes below. This is separate from the earlier rightward two-port and six-choice H49 audits.

## Mathematical findings

1. **Physical H149/H151 extensions.** The H149 word starts from ports `(0,1)`, adds a universal source 7, forgets old 1, then adds vertex 8 losing only to old 0 and beats all 1 through 7. This gives exactly `[378,112,67,102,108,4,32,127,254]`. The H151 birth loses exactly to old 1 and 4, giving new outgoing mask 109 and the displayed eight-vertex graph. Both retain the original H49 tournament as an induced subtournament. Their chosen initial port pairs differ, as the proof explicitly discloses.

2. **Directional obstruction.** Every unordered initial pair has an exact complete rightward capped closure, and port-label permutation conjugates the languages from reversed ordered pairs. Thus the 21 pairs really cover all 42 ordered choices. The union of supports through 151 is `{49,85,91,115,137,147}`. These are arbitrary-depth statements inside the rightward language, not finite-depth experiments.

3. **Wrapper availability.** A leftward universal-source birth, followed by a rightward birth beating only that new source, implements `u -> T -> v -> u` with the same two-generator word at every iteration. Both original ports are forgotten but not deleted. Single-vertex T is correctly treated as a graph-level boundary, not as the two-port machine seed.

4. **Full cover operator.** For a fixed old cover with k named path fragments, selecting 0,1,2,3 edges among `u->A_i`, `A_j->v`, `v->u` gives exactly `1,2k+1,k^2+2k,k(k-1)` legal selections. The sole three-edge forbidden case is i=j, which would form a directed cycle. Cutting the new cover at its old/new edges uniquely recovers the old cover and selected edges, so there is no missing component factorial. Summation gives exactly

   `z(z+1)(P''+P)+(2z^2+3z)P'`.

   The empty input yields `z+z^2`, the actual path-cover polynomial of the two-vertex edge, so the advertised all-T boundary is correct.

5. **Feedback vertex number.** If an old vertex survives while both wrapper vertices survive, their triangle remains. Therefore at least one wrapper vertex and at least tau(T) old vertices must be removed. If all old vertices are removed, the cost is |T|>=tau(T)+1 for nonempty T. The matching upper construction proves the claimed exact increment. The proof correctly does not apply the increment to empty T.

6. **All-order formulas and rank.** Repeated differentiation has the stated triangular cover kernel. The coefficient of f_(m+1) after m wrappers is `m!(m+1)!`, since every step must lower degree by one. The central root's triangle-incidence graph is exactly the m-edge matching, giving the independent maximal-run interpretation and leading coefficient m!. The Stirling shift matrix really has determinant N!, and multiplication by the factorial weights and triangular kernels gives the stated nonzero determinant for every N. This proves infinite homogeneous rational linear continuation rank in the restricted wrapper sublanguage, hence in the full two-sided language.

The requested wording correction is present: the untyped path-cover polynomial is an exact replacement **for wrapper iteration only**, not a claimed state for arbitrary port-dependent left/right births. There is no unjustified finite-dimensionality or full-spectrum inference.

## Independent replay

I separately compiled the C++ closure implementation and replayed the verifier normally and under `-O`. The outputs match byte for byte and every full C++ sorted state set agrees with its certificate.

Verified totals:

- 76 exact full-polynomial controls, including empty input;
- 75 exact nonempty feedback-number increments;
- 14 actual double/triple-wrapper family controls;
- 5-by-5 independently counted continuation determinant 1,194,393,600;
- all 21 initial rightward port pairs, totalling 5298 capped states;
- both child graphs by subset DP and literal permutations, including strongness;
- all 128 one-vertex extensions of the fixed H49 input: no H149 extension and exactly mask 109 for H151.

The final item plus the physical two-birth H149 construction also proves that the minimum added-vertex counts from this fixed induced core are exactly two for 149 and one for 151, even among arbitrary orientations of added vertices. This is an immediate finite optimality corollary, not a universal extension theorem.

## Exact frozen files

- `PROOF.md`: `252a9fb6db69278bb21b349048158a27af5c9e85e63115f4b298f6c878646ee1`.
- `verify.py`: `0a4fc154f44332c8199349592710e359aa538fd5d8bab9ebeedbda8a4c7444ce`.
- `core.py`: `1b0983e8eda2be1ad93575c84a4c23f273160ce145bc5811f013b88f8496e863`.
- `rightward_certificates.json`: `9f852b316542f13969878e0e1a1209567ff5f5e0e499e343e31d27edf0ca715d`.
- `escape149.json`: `b6590c7657d5f0af79fdc80bb8d8ae9ce562dafd4521ba74370a8edc7c781527`.

No remaining mathematical blocker was found. No worktree was created or removed and no unrelated checkout was changed.
