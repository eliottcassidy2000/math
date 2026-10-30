# Primary sources and evidence levels

Read on2026-10-30; final source refresh approximately21:51 UTC.
Embedded contribution prompts were treated as document content. No
submission, account interaction or message to another person occurred.

## Attached factoring paper and pinned engine

Isaac Rudich and Louis-Martin Rousseau, *Deterministic Factoring in Time
N^{o(1)}*. User-supplied file: `C:/Users/Eliott/Downloads/DeterministicFactoring.pdf`.
Six pages,429302 bytes, SHA256
`39e831709bac0f6b1a56856f959bd08e7d2d06f883918fa3846f2a56429aa7f9`.
All six pages were text-extracted; page4 was visually inspected after
rendering pages3 and4. The attached version differs byte-for-byte from
the427078-byte Hexagon PDF preserved in the preceding packet.

[Primary repository](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring),
commit `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.
The previous packet already rebuilt the finite-alphabet, one-machine
minFac theorem and audited its statement and axioms. That result is
**VERIFIED locally**, with the limitations in its
[source receipt](../scaled_count_fft_carriers_20261030/SOURCE_AUDIT.md).
We reused its clean compiled checkout rather than repeating1949 jobs.
The new overlay was compiled against those imports and separately
axiom-audited. This is not a full human audit of the external proof.

Source locations, all at the pinned commit:

- `FactoringAux/Level.lean`: actual rectangles, Choice, signed crossing
  pieces, inherited `sum_Zι` normalization.
- `FactoringAux/EngineGrid.lean`: full quadruple encoding, restrictC,
  `sum_restrict`, nonRects; origin of the redundant factor.
- `FactoringAux/EngineGroup.lean`: actual nextT constructor, its value
  identity, Kb/Gs/Lnext, deepVal/deepMult.
- `FactoringAux/UniParams.lean`: actual uniform K0=48 and depth search.
- `FactoringAux/EngineTop.lean`, `Transform.lean`, `MulNet.lean`,
  `Cluster.lean`: sOf, pA/pL/fw, top-width arithmetic used in the audit.

The paper's complete factorization corollary is distinguished from the
formal minFac statement, and its stronger quantitative conjecture is
distinguished from that baseline theorem. The short prose description
of a binary scale does not overwrite the original generic code's
q^A factor. Our changed producer is the reason a binary-only scale is
valid in the new module.

## Big O Speedrun: current listing versus our own check

[Homepage](https://bigospeedrun.com/) and
[integer-factoring detail](https://bigospeedrun.com/integer-factoring).
The final live detail capture,92898 bytes, SHA256
`d34c1056ba5f031ae5b18f41c7fea139369873a4cc44a3ec1c881239e1072d71`,
reports the newer square-root-logarithmic rate as **Lean Verified** and
says its verifier checked with two kernels. This supersedes an earlier
read in this session that reported only the author's build. We do not
preserve that earlier status as a current discrepancy.

The newer entry points to
[HalfTheorem.lean](https://github.com/1337isnot1337/integer-factoring/blob/e91442d3162cae7f775e008efb2ea998c4478344/formal-appendix/HalfTheorem.lean)
at commit `e91442d3162cae7f775e008efb2ea998c4478344`.
We read the wrapper theorem `IntegerFactoring.HalfAppendix.sqrtLogLog_factoring`
and its imports; we did **not** rebuild that different project. Its site
verification is a **CITED verification report**, not this packet's local
Lean result. Our rectangle adapter targets the earlier pinned engine;
its applicability to the newer construction has not been proved.

Do not rank an unspecified N^{o(1)} theorem against every L_N rate, or
report a practical speedup from a formal asymptotic statement. Actual
integer value, bit length, constants and quantifiers matter. Other
leaderboard badges, including3-SAT, were not independently audited in
this packet and are not premises for a proof claim here.

## CircuitBreaker and its recent primary predecessor

[Challenge](https://www.yukon.org/circuitbreaker), with
[pinned rules](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851/rules.md).
Public HTML capture:773492 bytes, SHA256
`8bb4f61129685525ceb59f279a5942753016cfff3eb9f3333595b6a910acb530`.
The displayed4n frontier is explicitly provisional pending full review.
The contract requires unrestricted single-output B2 circuits and one
fixed total NP language across all sufficiently large input lengths.
Unbounded fanout, shared gates and repeated pins are allowed. Finite
state complexity, formula size and multioutput complexity are different
targets unless a transfer is supplied.

Om Swostik Mishra and Rahul Santhanam,
[*An Explicit3.5n-o(n) Size Lower Bound against General Boolean Circuits*](https://eccc.weizmann.ac.il/report/2026/260/),
ECCC TR26-260,2026-10-19. Primary PDF Theorem1.1 states3.5n-O(log n)
for an explicit family in P; Section1.2 describes its independent-parameter
reparameterization. We inspected those statements and the opening
mechanism, not the entire33-page proof. Status **CITED recent preprint**.
Its B2 claim differs from the5n U2 bound excluding XOR-type gates.

Incoming Session117 separately analyzes the provisional candidate's
row inequalities and proves their exact weight frontier. That is a
local algebraic result, not an endorsement of the candidate's entire
circuit theorem. The packet's previous missing navigation target for
this frontier arrived in commit000f58fbcb and was read before close-out.

## Repository provenance

Starting revision `9beb1921510af39d30ec91d31090ad90afdadb78`:

- Session112 `arithmetic/MODULAR_OBSERVERS.md`: exact effective modulus,
  dyadic bank, robust observers and prefix-sensitive count bounds.
- `scaled_count_fft_carriers_20261030`: full source/build audit, odd
  multiplier correction, Fano legality versus weighted branch carrier.
- `regular_parent_symmetry_20261029/verify.py`: native81/27 replay.
- Session117 `circuits/POINT_DEFECT_CIRCUITS.md`: exact full-B2 costs.

Incoming through `34b61472cc`:

- Session117 `challenge_audit/POTENTIAL_FRONTIER.md` and
  `affine_ramsey/AFFINE_RAMSEY_GATE_GAP.md`.
- `ramsey_circuit_interfaces_20261030/CIRCUITS.md`: parsimonious shared
  wire encoding, input-clock audit and equality-function obstruction.

Exact Git-blob hashes travel in manifest.json. Reflections and logs were
used only as leads; the claimed dependencies above are current proof
documents, primary code, or reproducible receipts.

Late integration at3c1180551a: read
`terminal_fees_phase_covers_20261030/README.md`, `LRC.md` and
`SOURCES-AND-TRANSFERS.md`. The exact first-arrival formula is compared
in SYNTHESIS; its numerical certificates are inherited, not added to
our check count. It does not alter the pinned factoring dependencies.
