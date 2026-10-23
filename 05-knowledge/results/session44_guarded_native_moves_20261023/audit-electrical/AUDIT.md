# Independent audit: universal native electrical guards

23 October 2026. Verdict: **PASS after one corrected receipt-provenance error**. No mathematical counterexample to the scoped theorems or finite native target result was found. This audit did not modify producer files or perform Git operations.

## 1. A concrete error found and repaired

The first producer version reused the variable `source` in its coupled-separator loop. That loop ended at source 2246, after the literal hostile paths at source 1238 had already been authenticated. The final receipt therefore incorrectly reported `hostile.source=2246`, even though both nested path records correctly said1238.

The producer repaired the enclosing field to1238, added an explicit consistency check for the two path source records, and regenerated the receipt and hashes. The initial self-replays passed despite the error: equality with a self-generated frozen receipt is not an independent semantic check. The initial audit logs preserve that diagnostic history. The corrected final receipt is

    5f7dea23f728f685211dd5eac1dff0a389a12461d84946ec9f13b2d7f86faba2

The final producer script is

    c99b6a67f8eb209f7d9b9920c38f21dae38edc80d27bb1a4e01147129b202897

## 2. Mathematical proof review

- The annihilator identity `(ker M)^perp=im M^T` is valid over every finite ring `Z/mZ`, including composite moduli. The Smith-coordinate scalar proof does not substitute field rank for module structure.
- One root grounding is harmless for zero-sum charges because full harmonic vectors differ from grounded ones by constants. Additional constraints require the extra transpose-image term. The doubled-edge pinned hostile correctly distinguishes the two situations.
- The quotient `K/mK` and root-reduced Smith signature require the owner-marked left transformation, not invariant factors alone. The quoted order-eight decomposition is inherited from the pinned exact 2-primary calculation, not recomputed as a new Smith theorem here. Its3369-bit size and 3363 minimal module generators are arithmetically consistent.
- For zero-row-sum perturbations, `DeltaL=BL mod m` is equivalent to keeping every original harmonic vector harmonic. This is inclusion, not equality of old and new banks. The audit retains `DeltaL=-L` as a hostile: the new zero matrix admits additional harmonic vectors. Producer wording states the correct inclusion claim.
- The all-modulus vertex-alias theorem is sound for a connected positive-integer-conductance multigraph. If `Lf=e_u-e_v` with integral f, every threshold cut between f(v) and f(u) carries a total integer current of one, so exactly one unit edge crosses it. The nonzero gradient current is an acyclic integral unit flow and therefore a single source-to-sink path of such bridges. Conversely each unit bridge supplies its component-indicator image witness. Parallel edges of total conductance greater than one are not literal bridges. The trivial u=v case uses an empty bridge path.
- A single fixed modulus cannot replace all moduli: the unit triangle has only the constant harmonic bank modulo 8 but has a nontrivial bank modulo 9. This is included as an exact hostile.

## 3. Why the finite modulo-eight bank is complete for vertex comparison

The producer does not need to prove that its listed vectors generate the entire `Z/8Z` harmonic module. Every listed vector is checked against the full integer Laplacian and the actual root. The signatures separate every vertex pair except three pairs. Each remaining pair has a direct integer image certificate, so no harmonic character at any modulus can separate it.

This proves the exact equivalence relation on vertex charges, not a general-purpose coordinate representation of every charge. In particular it does not turn the 3525 abstract vertex labels into an additive group in which aggregate charges may be subtracted.

The independent audit uses a different binary construction. It reduces each vertex unit charge modulo the full binary Laplacian row space, rather than constructing the producer's grounded harmonic basis. It independently finds rank164 and the 3521/3524/3525 dyadic signature class counts after authenticating the higher lifts. The only final aliases are

    {0,3527}, {53,3492}, {109,3523}.

The independent code checks their integer Laplacian images and the absence of parent self-loops. Thus none can occur as two targets at one source.

## 4. Complete native universe and what remains inherited

The audit independently reads every one of the 2^21 normalized presentation masks from the frozen owner map. It retains both actual source and complement target, reconstructs all 438216 unordered pair rows with their exact multiplicities, and checks each source neighborhood against the edge export.

All 3528 source target sets are injective under the universal modulo-eight signature. The same replay confirms 3518 sources with collisions under the single original character, and 777396 false vertex-pair equalities relative to the complete character equivalence.

This consumes the previously authenticated tournament owner map. It does not redo isomorphism classification, claim to enumerate every separately labelled Hamiltonian path, or reconstruct target ownership from incomplete clock information. The equality of all characters identifies a target only after the native endpoints have been authenticated. The result supplies neither path connectivity nor an all-order native wiring classification.

## 5. Independent finite controls

`independent.py` uses only the Python standard library, exact integers and rational fractions. It imports no producer code. Its additional controls include:

- All 646 connected graphs on 2,3,4 labelled vertices with edge weights in{0,1,2}, checking3806 distinct vertex pairs. Solving the grounded integer-image equation over the rationals agrees exactly with bridge-path equivalence.
-512 zero-row-sum two-vertex perturbations modulo 8, including nonunit edge weights and nonsymmetric perturbations, checking the row-image criterion against direct evaluation on the grounded harmonic bank.
- A fixed-modulus cycle hostile, a strict-bank-expansion hostile, and direct checks of all higher native harmonic lifts.

Final producer normal and optimized replays are byte-identical: 887257 checks. Independent normal and optimized replays are byte-identical: 11404 checks. The independent normal output, including every retained class partition, is saved in `independent.normal.json`.

Both checkers were also tested in a layout omitting the producer's duplicated inputs, reading the pinned canonical inputs from the supplied repository root. Outputs agree exactly with the local-input runs. The independent invocation supports

    python3 -B independent.py --producer PATH_TO_PRODUCER --repo REPOSITORY_ROOT

and the same command under `python3 -O -B`. No optimization-sensitive assertion is used.

## 6. Scope and publication handoff

Publish this audit, independent checker, final receipts, and manifest without copying inherited atlas inputs. `AUDITED_ARTIFACTS.json` records exact producer artifact hashes. Initial logs are historical records of the corrected provenance bug and are not the current producer receipt. Target rigidity, all-character comparison, original one-character preservation, per-slot preservation, owned aggregate preservation, shared pins, actual kernel lifts, and full wiring remain distinct claims.
