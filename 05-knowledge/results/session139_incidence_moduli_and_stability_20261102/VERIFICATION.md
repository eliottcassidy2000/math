# Session139 verification and checkpoint lineage

**Current accepted scope:** exact weighted stability at n=3..9 and every
n>=526; exact all-modulus context function/clock/padding theorems; finite-exact
complete prime context states at3,5,7,11; the weighted n5 spectral boundary.
**10..525 remains OPEN.** No RH/GRH or integer-spectrum closure is claimed.

| Checkpoint | Scope | Evidence |
|---|---|---|
| 1, published `400dcd629598c30762eb9adb8bb85bea581698a7` | n>=748, weighted n5 spectral loss, corrected historical recovery | [Frozen scope](CHECKPOINT1.md), [original replay](REPLAY_CHECKPOINT1.json), [original source integrity](SOURCE_INTEGRITY_CHECKPOINT1.json); 42 published owned blobs verified |
| 2 | n>=526; all-modulus functions/clocks/padding; exact finite prime quotients and actual5/11 projectors | [Frozen scope](CHECKPOINT2.md), new source pins and independent audits, [fresh replay](REPLAY_CHECKPOINT2.json) |

All earlier scientific lane files, audit reports and replay receipts remain
byte-unchanged. Root README/SYNTHESIS/BOARD/VERIFICATION and SHA256SUMS are
maintained summaries/manifest. `CHECKPOINT1_SHA256SUMS` preserves the original
whole-packet manifest **at checkpoint-one's commit**; it is not a command to
validate later rewritten summaries against old hashes. The original
`PUBLIC_FILES.json` is likewise the historical checkpoint-one selector.
The root `SHA256SUMS` validates the cumulative current packet.

## Second checkpoint: light-mass dual closure at n>=526

[Proof](stability/PROOF_N526.md) · [Independent endpoint/full audit](independent_threshold526/AUDIT.md) · [Second structural audit](independent/N526_STRUCTURAL_AUDIT.md)

**ACCEPTED, including equality.** The proof retains actual light mass U,
uses a nonempty heavy matching, and pays the light-deficit budget with an
exact nonnegative quadratic dual. Endpoint-convexity and monotonicity—not
sampled orders—supply the uniform geometric tail bound. All divisions,
strict signs, parities and matching-coefficient directions were reviewed.
The final exact slack is `81/8416000`.

The producer and independent integer/Fraction certificates match their pinned
receipts in normal and optimized Python. The independent code imports no
producer implementation. A false extension to525 is rejected, as is a changed
dual multiplier `lambda=d/63`, in both modes. These are controls of the proof
certificate, not counterexamples at525. The older748 proof and its tests are
not rewritten to claim526. Accepted Session138 inputs/equality are inherited,
not all re-proved in this release.

## Second checkpoint: all-modulus operational arithmetic

[Proof](all_moduli/PROOF.md) · [Core/orbit audit](independent/ALL_MODULI_AUDIT.md) · [Padding audit](independent/ALL_MODULI_PADDING_AUDIT.md)

**ACCEPTED.** The all-M equivalence uses integer Newton coefficients, not
division by factorials modulo a composite. Exact transitive periods and
transients use prime-power valuation arguments plus CRT. Universal source
padding has sharp index A(M)-1 and period L(M); the individual delta-dependent
onset is a sufficient bound only. The actual image of Psi is a finite family,
not an asserted ring or a realization of all arbitrary functions.

The exact class counts5,49,685,322101 at primes3,5,7,11 combine all-order
upper bounds and complete finite constructor-quotient reachability. The
producer uses parity-tagged half-vectors; the independent orbit audit uses
full p-entry functions. Its raw unary counts are5,48,685,322100, not the
producer's tagged counts6,48,686,322100. Their final contextual counts agree
after identifying zero tags and adjoining a disjoint-support join at5,11.
All-depth unary zero exclusion at5,11 follows from invertibility, not a
bounded search. These results do not extrapolate to every prime/composite.

Independent executable coverage:

- 209 moduli for exact first-repeat clocks via Newton power sums rather than
  the producer's Stirling recurrence; literal order-nine H55 profile.
- 7,776 signed degree-five profile/modulus pairs, plus full-function orbit
  exhaustion at3,5,7,11.
- All1,099 labeled nonempty tournaments through order5 via literal vertex
  permutations and cuts, rather than producer Hamilton-subset/partition DP;
  70,336 actual tournament/modulus comparisons, 89,402 padding comparisons
  over132 moduli, and the actual E4 response117/172 modulo275.
- Producer separately checks519 clock moduli,70,588 period comparisons,
  35,168 Psi comparisons,225 actual operator/profile checks and its full
  constructor quotient orbits.

All five inherited source blob pins match the declared source commit. The
base theorem's proof/verifier/receipts/pins are verbatim. Two producer summary
files only lose references to the excluded later CONTEXT_CONTENT addendum;
[acceptance and curated selector](all_moduli/ACCEPTANCE.md) document that cut.

## Fresh second-checkpoint relocation and failure controls

```bash
# From this packet directory:
sha256sum -c SHA256SUMS
python3 replay_checkpoint2.py --output /tmp/session139-checkpoint2-replay.json
python3 replay.py --output /tmp/session139-checkpoint1-replay.json
```

The second replay copies six standalone kernels to a new temporary directory.
It performs **24 executions:12 positive runs and12 expected rejections**.
Each normal/optimized output pair must agree and match the frozen receipt.
False threshold, dual multiplier, dyadic period, orbit cardinality and actual
E4-response mutations are rejected by explicit runtime checks, including
under `python3 -O`. All second-checkpoint kernels use the standard library.
Source/output hashes are recorded in [the receipt](REPLAY_CHECKPOINT2.json).

The original first replay remains unchanged:22 executions, eight positive
runs and fourteen expected rejections. Its historical polynomial cross-check
requires SymPy; the two stability kernels and symbolic bridge audit use the
standard library. The archived `independent/SHA256SUMS` has original absolute
scratch paths; use the cumulative root manifest for portable validation.

## First checkpoint retained: n>=748, spectral loss and history

The [n748 audit](independent_threshold/AUDIT.md) accepts the original shared
wedge-budget theorem and equality;20 producer and22 independent exact checks
remain recorded. False747 extension and producer h=0 are rejected. The proof's
original review-pending sentence is its drafting snapshot, not current status.

The [bridge audit](independent/AUDIT.md) independently reconstructs the full
bivariate rational conjugacy, determinant and Hamiltonian polynomial with
standard-library arithmetic. The all-one-root conjugacy proves the full
rooted response is preserved. M23 repairs exact H at n5; binary first-blindness
remains6. Auxiliary binary cocycle controls were retained without promoting
the separate later integral/weighted cocycle continuation.

The [historical checks](riemann_five_eleven/verify_history.py) reproduce the
corrected transfer/tetranacci/pentanacci discriminants, Bernoulli denominator
repeat, Paley11 character sums/ranking33 versus35, and a separate Held–Karp
rational witness. Twenty historical source pins and six inherited input hashes
are preserved in [the first pin record](SOURCE_PINS.json). The exploratory
probe is non-load-bearing; no global priority, RH/GRH or spectral-proof
impossibility is inferred.

## Documentation, scope and publication checks

Native `agents/check_docs.py` needs baseline documents absent from this sparse
checkout. Its unmodified same-version checker passes with read-only overlays
on the clean materialized Session137 ancestor, with every intervening path
accounted for. No baseline mutation or hydration is used. The final receipt
is [checkpoint-two documentation verification](DOC_OVERLAY_VERIFICATION_CHECKPOINT2.json);
startup remains179966/180000 bytes and frontier450/450 lines. New links,
source hashes, frozen prior lanes, conflict markers and explicit staging scope
are independently checked.

No central mistakes ledger or older canon file is changed here. Prepared local
commits are not remote publication: the publishing agent separately verifies
the actual main commit and every owned blob. Cocycle, regular-critical-line,
context-content and local-factor-hostile work require their own audited release.
