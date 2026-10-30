# Proof and computation audit

Status: **PROVED scoped deductions**, **FINITE-EXACT controls**,
**CITED background**, **VERIFIED external replay**. Same mathematical
agent wrote and reviewed the proofs; the separate code paths do not
constitute an external review or a Lean proof.

## Quantifiers and mathematical gates

- C1 is for every odd carousel, subset, and both budget parities.
  The proof derives the at-most-two-adjacent forced marks, deals
  separately with at most three marks, and checks cyclic wraparound.
- C2's induction covers every nonnegative integer weight, including
  zero slots; its objective splits exactly, not in expectation.
- C3 depends on Liu's Theorem1.1 AND prime substitution closure. The
  second prime case uses superadditivity/C3 reduction, not a claimed
  chart for all U5-free primes. Equality is a literal recursive family.
- G1 needs a supplied proper b-fold coloring, a>=2b, an independence
  bound, and singleton as well as edge inequalities. An isolated vertex
  is the hostile for forgetting singleton constraints.
- D in G2/G3 requires all five children nonempty at a C5 root. Its
  complement closure is used; hereditary closure is not claimed. The
  exact first four rows include the exceptional B=1 and B=2 cases.
- G2 fixes slopes and a sufficient divisibility progression. It does
  not prove all-periodicity or the least period. G4 explicitly shows
  that integrality of the scalar bound is insufficient for equality.
- N3 is a necessary condition on a hypothetical object already excluded
  by the author's stronger theorem. It is not a new Ramsey improvement.
  The catalogue has two42-edge types, so replacing its maximum by41
  is explicitly rejected. Both full source and rooted occurrences remain.
- Frames used by the carousel chart are choices. No canonical
  automorphism-equivariant selector or unnamed metagraph compression
  is asserted. The81/27 control retains its distinct full groups.

## Reproduction from the repository root

```text
python -X utf8 -B 05-knowledge/results/ramsey_neighbourhood_lifts_20261030/check.py
python -X utf8 -O -B 05-knowledge/results/ramsey_neighbourhood_lifts_20261030/check.py
python -X utf8 -B 05-knowledge/results/ramsey_neighbourhood_lifts_20261030/audit.py
python -X utf8 -O -B 05-knowledge/results/ramsey_neighbourhood_lifts_20261030/audit.py
python -X utf8 -B 05-knowledge/results/ramsey_neighbourhood_lifts_20261030/seal.py
```

The producer has201681 explicit internal gates, plus a receipt-equality
gate on replay. The separate auditor has144912 gates. All survive `-O`.
The87376 complete charts themselves are checked by direct inequalities
inside `carousel.check`; those checks are additional to the gate counter.

The producer's exact universe is:

1. All2^(n+1) parity/subset charts at each odd n=3,5,...,15;
   42274 cases also check forced-mark geometry from actual intervals.
   Seed20261030 supplies250 integer-weight vectors at each of
   n=17,19,31,65,127, with coordinate range0..40. These1250 larger
   cases are sampled controls, not exhaustive assertions.
2. Every labelled tournament through6, filtered only by literal induced
   U5 membership (all relabelings). Transitive size uses a source
   recursion. Paley7 supplies actual U5 occurrences and fails the
   unrestricted bound. Sixteen Pascal-cell matrices attain G(k).
3. All recursive decomposition options with A=1..4,B=1..8, including
   join, union and every positive full-C5 budget tuple. DP values do
   not use the proposed four formulas. The311604 C5 candidates plus
   join/union choices give32 witnesses; separate native bitset clique
   searches test their actual relations and complements.
4. All feasible C5 integer vectors in{0,...,K}^5 for K=1..12,
   totaling130388, plus K2 superadditivity on1..64 squared.
5. Literal reconstruction of Session120's two N=231 spectral minimizers;
   naive autocorrelation reproduces their exact fourth moments, then
   induced five-vertex search produces U5 witnesses. No computation
   of their largest transitive subtournaments is claimed.
6. The inherited native81/27 runtime is executed, and the seven-type
   critical17 catalogue receipt supplies exact edge counts. The latter
   is an inherited exact bank, not a fresh exhaustive graph enumeration.

The separate auditor imports no producer module. It recovers maximal
carousel supports from absence of literal cyclic triangles through11;
tests G via its direct popcount sum; independently computes w by the
classical maximum recurrence through256; checks recursive slope
attainers by exact composition at A<=256; and checks residue and
nonuniform-root inequalities used in the proofs. Its large recursive
graphs are represented by their proved construction laws, not materialized.
Native materialized graph checks belong to the producer's32 witnesses.

## External certificate boundary

See N1 in [NEIGHBOURHOODS.md](NEIGHBOURHOODS.md). The author's exact
checkers were inspected and run at the pinned revision. We do not call
that an independent reimplementation. The replay wrapper explicitly
removes `PYTHONOPTIMIZE`, because the external scripts use assertions.
Its successful scope is17985 arithmetic checks and three exact flag
exclusions. Full source hashes and outputs are retained.

No external service was written to, no challenge submission made, and
no generated benchmark was passed off as a literature result.

## Packaging corrections during replay

The first receipt comparison exposed tuple-versus-JSON-list serialization,
not a failed mathematical gate; canonical JSON conversion now precedes
comparison, and both normal and optimized replays pass. Updating the old
saturated-port README also exposed its seal comparing a historical incoming
README pin with a subsequently amended live README. The seal still checks
the immutable historical Git blob, but no longer demands that a nonruntime
summary never change. Its existing live source checks remain. The old
packet's artifacts were resealed after this source-audit clarification.
