# Audit and reproduction boundary

2026-10-26. **PASS: independent mathematical reviews and finite exact controls.**
Proof review establishes the stated all-parameter deductions; finite runs
test the explicitly bounded universes. Hashes establish byte identity, not
mathematical validity. No Lean formalization or historical-priority claim
is made.

## Review separation

| Lane | Producer checks | Independent checks | Mathematical review |
|---|---:|---:|---|
|Tournament support and cuts|321237|205584|[P1–P4 review](audit_paths/PROOF_AUDIT.md)|
|Growing prime sieve|380043|490523|[W1–W3 review](audit_sieve/REPORT.md)|
|Paley cover and arithmetic collisions|18007|6499|[Field review](audit_field/AUDIT.md)|

The independent tournament program uses subset endpoint dynamic programming,
with literal permutation comparison through order5. It checks every labelled
tournament through order6, including the strong-only deletion lemma, exact
reconstruction, and the transitive and five-vertex hostile examples. Its
205584 checks are additional to the separate mathematical review; the root
manifest binds this program and its receipts, which the proof-review
submanifest deliberately does not claim to bind.

The sieve review independently checks direct primality, corrected self-hits,
wheel counts, recursive counts and first false-positive boundaries. The
field review independently constructs the pole tournaments, enumerates their
path counts and automorphisms, checks the modular cover by a different
residue calculation, and verifies the arithmetic collision factors.
Each of the six programs has byte-identical normal and optimized Python
receipts. Repeated execution is not counted as additional independent evidence.

The [integration review](audit_integration/REPORT.md) separately audits
README.md and METHOD.md, including their additional automorphism, language,
Paley-support and state-count consequences. It requested one wording repair:
augmentation counts monomials and prime evaluation computes a weighted sum;
neither scalar generally recovers the underlying shape or labels. The
repaired final files are bound by that review's manifest.

## Classical input and source versions

The tournament inverse uses the necessary separator condition in Arthur H.
Busch's 2005 thesis, Theorem4.12, printed pages45–47. The statement and proof
were read in the primary PDF, and the rendered statement inspected. This is
a CITED dependency, not a newly proved structural theorem or an iff test.
The [source record](paths/SOURCES.json) binds the PDF and identifies the
later Busch–Jacobson–Reid publication. The sieve is an elementary
specialization of classical inclusion-exclusion and the square-root cutoff;
the source and attribution boundary are in [SOURCES.md](sieve/SOURCES.md).

Inherited repository material is pinned to research base
`2c85cab14a5e3bf229d3c3993247ed17b2f92d97`. These historical pins remain
unchanged when a maintained theorem, index or navigation file advances.
The verifier reads their Git blobs at that revision, rather than silently
binding a later version. The packet's own files are bound at their final
reviewed bytes. Status-and-link-only promotions preserve the earlier reviewed
hashes in the lane manifests; no mathematical body was changed on promotion.

## Concurrent corrections and integration

The incoming S60 merge had restored unqualified all-k Mersenne-fibre
headlines while retaining a conditional proof. CITED provenance did not
discharge the missing distinct-prime Goldbach assumption. The
[regression record](REGRESSION_CORRECTIONS.json) specifies eighteen exact
replacements in six files, including duplicate removal and positive-domain
repairs. The [bounded correction audit](audit_sieve/REGRESSION_AUDIT.md)
replays them and preserves the older audit's other conclusions. It does not
claim to re-audit every assertion in THM-4805.

The separate integration record binds the START-HERE route and the promotion
of strong-atom C5 from OPEN to the linked all-order named-guard theorem.
Fresh, unnamed observations remain outside that theorem. Current maintained
files may subsequently advance; replaying the historical patch and checking
its result hash is distinct from requiring the live file to remain frozen.

## Reproduce

From the repository root:

```sh
python3 -I -B 05-knowledge/results/growing_sieve_path_20261026/verify_packet.py
```

This verifies the closed packet file universe, inherited Git source blobs,
review bindings, exact correction and integration replays, and reruns all
six programs normally and with `-O`, comparing output bytes with their saved
receipts. `--integrity-only` omits computational replay. `--live` additionally
requires the maintained files to equal this session's declared after-bytes;
use it for publication, not as a permanent ban on subsequent improvements.
The external primary PDF is not downloaded by the verifier; its reading
record and checksum are provenance rather than a network dependency.

The combined counts are 719287 producer and 702606 independent finite checks.
They support the stated universes and do not prove twin-prime infinitude,
abc, the unrestricted tournament spectrum, fresh quotient minimality or
black-core connectivity.
