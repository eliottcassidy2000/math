# Independent audit of current endpoint field codecs

## Verdict and scope

ACCEPT `phases/candidate_phases.py/json` and `phases/continuation_phases.py/json` as exact current-observer computations at the supplied physical tournaments. The formulas do not transport an ancestral spectral clock and do not establish tournament automorphisms or justify phase resets. Those exclusions are expressly stated in the files.

## Pinned source authentication

The audit fetched `05-knowledge/results/return6325_completion_20261103/transport.py` at commit `4ad839b20dd280cdcfbadf86cb62b2be3c9615ad`. Its blob SHA is `9ddbd5b5c90eb4dcf9121128ef34c63a76eb507d`. Its chart and cubic multiplication formulas agree exactly with the new codec implementation.

The same commit's `05-knowledge/results/triangle_return_170807_20261103/audit.json`, blob SHA `0f4e2a1881c07e2ac716b4c6f26f9e7a30b92e1a`, was fetched and compared structurally against `ports/input_audit.json`: the full objects, including all four historical states, agree. The local raw bytes have different whitespace; no byte-identical claim is made.

## Independent algebra

`audit_codecs.py` imports none of the author's codec code. It uses generic polynomial long reduction and linear repeated multiplication rather than the binary-power implementation. It checks:

- x^3+x^2+1 has no roots in F7, hence is irreducible;
- every one of the 343 coefficient vectors obeys z^343=z;
- both the inverse coordinate transformation and the F49/charge involution hold on all 343 vectors;
- the cubic Frobenius fixed set has 7 vectors while the F49-plus-charge fixed set has 49 vectors;
- all four historical saved charts agree exactly;
- all five 512423 charts and all six 1537271 charts agree with independent powering and with their named endpoint data;
- both candidates' recorded input SHA-256 values match their actual files.

The differing fixed-set sizes are an explicit control against identifying the two actions. At 512423, physical port50 is F49-fixed but has a three-element cubic orbit. At 1537271, physical root39 is genuinely fixed in both current charts, while the other five retained ports have cubic degree three. This is a computed endpoint fact, not an inferred reset or ancestral-clock relation.

## Physical-data binding

Each golden vector is reconstructed from the named endpoint fields as

    (H, (S2-E1)/2, (S2+E1)/2).

Integrality is checked. The candidate parser's actual `PORT` records agree with the freshly generated named fields in `ports/candidate_joint.json`. The continuation is bound directly to `ports/continuation_joint.json` and its physical-name candidate file. Both sets were replayed with the independent C++ endpoint DP through `ports/replay_joint.py`, including the historical fields and the 30 continuation directed joint responses.

The codec consumes the indicated marginal endpoint coordinates from the same physical graph. It does not reconstruct or replace the full joint tensor, which remains retained separately. The native guard, complement, and deletion proofs remain separate obligations.

## Replays

Both original codec scripts were run in normal and optimized Python. Their output JSON and stdout byte-match. The named joint replay also byte-matches normal and optimized output. The independent algebra checker normal and optimized stdout byte-match.

Artifacts: `audit_codecs.py` and canonical receipt `codec_audit.json`. The authoritative codec and joint input JSON files are in the sibling `phases/`, `ports/` and `returns/` directories. Duplicate replay stdout is omitted from publication.

Final provenance refresh: the candidate codec gained input SHA-256 metadata only. All pre-existing JSON fields, including all physical states and charts, remain identical to the independently audited version. The three new digests were checked against the actual inputs.
