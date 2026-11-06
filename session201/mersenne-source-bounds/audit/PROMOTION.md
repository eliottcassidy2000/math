# Publication promotion audit

Verdict: PASS. The release packet preserves the independently audited theorem
and is suitable for publication in its stated scope. GitHub publication and
remote readback remain the parent's responsibility and are not asserted here.

## Exact changes

1. The proof differs from audit/MERSENNE-before-promotion.md only in the status
   paragraph: RESERVED / UNREVIEWED becomes PROVED, independently audited in
   the stated bounded-source scope. The OPEN occupancy disclaimer remains.
   Regenerating the unified diff matches PROMOTION.diff byte-for-byte.
2. audit/check.py differs from the original independently authored checker
   in exactly one string: the MENUS.json input path now points to bundled
   dependencies/MENUS.json. No algorithm, test domain, or check was changed.
3. Author verifier, TWELVE.json, PROVENANCE.json, both author output files,
   original audit report, original audited hashes, independent RESULT.json,
   and original author-replay.json are byte-identical to the audited inputs.
4. Bundled MENUS.json and CORE_PROOF.md match the reviewed compact-central
   release byte-for-byte. RELEASE-SOURCES.json hashes verify. TWELVE.json is
   exactly MENUS.json key 5. The original audited proof hash resolves to
   audit/MERSENNE-before-promotion.md, not the status-promoted proof.
5. README.md accurately states the bounded-both-g-coordinate theorem,
   non-consequences, reproducibility commands, and computation sizes. There
   is no unsupported claim of occupancy or all-target nonrepresentation.

## Replays from the relocated release

With PYTHONDONTWRITEBYTECODE=1, all four commands passed:

- python3 verify_mersenne.py
- python3 -O verify_mersenne.py
- python3 audit/check.py
- python3 -O audit/check.py

Both author outputs match MERSENNE-normal.json byte-for-byte (189165 checks).
Both independent outputs match audit/RESULT.json byte-for-byte (13536024
checks). Fresh outputs are saved as PROMOTION-author-normal.json,
PROMOTION-author-optimized.json, PROMOTION-independent-normal.json, and
PROMOTION-independent-optimized.json in this audit directory.
The explicit RuntimeError-based checks remain active under Python -O.

## Scope

The result applies only to actual indices u,v>=2 in the stated twelve-template
family, n=2^k-1, and a simultaneous bound g(u),g(v)<=G with G>=1. It yields
n<=18630 G^8 and thus a necessary growing source g-coordinate along Mersenne
targets. It neither excludes every Mersenne target nor proves missing targets
when G grows, and it does not prove or refute uniform twelve-template occupancy.
The physical product-bound dependency is proved independently of occupancy.
Finite controls are not promoted into universal coverage evidence.

The release contains no unaudited exploratory inverse or selected-target scan.
This promotion audit changed no proof or executable code. Its additions are
this report, four exact replay outputs, and PROMOTION-HASHES.json. The latter
hashes every packet file other than itself, relative to the release directory.
