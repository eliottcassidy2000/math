# Independent audit: native-annulus arrival budget

Audit date: 2026-10-21 UTC.

**VERDICT: PASS for the frozen packet's stated algebraic identities, its independently replayed finite source universes, and its explicitly limited conclusions. No universal annulus payment theorem is proved.**

Audited producer manifest SHA-256:
`8423a629603147141ef1a5f5570c7e74aecc6f55505f201c0d598d9aac6ad00a`.

The producer packet was not modified. This audit made no repository or remote changes. Publication is a separate action by the parent's sole authorized publisher.

## Mathematical target and dependency audit

The target is the Eisenstein norm of original sources already first-arrived in the inherited 627-point frontier, at each finite native P-edge horizon. It is neither a single real-row energy nor the quadratic norm of coalesced moving endpoints.

The original consumers were inspected directly:

- `05-knowledge/results/native_cancellation_105_20261018g.md`, section 3, equation (4), explicitly defines W and H from old/new **arrived** sources and states this energy increment. Its equations (2) and (3) give fixed named-source transport and exact new-source norm 153N. Its final paragraph leaves the increment bound and new-source entry unproved.
- `05-knowledge/results/root_eisenstein_first_arrival_20261018j.md`, section 3, defines the original-source, terminal-labelled arrival polynomial at finite depth, with full integer clock; its trace is six times the Eisenstein norm sum. The original note uses E for a different real-row energy, so the packet correctly specifies the norm consumer rather than relying on a shared letter. Section 5 already contains the three-endpoint source sweep and separate-scale all-prefix contraction census. Those are inherited mechanisms, not new results of this packet. Separate inequalities F(N,d)<=0 do not order F(64N,d) and F(N,d).
- `05-knowledge/results/root_sparse_boundary_module_20261017i.md`, R5, equation (18), has both frozen and moving endpoint norms and additional merge terms. It is a different consumer. The annulus packet retains the original, unaggregated norm of moving sources as collateral and does not claim to pay the coalesced moving-endpoint Gram.
- `eisenstein_circuit_completion_20261018g.md` and `guarded_children_105_20261018g.md` were authenticated as inherited dependencies. Packet completion and grounding 63 named children do not supply an unbounded energy or universal-entry theorem.

All 15 manifest members and all 10 listed original source files authenticated byte-for-byte. The copied frontier and rank data also equal the authenticated original receipt, not merely a claimed hash. All 56,219 inherited source rows at N=4096, 5997, 11002, and 16384 were independently rebuilt and compared literally, including source ID, compatible mod15/mod35 pair, terminal, full integer clock, P-depth, two charge coordinates, and norm.

## Algebraic checks

1. **Exact horizon ledger A/B:** expand Q(H) into distinct source diagonals and unordered same-terminal new/new pairs; beta(W,H) supplies each old/new pair once. Subtracting the exact source bank gives remaining original-norm collateral minus cross debt. Therefore the full 153N bound is equivalent to cross debt <= moving collateral. Cross debt <=0 is stronger, not necessary.
2. **Event ownership C:** every new diagonal is owned by its actual arrival depth; every admissible pair is owned by the maximum of its two depths. Equal-depth arrivals are one batch, not an arbitrary order. Explicit pair enumeration at N=7 agrees at all 257 cuts. The simultaneous-batch formula was also checked in 6,561 exact lattice cases.
3. **Mixed/partial packet formula D/E:** aggregated packet boundary vectors retain all between-packet cross terms. Complete zero sum in the enlarged ledger does not erase an old baseline subvector. The correction -Q(O)-beta(R,O) is necessary; the stated algebraic hostile gives -3, not 0. Formula D was checked in 6,561 exact lattice cases. Disjoint source ownership is essential.
4. **Endpoint sweep:** the only changing original odd points in M -> M+1 are 2M+1, 4M+1, 4M+3. Their primitive ancestors must be grouped before a source-norm change is computed. This is the previously established sweep, now applied to the source-specific cross-scale debt.
5. **Local coherence G/H:** both grouped-terminal and pairwise source formulas match exact F differences. The 64-step telescoping formula is valid but is not a closed scalar recurrence or a sign estimate. It still consumes native source charges, actual terminal vectors, arrival cuts, and integer clocks.
6. **Corrected scale notation:** a new source is charged with A_(64N), not zero-extended A_N. Equality of these charges is invoked only for old sources. The frozen text and calculations retain this distinction.

## Exact independently replayed universe

`independent_audit.py` imports no producer code. Unlike both producer charge routines and its incremental sweep, it scans every original odd integer m<4M and assigns the appropriate lower/upper shell weight to its primitive ancestor rho(m). It constructs its own array-backed first-hit route DAG and reconstructs all original source rows, including zero-charge sources.

The 16 tested **base scales** are:

`1,2,3,4,5,6,7,8,448,512,4096,5997,11002,16384,28672,32768`.

For each base it independently rebuilds the base and its 64-fold enlargement, yielding 28 distinct **complete source universes**:

`1,2,3,4,5,6,7,8,64,128,192,256,320,384,448,512,4096,5997,11002,16384,28672,32768,262144,383808,704128,1048576,1835008,2097152`.

This rebuild covers 9,648,115 original source rows across these universes and every horizon 0..256 at every tested base. The four detailed production ledgers and all six hostile-ray horizon ledgers match exactly. Every source arrives within the cap, so subsequent horizons are unchanged within these finite universes.

The independently reconstructed largest native route DAG has 5,081,980 nodes, maximum requested first-arrival depth 233, and SHA-256:

`f05865f56cca385744e01a861bd5fb01a8b0d5ab23b75eaaca858a0ca65daaeb`.

This equals the largest producer DAG. Every nonfrontier node is independently checked against its canonical edge, full integer clock addition, and strictly smaller suffix depth. Frontier nodes stop immediately at depth/clock zero. Cycles, missing suffixes, and cap violations fail explicitly; no unresolved route is dropped.

## Replayed consequence objects

- N=1..6 pass the full bank at all finite horizons. The N=7 first failure is d=8: old energy 19, new-only energy 1086, mixed term -2, arrived-new norm 1054, moving-new norm 17, total cross debt 30, increment 1084, budget 1071, bank -13. This refutes the unrestricted all-positive-integer statement, not the intended N>=4096 claim.
- At N=4096,d=1: new-only energy 2064, arrived-new norm 1860, new-only cross +204, mixed cross -1129, total cross -925, increment 935. Exactly 204 terminal fibres have positive local cross debt. Terminal 9's old sources are 9 and 6371; its new sources are 203889 and 407779. These source IDs and vectors are directly reconstructed.
- The all-new partial zero pair is literally sources 276145 and 92051, terminal 1, both clock 24, depths 4 and 5, charges (-2,-1) and (2,1). Its exposed depth-4 boundary has norm 3 despite completed norm zero.
- N=448,d=1 has cross debt +47 with moving collateral 66707. N=512,d=1 has cross debt +60 with moving collateral 76477. Thus the stronger arrived-only statement can fail while the full bank remains valid.
- All six printed ray peaks independently reproduce: `(N,peak,depth,budget)` = `(7,1097,22,1071)`, `(448,54514,92,68544)`, `(28672,2566744,156,4386816)`, `(8,1301,35,1224)`, `(512,62145,92,78336)`, `(32768,2858748,90,5013504)`.
- `audit_recurrence.py` separately checks 65 elementary source edits at 257 horizons using both grouped and pairwise formulas, and independently rebuilds the complete final source universes 4097 and 262208. Its G values and derivatives match the saved producer recurrence receipt. Complete-arrival G increases from -201674 to -201616, an increase of 58. A scalar monotone-in-N inference is therefore unavailable.
- `audit_printed_witnesses.py` parses the literal README witnesses and all six JSON log lines, verifies their receipt hashes/check counts, and checks the consequence fields in all 12,289 saved contiguous-sweep records. It recomputes the saved-record worst ratio by integer cross multiplication, obtaining N=5332, peak 581462 at depth 150, budget 815796.

## What was not independently rerun

The entire producer contiguous sweep over every base N=4096..16384 was **not** rerun by this audit. Its source, declared universe, deterministic receipt, full record coverage, scoped consequences, hashes, and normal/optimized logs were inspected and authenticated. Four selected complete base/enlarged source universes were independently rebuilt, not just sampled paths. The full six-ray universes were independently replayed. The global maximum over the entire contiguous interval remains the producer's reproducible FINITE-EXACT result, not a new independent exhaustive census.

The producer's recurrence program was additionally rerun unchanged in both normal and optimized Python, reproducing its exact receipt. Independent replay does not extend the finite universe or replace the need for a uniform proof.

## Reproduction and robustness

All checks are explicit exceptions, not Python assert statements. Independent normal and `-O` runs produced byte-identical receipts and logs:

- Main independent audit: 29,471,605 checks; receipt SHA-256 `0e180e1c0efd154abeb98c968ff37791b33258df68252a221bfc341467820859`.
- Independent recurrence: 2,272,117 checks; receipt SHA-256 `668fc8626275e9e8425c4b1f19b2416f874ed7317e86cf8231432659453b175a`.
- Producer recurrence replay: 6,929,714 checks; receipt SHA-256 `db9e5b12dfc8c59dd7c234e03865437c91ce992c53602398e6a1f26ff13bf7e7`.
- Printed-witness parser: byte-identical PASS logs normal and optimized.

Example commands, from any directory:

```sh
python3 -B /path/to/audit/independent_audit.py --package /path/to/producer-package --original-root /path/to/math-checkout --out /tmp/audit-normal.json
python3 -B -O /path/to/audit/independent_audit.py --package /path/to/producer-package --original-root /path/to/math-checkout --out /tmp/audit-optimized.json
cmp /tmp/audit-normal.json /tmp/audit-optimized.json
python3 -B /path/to/audit/audit_recurrence.py --package /path/to/producer-package --out /tmp/recurrence-normal.json
python3 -B -O /path/to/audit/audit_recurrence.py --package /path/to/producer-package --out /tmp/recurrence-optimized.json
cmp /tmp/recurrence-normal.json /tmp/recurrence-optimized.json
python3 -B /path/to/audit/audit_printed_witnesses.py --package /path/to/producer-package --audit /tmp/audit-normal.json
python3 -B -O /path/to/audit/audit_printed_witnesses.py --package /path/to/producer-package --audit /tmp/audit-normal.json
```

No third-party packages, network, or producer imports are needed. The optional `--original-root` is necessary to repeat authentication against the original ten source files and literal original-row comparisons. Without it, the complete independent source rebuild and packet-hash comparisons still run, but the report correctly omits original-file authentication and literal historical-row comparisons.

## Remaining frontier

**OPEN:** uniform full-bank payment for every N>=4096 and every finite arrival horizon; the stronger uniform source-specific coherence nonincrease; control of the separate coalesced moving-endpoint Gram; unbounded child grounding; universal Collatz entry.

No finite interval, finite set of 64-scale rays, exact source-bank identity, packet completion, or local recurrence closes these obligations. The frozen packet states this boundary correctly. No blocking mathematical or provenance discrepancy was found within its declared scope.
