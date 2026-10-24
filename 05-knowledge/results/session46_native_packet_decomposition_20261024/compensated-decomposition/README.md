# Compensated native packet decomposition

Main result: authentic added charged source options reduce the original sixteen-source circuit's optimal maximum move support from 16 to 6, while restoring every compensator. One additional source 3015 suffices. Two cost-two buffers 1394 and 2905 achieve the sharp twelve-clock budget. Both new eighteen-source circuits also admit four six-slot moves using one cost-two buffer.

Important boundary: these are unit-capacity routes permitting later reoccupation of a released resource. They do not satisfy strict temporal nonreturn. Under strict temporal nonreturn, restored buffers cannot change, and the primitive sixteen-/eighteen-slot atoms remain unavoidable. The distinction and exact reentry ledgers are part of the result, not a side condition omitted from the claim.

Read `THEOREM.md` first. `packet.json` contains the original-source matrices, literal paths, complete complement/converse lifts, shared/source-private owner symbols, and three main routes. `packet_eighteen.json` and `packet_eighteen37.json` add the cost-43 and cost-37 eighteen-source endpoint pairs.

## Verification

Run from this directory against the original repository data:

    python check.py --repo /path/to/math --out normal.json
    python -O check.py --repo /path/to/math --out optimized.json
    python check.py --repo /path/to/math --include-eighteen --include-eighteen37 --out all-normal.json
    python -O check.py --repo /path/to/math --include-eighteen --include-eighteen37 --out all-optimized.json
    cmp normal.json optimized.json
    cmp all-normal.json all-optimized.json
    python variants/check.py --repo /path/to/math --out variants/normal.json
    python -O variants/check.py --repo /path/to/math --out variants/optimized.json
    cmp variants/normal.json variants/optimized.json

The checker does not import or trust the search scripts. It reconstructs native matrices and full fixed-clock domains, verifies actual ordinary/converse classes and black multiplicities, checks all normalized resources, and multiplies every complete image witness over the original graph. It also enumerates exact guarded fibers and their move graphs, all twenty-seven four-slot trades, all eight single-cost-two-buffer menus, and all twenty-four single-buffer weighted menus. Explicit `need` checks remain active under `python -O`.

The twenty-nine full image vectors are in four files: twelve main vectors, four cost-43 extension vectors, four cost-37 extension vectors, and nine improved-route vectors. `variants/THEOREM.md` proves the exact event/work/peak-clock trade-offs; `variants/check.py` includes all eight certified route witnesses. `cost37/THEOREM.md` and `cost37/check.py` provide a focused entry point for that extension. To regenerate an image file, use `solve_images.py --repo /path/to/math --charges <charges file> --out <image file>`. The solver is untrusted; the independent full matrix multiplication is decisive. Its back substitution uses exact bit-plane popcounts, with no floating point or external solver dependency.

`inventory.json` and `circuit16.json` pin the inherited finite supply and original endpoint relation. The optional search scripts use these local inputs. `search_four.py` performs a redundant full-buffer reachability exclusion; the short integer invariant in the theorem makes this enumeration unnecessary for the principal lower bound. `SOURCE_PINS.json` records original inputs and the inherited exact dyadic theorem/audit. `MANIFEST.json` records the final artifact hashes and receipt equality.

## Scope

The positive routes require genuinely available resource alternatives and an unchanged, compatible, disjoint remainder balancing the enlarged initial packet's current. No such global remainder is constructed. The full canonical atlas is not rewired; its alternative resources are already occupied. No common-host amalgam, all ordinary per-port phase preservation, odd-primary preservation, or all-size native-supply/connectivity theorem is asserted.
