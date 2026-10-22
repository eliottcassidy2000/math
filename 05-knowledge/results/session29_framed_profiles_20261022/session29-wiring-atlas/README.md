# Tournament boundary atlas

Open `atlas.html` in a modern browser. It is a complete portable file: no server, external library, network request, installation, or data download is needed. Choose an address in the upper-triangular grid or the native selector. Cell numbers are Hamiltonian-path counts H, not incidence counts. The inspector gives each source's exact blue/black multiplicities, target automorphism order, and full merged tiling degree versus the captured family incidence. Expand the witness section for rooted tile IDs and individual edge IDs.

Eight ochre-bordered cells are leaves in the original restricted boundary. They are **not global bridges**: the subsequent extension supplies alternate paths. All 324 selected lines have audited global nonbridge certificates. This atlas does not claim that the full order-10 metagraph has no bridges.

The all-r formulas enumerate the framed source family's targets, not all tournaments. Relative proof links work when this folder is beside `session29-target-wiring`, `session29-bridge-extension`, `session29-reconstruction-audit`, and `session29-detour-audit`. Standalone interaction works even without those folders.

## Reproduction

- `python build.py`: reads audited sibling JSON, checks 737 consistency conditions, regenerates `atlas-data.json`, `atlas.html`, and `validation.json`. No Git operations.
- `node test_structure.js`: runs all 91 selection handlers and the native selector's change handler in a minimal DOM harness. Writes `interaction-validation.json`.
- `python test_browser.py`: optional real Chromium/Playwright interaction and screenshot test, provided an environment that permits browser startup. Requires Chromium at `/usr/bin/chromium` and Python Playwright.

## Testing disclosure

Data validation and executable JavaScript interaction checks pass. A genuine Chromium run was attempted, but the sandbox prohibited its singleton socket; an escalated attempt failed in environment startup on a read-only mount. Therefore no visual screenshot, actual keyboard behavior, or responsive rendering has been verified. The minimal DOM test is not a substitute for browser rendering. The retained browser test script can perform that remaining check in a working environment.

## Files

`atlas.html` is the user-facing artifact. `template.html` and `build.py` are its reproducible source. `atlas-data.json` is the normalized embedded data. Validation results record source SHA-256 hashes. `SHA256SUMS` covers the deliverables.
