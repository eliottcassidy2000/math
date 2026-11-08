# Verification scope and release record

**Mathematical audit PASS for all four result groups.** The universal claims
are proved in their packets; the following finite universes are controls.

## Producer and independent universes

- Strong graft criterion: every labelled root-zero source and gadget of orders
  1–4, all 75×75=5,625 pairs. Ten nine-vertex H=9 separator outputs are counted
  by subset DP. Independent Warshall checks cover 2,175 source/root/gadget
  cases; independent permutations reproduce all ten output counts.
- Graft determinant/rank/inverse: 1,099 root-zero gadget presentations through
  order five, with all nine statistics independently counted by permutation
  splits; 4,392 literal grafts and 2,080 inverse checks. The separate audit
  covers 5,405 all-root presentations, 1,140 literal grafts and 10,304 rational
  inverse basis-vector checks. Minimum order four is proved, not inferred only
  from the census.
- Physical seed/descent theorem: 5,405 labelled rooted profiles through order
  five, 13,120 finite words and 167,895 exact congruence tests. Separate code
  tests 855 literal O/Z/J transitions, 812 state/cap envelope cases, and the
  five-line hostile on 500 odd values below 1,000.
- Modular theorem: 273 producer adjacent likelihood-ratio checks; direct branch
  enumeration at moduli 3,5,7,9,25,27; strict MTP2 and monotonicity at 15,45,1575.
  Independent conditional-sum code checks 385 adjacent steps, 324 valuation
  pairs, 20 conditional strict covariances and a direct mod-15 CRT enumeration.
  Independent integer interval arithmetic certifies the displayed Euler-product
  enclosure. Floating-point estimates are descriptive only.

The initial modular producer used assert-only checks. The auditor caught this
verification weakness; all were replaced by explicit exceptions, and both
normal and optimized versions were rerun. No mathematical statement required
correction. The final combined replay uses only optimization-safe checkers.

## Reproduction

From this directory:

    python3 verify_all.py > replay.json
    cmp replay.json REPLAY.json
    sha256sum -c SHA256SUMS

The driver runs eight scripts, each in normal and optimized isolated Python
mode, and compares every stdout byte with a committed reference receipt.
All code uses the standard library and requires no network or external data.
`graft/inherited_verify.py` is the unchanged Session225 graft verifier, copied
with attribution to make the packet portable. Other source dependencies are
proof references, not missing computation inputs.

Audit input hash ledgers retain the exact producer paths/bytes originally
reviewed. Final proof status headers and an explicit prior endpoint-defect
attribution were edited after review; these editorial changes are recorded in
the corresponding audit. The top-level SHA256SUMS covers final publication
bytes and is the portable integrity check.

## Repository validation limits

The session used an isolated sparse clone, leaving the older shared checkout
untouched. Explicit owned paths only were staged. The repository-wide
`agents/check_docs.py` was invoked but could not complete against this sparse
checkout: omitted theorem, hypothesis and reflection files triggered missing
file/status checks. No repository-wide documentation or CI pass is claimed.
The release instead checks the owned files for whitespace errors, Python
syntax, relative-link targets and byte-identical mathematical replay.

No Lean proof, global repository build, or enlarged scalar coverage census was
attempted. Original Session225 artifacts were not restored.

## Publication checkpoints

1. `54b1b8a9c93e60b0c1b47977978342ad1659fedf`: strongness theorem and finite
   first-birth separator; all seven files fetched back and byte-verified.
2. `368fb96114cb461c4c533a272823a1d6a483610d`: arbitrary physical seeds and
   bounded-intercept descent plus separate structural audit; all sixteen
   changed files fetched back and byte-verified.
3. The final checkpoint contains the audited determinant/inverse and modular
   results, combined replay, checksums, and final navigation. Its immutable
   commit is recorded in the delivery message and repository history; the
   file cannot contain the hash of the commit containing itself.

Every publication uses the current main as parent and a nonforce expected-head
update. Unrelated changes and prior recovery packets are preserved.
