# Audit and reproduction

**Status:** exact algebraic proofs self-reviewed; finite certificates checked
by distinct implementations. No independent human/subagent proof audit, Lean
formalization, novel-priority claim, or universal Collatz proof is asserted.

## Universes and controls

- `content_renewal.py`:506587 explicit gates. All21844 valuation words with
  letters1–4 and lengths1–7 check two carry constructions, every split's gcd
  law and the carry envelope.5001 positive odd sources through10001 each
  retain24 actual odd steps, endpoints and clocks.5096 additional words
  (all positive compositions of total clock at most12, plus1001 actual
  48-step source words) check fixed-prime packing and four-layer linear debt.
  Four unrestricted minimum-cost certificates retain source weights. The
  independent lower-cost audit enumerates only necessary scalar clock pairs
  and their357+1 possible words.
- `five_potential.cpp`: warning-clean C++17, O0 and O2.325520 total affine
  states across the four moduli, with3879120 checked Bellman inequalities.
  The bound for larger letters is checked symbolically from integer maxima.
  The positive-word equality control proves the reported coefficient sharp.
- `audit_potential.py`: separately constructs Python tuple states, modular
  inverses, rewards and transitions, then checks all3879120 inequalities on
  the C++ vectors. The four vector hashes are in `potential_audit.json`.
  C++ O0 and O2 vectors and receipts match byte-for-byte.
- `flat_compiler.py`:5428 gates. Every indicator is compared pointwise on the
  full1024-element five-vertex orientation cube. Independent literal legal
  moves and full relabeling orbits audit the family image and class quotient.
  A direct arrow-based module oracle audits the rank-based prime filter.
  The normalizer/full-action discrepancy and the retained-port count are
  hostile and positive quotient controls, respectively.

Python normal and `-O` receipts match byte-for-byte for both main scripts.
No check uses a Python assertion that optimization could remove.

## Commands

Run from the repository root. C++ writes the temporary potential vectors;
Python audits those exact vectors and reports their content hashes.

```sh
python3 -B 05-knowledge/results/anchor_content_flat_quotients_20261030/content_renewal.py
python3 -O -B 05-knowledge/results/anchor_content_flat_quotients_20261030/content_renewal.py --output /tmp/content-renewal-optimized.json
python3 -B 05-knowledge/results/anchor_content_flat_quotients_20261030/flat_compiler.py
python3 -O -B 05-knowledge/results/anchor_content_flat_quotients_20261030/flat_compiler.py --output /tmp/flat-optimized.json
c++ -std=c++17 -O2 -Wall -Wextra -pedantic 05-knowledge/results/anchor_content_flat_quotients_20261030/five_potential.cpp -o /tmp/five-potential-check
/tmp/five-potential-check
python3 -B 05-knowledge/results/anchor_content_flat_quotients_20261030/audit_potential.py
```

For the compiler-optimization comparison compile again with O0 and run with
the optional output prefix `/tmp/five-potential-O0`. Compare all four text
vectors and the JSON stdout. Temporary vectors are reproducible from the
included source, rather than adding312500-line artifacts to the repository.

## Mathematical validity gate

1. `q|content` means the *full* affine identity, not one endpoint returning.
   The original n and actual endpoint are linked by equation(3).
2. The word cost is `sum a_i`; reducing a transition modulo a finite group
   does not reduce this integer cost. The unlimited-letter tail proves this.
3. A potential with identity at its maximum gives a nonpositive telescope;
   reversing that inequality would invalidate the bound. Both audits check
   the required direction explicitly.
4. The first four valuation layers are not all layers. The proposed uniform
   lift is named OPEN, and even it would leave other-prime and endpoint debt.
5. Union replaces descriptions by objects. Burnside is applied only after
   full-action saturation; marked ports restrict that action. Histories and
   noninvariant source weights cannot silently use this quotient.
6. Prime filtering uses every nontrivial proper module, not degree regularity.
   The coalescing operation preserves membership, not construction history.
7. Existing all-order reset/root-switch and Burnside methods are credited.
   No claim of a new scalar A000568 recurrence or a performance improvement
   follows from this five-vertex family test.

The failed generic clock formula at25 is recorded and repaired in
`01-canon/MISTAKES.md`. The documentation route uses the preceding topic
packet instead of enlarging the bounded startup surface.
