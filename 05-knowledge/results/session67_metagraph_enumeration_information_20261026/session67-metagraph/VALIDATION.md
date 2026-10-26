# Independent verifier receipt

**PASS.** `verify_core_families.py` is a self-contained Python standard-library
implementation, with no imports from inherited repository code, no external
executables, and no downloaded binary execution.

Commands run:

```
python3 verify_core_families.py --max-order 40 --class-order 7 --out receipt.json
python3 -O verify_core_families.py --max-order 40 --class-order 7 --out receipt-optimized.json
cmp receipt.json receipt-optimized.json
python3 -m py_compile verify_core_families.py
```

Both receipts are byte-identical, SHA-256
`b34005c58804886147298ccb2e152b1522f82ba8ee9eb5e7e8c5732f4a5d8c5e`.
All semantic checks raise explicit exceptions and survive optimization.

- Arithmetic: every odd permutation cycle type through order 40; complete
  core-size triangle, diagonal pair-free counts, and row totals.
- Independent arithmetic: direct expansion of the bivariate substituted
  Davis series through order 14, versus signed fixed-core inversion followed
  by positive weight assembly.
- Structural: all classes through order 7, produced by arc-flip BFS.
- Controls: canonicalization against full permutation minima through 6;
  labelled orbit mass at each order; actual automorphism cycle indices of
  all pair-free cores through 7; canonical chain extraction and inflation;
  complete positive-weight orbit coverage and no duplicate families;
  weighted-stabilizer automorphism law on every generated class.
- Inequalities: all 780 admissible `(n,d)` core-deficit tail bounds and all
  820 `(n,q)` uniform row bounds through 40, using exact rational arithmetic.
- Hostiles: transitive-triple pair relation is nontransitive; nonreduced TT2
  quotient weights merge; equal group order but different vertex-cycle action
  gives distinct weight counts; TT_n line partners have full pair-free core
  for every checked n=5..30.

Known external sequence entries through order 12 are secondary controls.
Completeness of the structural universe follows independently from labelled
arc-cube connectedness and equivariant transport of arc flips, not from
matching these entries. Arithmetic uses inherited Davis fixed-point counts;
there is no claim of an independent proof of that classical formula.

No actual class census above 7 was run here. The order-40 outputs are exact
formula evaluations, not individual tournament generation. The older n=8
structural results and 52-attachment growth hostile are recovered source
claims, clearly marked as inherited rather than rerun.
