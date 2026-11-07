# Session216 structural paired-source extension

New work, separate from the frozen carriers and complete order9 packets. Independent audit requested.

- LEFT_RESET.md: positive physical Q-reset paired ray; broader five-state binary criterion; nonnegative conserved defect.
- TWO_ADIC_AND_PARITY.md: stronger necessary universal-pair congruences and an exact opposite-parity global-moment counterexample.
- FINITE_SEARCHES.md: precisely bounded periodic-macro and counter-rotation exclusions.

No eventual spectrum theorem, F(8,137) certificate, universal paired seed, or full binary Q/Z paired seed is claimed.

Reproduction: Python3 standard library; C++17 compiler. Run from this directory. The frozen sibling packets provide the audited BANK.txt and BANK9.txt data. The provenance receipt records hashes of external inputs and local deliverables. Avoid interpreting raw search output files as proof beyond the documented finite scopes.

Core controls:

python verify_left_reset.py > LEFT_RESET_RESULT.recheck.json
python -O verify_left_reset.py > LEFT_RESET_RESULT.opt.recheck.json
python verify_structural.py > STRUCTURAL_CONTROL.recheck.json
g++ -O3 -std=c++17 verify_global_product.cpp -o /tmp/verify_global_product
/tmp/verify_global_product

The last program's input bank paths are relative sibling dependencies; its direct graph-count certificate is already retained. All discovery code is retained for reproducibility, but large regenerated intermediate source lists need not be published.
