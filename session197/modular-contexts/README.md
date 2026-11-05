# Finite modular whole-module contexts

The reviewed proof gives an exact finite congruence quotient for every modulus and every finite alphabet of intact-module tournament substitutions. It includes universal triangle recurrences and exactly49 global behavior classes modulo5 and685 modulo7. Nonzero universal-content classes can reach any scalar residue with at most4/5 strong D/E wrappers and10/13 added vertices respectively. These are modular statements, not exact-integer or prime-coverage results.

Run `python3 verify.py` and `python3 -O verify.py`. Both write only stdout and reproduce normal.json. The unchanged core.py is the earlier triangle-content verifier, used only as a helper here. No external Python packages are needed.

The independent checker is self-contained: run `python3 audit/check.py` and its `-O` form. Both reproduce audit/independent-normal.json. The independent review is audit/REPORT.md. Its source/ paths refer to the proof/verifier/helper at this directory's root; its source-replay receipts are retained under audit/.

The four mathematical dependencies are copied unchanged under dependencies/. The frozen proof's path-polynomial PROOF.md is dependencies/CENTRAL_OPERATORS.md; its triangle-content/PROOF.md is dependencies/TRIANGLE_CONTENT.md. No exploratory prime-power probes or later addenda are included in this release.
