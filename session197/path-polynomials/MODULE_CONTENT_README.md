# Exact module contents and a physical9-to27-to81 hierarchy

Two independently reviewed theorems compute exact gcds of physical Hamilton counts. Universal intact-module content is gcd (r-1)!r!c_r and is certified by an explicit bank of only ceil(n/2) actual matching contexts. Restricted cyclic-context content is the gcd over all physical T,U of H(C3[T,U,V]).

In central-factorial coordinates, ordinary bin content, universal content and cyclic content are respectively gcd gamma_k u_k, gcd k!gamma_k u_k and gcd (k!)^2 gamma_k u_k. Their exact chain is g|G|C|H, with the sharp mixed relation gC|G^2 and bound G^2|g^3. A physical eight-vertex graph has g9,G27. A ten-vertex graph has g9,G27,C81; a directly counted16-vertex matching context escapes the extra cyclic divisor.

These results distinguish a true universal divisor from one inferred only from a restricted testing language. Gcds need not be attained counts. The proofs explicitly show that these valuation refinements add no prime/7p filter beyond ordinary content and induced-count monotonicity. They do not establish the full spectrum or an inverse for every multiple.

Both proofs, source verifiers and normal/-O receipts are included. The independent reports and standalone checkers are in `universal-content-review/` and `cyclic-content-review/`. The latter report binds the former's completed PASS, resolving the provisional review status recorded in the frozen cyclic proof. Run `python3 verify_universal_module_content.py` and `python3 verify_cyclic_module_content.py` from this directory, also with `-O`; independent reproduction commands are in the corresponding reports.

The source imports are listed in `MODULE_CONTENT_DEPENDENCIES.json`, including the already published Gram theorem and its audit. `MODULE_CONTENT_MANIFEST.json` binds every new file. Frozen proofs, verifiers and independent reviews are preserved byte for byte.
