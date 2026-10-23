# Independent review of the multivariable filtration supplement

PASS for the regular representation over Z/(3^k), with r,k>=1.

The root independently checked the exact integral identities D^(2s+1)=(-3)^s T^s D and D^(2s+2)=(-3)^s T^s D^2. Their tensor contents add exactly because each factor has an entry attaining its minimum3-adic valuation. A surviving monomial has total degree at most2(k-1)+2r; the exhibited exponents(2k,2,...,2) attain this bound. This proves both directions for the ideal, not only separate coordinate powers.

At maximal surviving degree every coordinate must be active with positive even exponent. The common factor3^(k-1) allows the remaining product to reduce modulo3, where D_i^2 is the cyclic norm. All survivors therefore yield the same signed full norm. Translations fix this norm, so its generated ideal has precisely three elements, not one dimension for each boundary monomial.

The proof correctly warns that sums of monomials may cancel, that faithfulness is needed for sharpness, and that a trivial representation has smaller nilpotency index. It does not use the nonunit denominator3, replace the original mod8 native guard, or infer legal-move connectivity from table equality.

The separate multivariable checker constructs full regular delta columns both by tensoring integer one-axis columns and by literal coordinate shifts, covering6882 boundary monomials for1<=r,k<=4. All20864 checks pass normally and with optimized Python. These finite tests are supplementary to the proof above.

Theorem SHA256: 25d661f4c42eb7eb71f52de1c12360f493145927ddc1fe7310bf0e4f84baacfe
