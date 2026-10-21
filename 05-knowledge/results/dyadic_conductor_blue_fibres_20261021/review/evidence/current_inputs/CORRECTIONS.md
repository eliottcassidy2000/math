# Retained scope correction

Date: 2026-10-21. Status: corrected before publication.

The first draft of section4 began with every conductor divisible by8, while its following square-kernel and genus-completeness deductions used the dyadic group structure from Theorem1. That broadened quantifier was false. The corrected scope is f=2^a, a≥3 throughout those deductions.

Exact hostile: at f=40, the class group has24 elements, the norm-modulo-8 kernel has6, and the square subgroup has3. The class of α=1+8ω has norm57≡1 mod8 but negative character (57/5), so cannot be square. Its reduced form is (25,20,52).

A related generic wording issue was corrected in section2: the residue-to-class map requires gcd(N(α),f)=1; odd norm alone is sufficient only when f is dyadic.

The intended dyadic theorem, all form-extension tables, blue-predicate bridge, and conductor16 obstructions are unchanged. The checker now guards its dyadic quotient helper and separately exhausts the f40 hostile. Both Python modes pass28,754 exact conditions with identical output.
