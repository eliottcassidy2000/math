# Independent audit of the #P gap certificates

**Verdict: ACCEPT** the counting-complexity consequence in
`global_referee/COUNTING_GAPS.md`, SHA-256
`b665291d2fd022d8019784fea2a279249eb37755ba6a339e9df3fce71be9ed87`.
The implementation reviewed and tested is `global_referee/counting_gaps.py`,
SHA-256 `4fa3b805dd247bd7d6219f93e3bf4ee4c1b20209b123bf6033a7e832ae4c1ca7`.

This is a distinct structural certificate, not an appeal to closure of #P
under subtraction. The four functions are delta=H-u, e=H-u-C_all,
g=H-u-|v|, and Gamma=H-u-|v|-C_proper on nonempty binary tournaments.
Their specified value at malformed/empty input is zero; in particular the
empty g/Gamma value is an extension, not their negative subtraction formula.

## Mathematical checks

The D/E/G/Q recurrences follow exactly from the independently audited
vertex-cycle identities. The e recurrence subtracts one per new cycle;
Gamma subtracts one only when its cycle remainder is nonempty. Hence its
leaf weight is `u_R-epsilon*v_R-1_[R nonempty]`, with empty weights 0 or 2,
not a negative coefficient. Nonempty leaf weights are nonnegative by the
already established integer algebraic positivity.

The sign epsilon is computed once from the original input. A top-level g
skip enters delta mode; a Gamma skip enters e mode. Every colored recursion
enters delta. Consequently no later signed recursion silently reselects the
sign. The two recursion colors represent a genuine factor two.

Every recursion removes vertices. A trace therefore has at most n records,
and the fixed-width n-record encoding is injective: all inactive fields,
unused path slots, and post-leaf records are forced to zero. The root is
canonical, and a directed cycle through it has a unique old-path ordering,
so no cyclic rotation factor or alternative padding encodes the same trace.
The initial mode is fixed by the function, not guessed nondeterministically.

Each determinant has absolute value at most n! by Leibniz; every leaf bound
is at most 2n!. The chosen fixed label width therefore suffices. Exact
determinants/signs and all syntactic checks are polynomial-time, and the
certificate length is O(n^2 log n). The separate diagnostic enumerator is
not claimed to be polynomial-time. These facts give actual #P verifiers.

## Independent implementation-under-test check

`audit_counting_certificates.py` intentionally loads the proposed verifier
and encoder as the target under test. It does **not** use the target's
expected-count function or determinant. Expected values come from separate
literal Hamiltonian permutations, minimum-rooted cycles, and Fraction
Gaussian elimination in this independent packet.

Complete controls: every labelled input of order 0..4, plus C3=>C3 and
source=>C3=>C3, totaling 78 inputs and 184 valid certificates. Every
single-bit mutation of every valid certificate was tested: 23,834 mutations.
Exactly 38 were other genuine valid certificates; every other mutation was
rejected. No duplicate encoding occurred. Wrong length, nonbinary symbols,
malformed matrices, and empty input were rejected as specified.

The two positive-Gamma controls independently give respectively
`(delta,e,g,Gamma)=(5,3,4,2)` and `(5,3,5,3)`, exercising recursive colors
and Gamma's skip-to-e branch. Normal and optimized audit outputs match.

Reproduce with the sibling referee packet present:

```
python3 -B audit_counting_certificates.py ../global_referee/counting_gaps.py > /tmp/counting-audit.json
python3 -B -O audit_counting_certificates.py ../global_referee/counting_gaps.py > /tmp/counting-audit-O.json
cmp /tmp/counting-audit.json counting_certificates.normal.json
cmp /tmp/counting-audit-O.json counting_certificates.normal.json
```

Required additional files beyond the global determinant audit subset:
`COUNTING_GAPS_AUDIT.md`, `audit_counting_certificates.py`,
`counting_certificates.normal.json`, `counting_certificates.optimized.json`,
and the already separately published target referee script. The determinant
and tournament-kernel dependencies are already in `GLOBAL_FILES.txt`.
