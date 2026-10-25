# Session55 independent integer-search audit

Verdict: PASS after verifier repair. No mathematical overclaim detected in the audited bounded results. No general MSOS nonexistence or search-record claim is supported or intended.

## Found and repaired

Initial documented --verify failed: integer histogram keys in memory did not equal string keys after JSON parsing. Producer changed histogram keys to strings. Re-ran normal and optimized full-bound --verify successfully. Initial failure preserved in verify.log, repaired runs in verify-fixed*.log. Independently generated normal and optimized certificates coincide, and structurally equal the producer certificate (JSON ordering is immaterial).

## Independent checks

independent.py enumerates actual square values and ordered top-corner choices at every center root 1 through 2000, constructs remaining entries from row equations and opposite pairs, and tests actual square membership and nine-value distinctness. It checks 22,532 corner pairs; no solution. Every opposite-pair bank agrees with the Euclidean sieve. It also tests all eight literal dihedral transformations for 320 distinct signed parameter examples, finding exactly one x>y>0 representative.

Both full bound-100000 searches were rerun, including every Gaussian count and brute root comparisons through 2000. The reported 15,919 primitive parameters, 161,436 insertions, 347,561 eligible gap pairs, maximum 40 gaps, and zero solutions reproduce.

Atlas checker rerun normally and with -O gives identical outputs: exact unoriented Q15..Q25 counts 1,1,1,0,0,0,0,0,3,0,10; 13 gaps at root 1105 with no compatible diamond; Q23 endpoint hostile and Q25 forbidden edge hold.

## Proof review

Dihedral action on corner deviations is the full signed permutation group. Nine distinct entries imply nonzero distinct absolute corner deviations; exactly one chamber x>y>0 remains. Four-gap distinctness reduces precisely to x!=2y.

Primitive opposite-square triples reduce to odd roots and coprime opposite-parity U=(b+a)/2,V=(b-a)/2. Both choices of which Pythagorean leg is larger give the same positive gap; common scaling k is unrestricted for each pair. Thus the sieve does not discard nonprimitive full arrays. m<=floor(sqrt(B)), n<m and k<=floor(B/(m²+n²)) include all centers <=B; impossible equality m²=B merely creates no admissible n.

Zero entries cannot occur in a distinct square solution because the opposite square would be 2e². Count formula subtracts exactly four diagonal representations and divides the remainder by eight, including all common factors.

Weak no-wrap endpoint equality T=p implies p=3,c=1 and contradicts nine-value distinctness. Integer no-wrap p<=3B² forces e<=B, so the claimed prime exclusion through 30 billion is correctly inclusive. This is not a bound on existence above that range. The first weak canonical PRIME witness 59 says nothing about genuine integer squarehood or extension fields.

Primitive mod24 proof is valid: modulo4 all roots have center parity; modulo3 rows are all zero or all units and columns force the same choice. Primitivity rules out either global divisibility. Odd nonmultiples of3 square to1 mod24.

Literature status review is limited to the atlas's careful separation of inherited PROVED/audited items from CITED all-order constructions and unaudited preprint computations. This audit does not independently replay external manuscript code or promote those citations.

## Scope and hygiene

Only this audit directory written. Shared math checkout was dirty and older; no synchronization or remote actions attempted. Startup script ran read-only; processor status blocked on missing .machine-id, and no registration was made. Parent owns publication. No external source dependency needed for independently checked elementary statements.
