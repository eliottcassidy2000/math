# Independent review: finite modular whole-module states

PASS, no correction. This report is bound to the exact frozen MODULAR_CONTEXTS proof and verifier below. No source file was edited. This is a new audit of the finite congruence theorem, not an extension of the earlier triangle-content audit.

## All-order review

The free action on distinct internal exterior paths contributes precisely (r-1)!, while ordering the input path cover contributes r!. The actual matching contexts give the lower unit-triangular scalar response bank, so the equivalence is necessary as well as sufficient. The factorial cutoff is uniform in input order. Integer inversion followed by composition recovers every output weighted Newton coordinate. Each composition still has the entire original input as one module; finite row support and integrality survive parity changes. This establishes a genuine continuation congruence, not merely equality of the present scalar readout.

The central-coordinate change is integer unit triangular in the relevant direction; all later weighted Newton coordinates are integer combinations of the weighted central coordinates. Thus the weighted central transition columns have integral quotients before modular reduction. The proof never divides a nonunit in a residue ring. Prime-square cutoffs and the stated general prime-power upper bound have the correct odd/even boundary at k=p-1.

The Jacobi matrix orientation is correct: weighting turns the upper diagonal into1 and the lower into beta. The characteristic and sixteenth-power arguments give uniform modular recurrences for arbitrary actual inputs, without a Stieltjes or Hankel-sign assumption. The norm calculation in F49 gives alpha^8=-1 for the even quadratic factor and hence alpha^16=1.

The finite searches distinguish tagged states from global behavior classes. The two parity function spaces on Fp intersect only at zero. Below p all relevant factorials are units, and terms at r>=p vanish as functions on Fp; this justifies the parity-free function comparison. The C5 automorphism acts freely on Hamilton words in every whole-module context and supplies the missing mod5 zero class. Strong representatives follow from nonempty cyclic-wrapper words; the empty initial-state presentation is replaced by a nonempty return. Positive word length and added-vertex bounds are checked jointly, not optimized separately and then conflated.

## Independent computations

check.py imports no submitted code. It enumerates directed odd cycles as literal vertex words rooted at the least vertex, then enumerates disjoint cycle collections with the appropriate factor2 per cycle. These actual polynomials are independent of the submitted path-cover implementation.

- All29 rooted contexts through order3 and75 actual input tournaments through order4 give2175 physical context-profile controls.
- A total of203 weighted Newton columns are reconstructed from actual transitive substitutions. Their integer divisibility is checked before reducing modulo any modulus.
- These controls give19575 exact transition checks modulo2,3,4,5,7,9,16,25 and49.
- All four D/E transition matrices at each of mod5 and mod7 are reconstructed from actual graph profiles, not copied from the source tables. The integer three-state Jacobi matrices are recovered the same way.
- A separate complete finite graph search reproduces48 and686 tagged states, maximum shortest lengths7 and14, and the full nonzero communicating components of sizes48 and684.
- All240 mod5 and4788 mod7 nonzero-state/scalar targets have a positive word satisfying both claimed bounds. Independent polynomial-function inventories give49 and685 global behavior classes.
- Actual cycle collections verify the C5 zero profile. Exact matrix arithmetic checks both characteristic cubics, the odd mod7 recurrence and the even sixteenth power.
- In total781 distinct actual graph profiles are cached. Independent normal and -O receipts agree byte-for-byte. The copied source normal and -O runs reproduce the frozen source receipt exactly.

The result is modular and whole-module-specific. It does not prove exact integer coverage, prime production, or closure under arbitrary port-dependent births/deletion. Zero absorption refers to universal content, not merely scalar divisibility. Formal residue vectors are used only as an upper bound until the finite physical orbit establishes realization.

## Reproduction and exact hashes

Run python3 check.py and python3 -O check.py. The independent checker has no file inputs. Source replay uses only source/verify.py and source/core.py. The named mathematical dependencies remain the separately reviewed universal-content, central-operator and triangle-response packets identified in the submitted proof.

- source/PROOF.md: 1e7122c4573fa725e38dec08dbc811d7ac3f0b3dc91edd54c9851337ddbc16a7
- source/verify.py: 35c97f8a6a860f756b93f4a303390437919bcb061465042a07eca6c72c7034ef
- source/core.py: bbd76a9897f30d10bd846e39468ab644098d19d7ed9f36c90241f6bcb357fbcc
- source-normal.json: 13ac33a469f8cec306f8203d1b057f513c7a48dd529fe4d0cbbe879b97574ef7
- source-optimized.json: 13ac33a469f8cec306f8203d1b057f513c7a48dd529fe4d0cbbe879b97574ef7
- check.py: 64a09c3aacae23301bbc6bc8001f246e08095ec8e34f2c3445b34a4d5413788f
- independent-normal.json: 3d63b234dcb830c2b7b5624efebadad30d4f1d16250203aa127ad4bf3dd5d4d9
- independent-optimized.json: 3d63b234dcb830c2b7b5624efebadad30d4f1d16250203aa127ad4bf3dd5d4d9
