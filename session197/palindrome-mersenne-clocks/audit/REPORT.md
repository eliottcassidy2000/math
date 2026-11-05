# Independent audit: palindromic clocks and Mersenne targets

**PASS, no correction.** The packet proves an exact infinite family, a literal fixed-three-port returning operation, and a sharp order bound within the strong feedback-one class. It does not claim unrestricted minimum order or general spectrum coverage.

## 1. Frozen snapshot

Source directory: `session197-prime-atoms/palindrome-mersenne-clocks/`.

- `PROOF.md`: `ad564afb46316206d15a7a4e5bc5fd239d669eff919505a71b2b4b9078c3a9c4`.
- `verify.py`: `a65401af1033f4fb418e1ab296295c936821e0a2448053401befaa435913bc87`.
- Both source receipts: `755fa0fbc2b5de78c43e44b57bdd5e2bd9d6b0890030d2f5901df6320c9a1b19`.
- `witness.json`: `438380fc1fd21ab6d517cc344ade97db4bb4342ea0fce30bd3d4229b8eb867f7`.

The physical integer-gadget recurrence is explicitly inherited from the reviewed ROOTED proof. After the completed review, the author changed only the header phrase “pending independent review” to “independently reviewed”. This report binds that final status-only snapshot; byte comparison confirms no mathematical change. The submitted source files were not edited by this reviewer. All replays were run from exact copies under this audit directory.

## 2. All-order proof review

The two-block identity is valid for every odd positive c, p>=1, q>=0 and t>=0. Reading the central word after the even power-of-two prefix contributes the reverse-complement coefficient with the correct orientation. If trailing zeros exist, their first predecessor is odd and supplies precisely one parity correction. If q=t=0, that correction is instead supplied by the final odd-n term in H=n+f(n)+(n mod 2). Thus the same polynomial includes the empty-trailing-block boundary.

The two physical births really preserve the same three old ports. The left birth loses to a and b and goes immediately after a in the transitive deletion. The right birth beats only the fixed final vertex v and goes immediately before v. The apex b beats both new vertices. The two zero blocks grow by one in each period, and neither named port is replaced. This is a two-module mixed-port clock, not a single-module transitive-clone clock.

The difference H_(t+1)-4H_t is -2B*2^t+3delta. An eventual constant difference therefore forces B=0. The oddness of both c and its reversal first forces p=q+1, then forces literal binary palindromicity. Conversely those conditions cancel the binary mode at every t. The offset inverse is finite and complete: delta+2<=c<=2delta+3 for positive delta, and c must be an odd palindrome. The opposite unit offset occurs only at c=1; delta=1 gives exactly c=3 and c=5. No broader classification of mixed-port clocks follows from this language-specific iff.

Those two positive-unit-offset seeds realize every 2^k-1 for k>=4, with exact binary width k and graph order k+2. Strong feedback-one completeness is correct: once the deletion is put in its unique transitive order, strongness forces the apex to beat its first vertex and lose to its last, leaving exactly the stated binary gadget. The parity-separated bounds on g then exclude every graph of that class with at most k+1 vertices. The regular five-tournament at H=15 is correctly retained as a hostile to an unrestricted minimum-order interpretation.

Finally, 7 divides 2^k-1 exactly when 3 divides k. Factoring at k=3j into (2^j-1)(2^(2j)+2^j+1), both factors exceed one for the relevant j>=2. If the total is 7 times a prime, one factor must equal 7, yielding only j=3 and 511=7*73. This also excludes the repeated-prime case p=7. The proof makes no unproved Mersenne-prime infinitude assumption.

## 3. Independent controls

The standalone `check.py` imports no submitted implementation. It evaluates f by a direct sum over zero positions, constructs graphs with the apex labelled first rather than last, and counts Hamilton words independently by splitting the transitive deletion before and after the apex. Cycle counts use literal anchored word search. Fixed-port controls repeatedly retain the original physical port labels across five periods instead of resetting the ports after each canonical relabelling.

Normal and optimized runs agree on:

- 106,496 exact two-block identities, including nonpalindromic and nonpure cases.
- 2,560 literal fixed-port period returns and every directed incidence under the independent base-order mapping.
- All 1,023 canonical strong feedback-one incidence words through order 12, plus the five named witness graphs, for 1,028 independent physical counts.
- Complete smaller-order checks for exponents 4 through 11.
- 500 complete bounded positive-offset inventories, compared against every odd palindrome below 65,536, together with both unit-offset classifications.
- 1,025 exact Mersenne exponents, including exponents through 4,096.
- Every stored witness arc, its exact H and unique cycle, and all 110 trial divisors in the four prime certificates.

Finite controls test the formulas and code. The arbitrary-order construction, minimum-order conclusion and unique 7p classification follow from the proof, not extrapolation from these controls.

## 4. Reproduction and hashes

Run from this directory:

```
python3 source/verify.py --witness witness-normal.json > source-normal.json
python3 -O source/verify.py --witness witness-optimized.json > source-optimized.json
python3 check.py > independent-normal.json
python3 -O check.py > independent-optimized.json
cmp source-normal.json source-optimized.json
cmp source-normal.json source/receipt-normal.json
cmp witness-normal.json witness-optimized.json
cmp witness-normal.json source/witness.json
cmp independent-normal.json independent-optimized.json
```

- `check.py`: `5277d840859f143ebaee6c541a158ae7c8fb04ae824076141dbb38ede8823fe3`.
- Both independent receipts: `ab567351c6ee1814b0e494cc08e10f75869c630451c70a2a3d0a74d218aa8255`.
- Both copied-source replay receipts: `755fa0fbc2b5de78c43e44b57bdd5e2bd9d6b0890030d2f5901df6320c9a1b19`.
- Both regenerated witness files: `438380fc1fd21ab6d517cc344ade97db4bb4342ea0fce30bd3d4229b8eb867f7`.

The `source/` subtree contains only exact copies of the five named frozen inputs. No external package, network access or repository checkout is needed for reproduction.
