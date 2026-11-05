# Independent audit: monotone C1 cutoff and certificate complexity

Verdict: PASS.

Frozen source: session198-c1-monotone.

- PROOF.md SHA256: 6fbdadfe70609af384a6277f02262822f5bb00225d24c2abb2b60a044f71973c
- CERTIFICATE_PROOF.md SHA256: 2253ec47fd5dcb7a381aab20600f02096809e3325e4e635f88b76e1bdd9b2f64
- monotone.py SHA256: a894482c346539c627616c16e7f33197975f5db51abd1a04d109b8e00f55bb1c
- verify.py SHA256: fd2538076611ceaf1c076299198d03f32d4aa275c2dcadece417ca24d60f4874
- Release manifest SHA256: ecc5ef729c2768e2750ca109f9187a67c2b76cc6cb196f583004c6e55deef71b

The readout F=D-E+U gives the exact strong C1 image on nontrivial physical
sources. The closing construction and converse through deletion of a
cycle vertex with one outgoing edge establish both inclusions. Singleton
and nontrivial boundaries are distinguished correctly.

The displayed increments are exact. O raises F by U+E; every indexed M_L
raises F strictly, including L=1 and the singleton input. F(M_L)>=L U
bounds every index in the F<=N search by N/U. Also U<=F and all source
coordinates are bounded by U. Thus the search is finite; strict monotonicity
ensures every prefix of a retained recipe remains within the cap. Profile
merging is legitimate because all successor rules depend only on the profile.

The certificate addendum is valid at all orders. Each move raises F by at
least U>=F/2. A recipe ending with U<=B therefore has O(log B) indexed
moves, each index is at most B, and the recipe takes O(log^2 B) bits.
The gadget recurrences use O(log L)-bit constants and O(log L) recursion
depth. Source and template evaluation is polynomial in log N for the
restricted C1, C2 and C3 scalar languages. This is not a complexity claim
for unrestricted graph counting or unrestricted realizability. The earlier
finite-grammar C2/C3 bounds are stronger but do not invalidate this common
certificate bound. Constructed witnesses have O(log^2 N) vertices.

Isolated normal and optimized runs pass with identical receipts and all
release hashes match. At 50,001 the exact search has 2,074,518 source
profiles, 24,971 positive C1 values and 29 C1 gaps. All selected positive
recipes are independently replayed. All 27 supplied repairs are literally
counted: 26 C2 and one C3, maximum graph order 14. The resulting C1/C2/C3
union has only 7 and 21 missing through this finite bound.

As an additional independent control, independent_three_letter_check.py
uses only the one-vertex O/Z/J physical grammar, not indexed gadget
constants. It obtains exactly the same 68,750 profiles at bound 5,001,
then independently obtains all 2,074,518 profiles and the same complete
positive C1 image and 29 gaps through 50,001. Its receipt is
INDEPENDENT_THREE_LETTER_AUDIT.json.

No finite result is promoted to universal odd coverage or to an
unrestricted negative realization theorem.
