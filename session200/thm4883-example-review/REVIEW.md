# A scoped correction to the J-born example in THM-4883

**Status: VERIFIED COUNTEREXAMPLE to one illustrative example.** The exact T4 formula and its gcd qualification are preserved. No other THM-4883 statement is audited or refuted by this packet. Original files and their ownership are unchanged.

## 1. Exact claim and source snapshot

At main commit `45766d058d9507f97e87336c09a15cd3ef65796d`, the example occurs in two places:

- [THM-4883, T4, line 64](https://github.com/eliottcassidy2000/math/blob/45766d058d9507f97e87336c09a15cd3ef65796d/01-canon/theorems/THM-4883-signed-one-apex-responses-exact-2-adic-structure.md#L64): “any J-born source over a strong source”.
- [The proof packet, Section 7, lines 271–273](https://github.com/eliottcassidy2000/math/blob/45766d058d9507f97e87336c09a15cd3ef65796d/05-knowledge/results/signed_response_2adic_structure_20261105/PROOF.md#L271-L273): “every J-born source over a strong source” is listed as having d_v=E_v+1 and yielding all odd offsets.

The source blobs are respectively `c1d0e5efd5a242a156e16a87fd4f7d4645e28e88` and `39a3e9c0c226e06b2f61235ee07fcd25124c2349`. Named-file histories end at the original pre-audit checkpoint `c4e9654f72295da4562f1f94bb193f6cb8754989`. A commit search for 4883 found no later relevant correction, and the referenced `audit/AUDIT.md` returned 404 at the pinned main snapshot. This is a bounded freshness finding, not a claim about future revisions or every possible note elsewhere.

## 2. The first failed implication

For an actual rooted start-one source s=(D,E,U,d), the physical J operation is

    J(s)=(D-E+U,D-E,D+U,D).

Consequently

    d(Js)-E(Js)=D-(D-E)=E(s).

Strongness of s gives E(s)>=1, not E(s)=1. Replacing that endpoint count by one is the error. The deletion and endpoint coordinates must both be retained.

## 3. Actual graph counterexample

Start with the transitive-edge source t=(1,0,2,1). Chronological source words give

    u=t.ZOJ=(9,3,9,3),
    v=u.J=(15,6,18,9).

Both are literal strong tournaments with exactly one Hamilton cycle and exactly one path starting at the designated root. The verifier constructs their adjacency rows and independently counts the full profiles. Here v is J-born over the strong source u, but d_v-E_v=3.

For T4's coefficients,

    alpha=D_v+4U_v-2d_v=69,
    beta=alpha+E_v-d_v=66,
    c=-D_v+U_v+d_v+2E_v=24,
    gcd(alpha,d_v-E_v)=3.

For every strong unit source q=(D,1,U,1), the Hamilton-path count D is odd. Therefore

    M(v,Jq)-2N(v,Oq)=69D-66U+24 = 3 (mod 6).

Indeed, writing D=2k+1 gives the identity

    69(2k+1)-66U+24=6(23k-11U+15)+3.

Offsets +1 and -1 are impossible. This is a universal obstruction for this fixed actual left source, not merely a failed bounded search; no parity assumption on U is needed.

Taking the triangle q=t.Z=(3,1,3,1) supplies a small literal control:

    N(v,Oq)=417, with exactly two Hamilton cycles, order 11;
    M(v,Jq)=867, with exactly two Hamilton cycles, order 12;
    867-2*417=33.

Both graphs and all source rows are recorded in the verification receipts. The source graph strongness is certified by its Hamilton cycle, rather than assumed from an unrestricted numerical profile.

## 4. Strongest survivor and repair

The full T4 algebraic identity remains correct. The verifier checks it on a complete coefficient basis, independently of the literal counterexample. Its gcd-qualified conclusion predicts exactly the observed class modulo six, so the counterexample does not refute that conclusion.

For this witness the signed readout is

    24-3(22U-23D).

The previously audited [even-numerator response theorem, Sections 1–3](https://github.com/eliottcassidy2000/math/blob/45766d058d9507f97e87336c09a15cd3ef65796d/session198/rational-response-coverage/PROOF.md#L10-L151), applied at slope 22/23, makes 22U-23D range over all odd integers on actual strong unit sources. Thus this example's offset image is exactly 3+6Z, consistent with T4's qualification. No new number-theoretic conjecture is needed.

A sufficient replacement for the false illustrative phrase is:

    J-born sources whose predecessor has E=1.

Alternatively, retain the exact general condition

    gcd(D_v+4U_v-2d_v,d_v-E_v)=1

in the E_v<d_v branch. The condition d_v-E_v=1 remains a valid sufficient special case. There is no claim here that the positive parent count can be independently prescribed.

## 5. Independent checks and scope

`verify.py` uses separately implemented physical graph births, N/M grafts and subset dynamic programming for Hamilton paths, endpoint counts, reserved two-path covers, deletion counts and Hamilton cycles. It reproduces u, v, 417, 867 and the offset 33. It also checks the correct J-gap identity and the general T4 offset identity on formal coordinate bases. Normal and optimized receipts are byte-identical.

A second implementation, `audit/independent.py`, reads only the recorded adjacency rows. It enumerates permutations for the two complete rooted-source profiles and directly traverses all simple directed paths for the two larger grafts, without subset DP, memoization or imports from `verify.py`. It checks tournament adjacency, recomputes the same Hamilton counts, and normalizes cycles by all possible starting vertices. Its normal and optimized receipts also agree.

This packet records the exact false implication, actual witness, missing coordinate and valid repair. It does not modify the original researcher's files, claim ownership of the surrounding theorem, or provide a full audit of its 2-adic classification, half-line, clock-uniqueness or witness-size claims. Session 200's previously audited signed-fiber and canonical-descent results do not use the false example.
