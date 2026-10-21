# Independent twisted-palette audit

Status: proofs APPROVED for the displayed families. All-family exhaustiveness inherits the structural classification and requires its approved source binding. Computational replay pending at initial writing.

The C3/R5 transfer divides ordinary regular words by p and counts reflection-fixed words after removing f constant words. The resulting F=f(q^((p-1)/2)-1) is correctly normalized. For P7, the seven order-three subgroups have pairwise common fixed locus consisting of constants. A twisted-negation-fixed nonregular word belongs to the H0-fixed locus: Ht and H−t for nonzero t imply constancy. That locus has f q words, not f q² or f³. Thus the reflection count f(q³−q), seven reflections, and 21 ordinary representatives yield F=f q(q²−1)/3. Anti-elements of order six cannot stabilize ordinary-regular words because their square would be a nonidentity ordinary stabilizer.

The wreath/palette reduction needs canonical modules and local ordinary regularity; these conditions are correctly stated for strong B and prime X. Q−F counts signed-free ordinary orbits; signed-group orbits halve this. Raw coloring count a(Q−F) has the correct factor.

The q,f recurrences, f/q ratio, strict q growth for c>=3, raw recurrence R_next=R³(1−q^−2), and logarithmic tail all follow directly. The tail bound uses sum_(j>=k)3^(−j−1)=1/(2*3^k). It proves p_k/exp(lambda*3^k) tends to one for c>=3, while the c=2 fixed point correctly retains prefactor 2/sqrt(3). beta positivity follows from the stated uniform lower bound. Conditional signed-freeness and unconditional exponential rarity are different statistics and properly separated.

Automorphism order a_X*3^((n−d)/2), labelled count n!/a, and L(z)=z³/3+z⁵/5+z⁷/21+L(z³/3) are correct for the explicit family. The support has arbitrarily long gaps before later nonzeros (between 7*3^k and 9*3^k). A finite-order polynomial-coefficient recurrence with eventually nonzero leading coefficient cannot have such gaps without becoming eventually zero. Thus both ordinary and labelled generating series are non-D-finite over characteristic zero. No natural-boundary or Hamilton-gap claim follows.

The elementary noncommuting hostile is explicit: binary R5 has palette (6,6), so C3[R5] has (Q,F)=(70,30), whereas R5[C3] has (6,6). Therefore only the latter remains exceptional. This confirms the tower direction cannot be reversed.

Script review: full base automorphism/anti groups are enumerated independently; canonical regular palette representatives and direct twisted reflection incidences agree with distinct formulas. q=0..6 and every feasible f type are covered. Every anti-map is checked for q<=4. No assert-only validity gate. Dynamic/labelled exact checks are algebraic consistency controls, not numerical tests of the infinite-limit proof.

## Replay and frozen version

Normal and optimized replays completed; both match supplied receipts byte-for-byte (9ac0200e62204ab5ae1044bdaaaabd9cbba004edda87751aad97bd753f216916). Verified every MANIFEST.json entry and every SOURCES.json entry. Final document SHA256 737fcf8ac5987b72ecd817a8a3ebc4a4696c4a33ddccd21f421051893a4e43a3 includes the further all-three-seed common asymptotic rate: each ordinary probability is p_k^d*h_X(q_k), the three h_X factors are correct and tend to one, and F/Q tends to zero. Hence asymptotic equivalent exp(lambda_c*n) for both probabilities holds for c>=3, not merely equality of logarithmic rates.

Structural dependency provenance: signed-review/REVIEW.md does approve the structural proof and both finite replay modes, but its last frozen hash is 9a6d9b06ed2e6724c106c5cb101fa8ea95aeeefffc4c0cbb37ed877b0b6b9a1f. The bound current structural source is 675e8943bbf39317af11dce4ab3f64f97dc07375871bdb0ad5e54aca504432bd. Parent asked to confirm intervening approved prose-only changes. This is not a defect in any new palette proof; exhaustive ALL-family wording inherits the dependency approval.

## Dependency resolution and final approval

Read the published classification review's 'Frozen publication revisions' section at math/05-knowledge/results/converse_exception_classification_20261021/session24-signed-review/REVIEW.md. It explicitly approves pre-promotion source 675e8943bbf39317af11dce4ab3f64f97dc07375871bdb0ad5e54aca504432bd after requiring a nontrivial block system. Independently diffed that exact bound source against published STRUCTURE.md (62d039ad51079d719d9a8fb01a2bec17e47788103f24d3a0f7f7971a69f492c7); differences only promote audit status and qualify historical scope, with no proof changes. Thus source lineage is resolved. ALL-family labelled enumeration is approved relative to the classification's stated classical dependencies, not conditional on any unfinished audit.

FINAL VERDICT: APPROVED, no corrections required. APPROVED_SHA256SUMS lists six exact producer files including the manifest. Safe publication may change PROOF CANDIDATE to PROVED relative to the cited and independently audited dependencies, but such document bytes need a refreshed hash.

## Publication status-only promotion

The publication header was promoted from PROOF CANDIDATE to PROVED relative to the independently audited classification and cited classical dependencies, with independently replayed finite controls. No mathematical prose changed. Document SHA256 changed from 737fcf8ac5987b72ecd817a8a3ebc4a4696c4a33ddccd21f421051893a4e43a3 to 1a6b679bf8900c38baa8377cff39ca86931dd9ecc02a2ecb73f25bba53150ae3; MANIFEST.json changed from 9568322216c0b2f8b7cdb566fb1524bcbf06bbd88873e73207390b5a84ffbe42 to 182304f1351ac671ce3e743ad0feebeb91ab83640603d12b96733b82bbd10ae8. The original audit allowlist identifies the pre-promotion lineage; the current allowlist is updated to these publication bytes.
