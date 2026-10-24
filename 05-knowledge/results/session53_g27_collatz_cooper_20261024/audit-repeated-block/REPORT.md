# Separate audit: repeated 27 block and OPEN tower appendix

Verdict: the repeated-block mathematical statements pass independent proof review. The tower remains an OPEN higher-depth hypothesis, with its certificate treated only as finite evidence. This audit does not modify or supersede the earlier birth/return audit.

## Repeated-block proof

For the repeated word w^k, composing F(n)=(an+b)/q gives numerator a^k n+b(q^k-a^k)/(q-a). Its odd source coefficient gives one integrality residue modulo q^k. The parity-cylinder theorem supplies a realized word in exactly one residue class; necessity of this integrality congruence identifies the classes. Thus the argument does not assume that arbitrary unguarded rational iterates are native Collatz paths.

With d=q-a odd, the residue is b/d modulo q^k. The rational fixed point is reduced and strictly between 8 and 9. It is not an integer, so v2(dn-b) is finite for every positive integer n. Membership in the nested cylinders gives K(n)=floor(v2(dn-b)/59) exactly. For K>=1 the identity dF(n)-b=(a/q)(dn-b) lowers this valuation by exactly 59 and the fuel by one. No decrease is claimed for K=0 or for unrestricted Collatz steps. The inherited pumping note, section 5, already states the general binary-fuel mechanism; the packet correctly credits it rather than asserting novelty.

Every block source lies in the original positive 27 residue class and hence is at least 27. The prior source-relative slope certificate applies afresh at that source. Consequently every prescribed block has first shortcut descent at relative time 59, ordinary descent at 96, and the complete prefix uses 96k ordinary operations. This is not a total stopping formula. Inserted ordinary odd peaks exceed their preceding values and cannot create an earlier descent.

Modulo 7, F is 1-x, so the coefficient reduction has image C2. The real slope a/q is strictly between zero and one, while the 2-adic norm is 2^59. Its square is not the native identity. Every guarded integer source exceeds the nonintegral fixed point, so its endpoint strictly decreases even when its residue returns. Nested finite cylinders have no common positive integer; their unique 2-adic intersection does not imply an integer cycle.

Each cylinder is one arithmetic progression. Subtracting adjacent progression counts gives the geometric initial-fuel density, and CRT with the odd modulus 7 gives the joint density divided by seven. Each progression's discrepancy is at most one, so the difference discrepancy is at most two, uniformly in level and interval. This is an initial-word distribution under natural density, not a prime distribution or convergence prediction. Independent exact rational checks cover interval boundaries at large cylinder representatives and X=2^2000+7, rather than relying only on the producer's small-X enumeration.

## Code and package checks

The producer check.py computes closed affine residues and independently reconstructs them by one-bit parity lifting. Its direct ordinary replay checks each fresh-source first descent and 96-step block; its prefix comparison tests exact exhaustion of K rather than only the lower bound. Normal and optimized modes must reproduce both saved fixtures. The audit verifies every producer manifest hash and copies only declared payloads plus the manifest into a fresh work directory.

Source mutations test the period exponent, affine carry, fuel denominator, native-involution confusion, clock confusion, fixed-point interval and parity of the quotient. Explicit exception gates remain effective under Python -O. A genuine fixture-coverage gap was found: changing the fuel divisor from 59 to 58 survives the producer's bounded samples because their tested valuations do not distinguish the two quotients. Independent sources constructed with exact valuations 58, 116, 117, and higher boundaries distinguish the wrong formula through direct parity replay. This is a test-coverage limitation, not a defect in the theorem or the stated correct code. Changing the divisor to 60 is rejected by the producer checks. Detailed outcomes, including the surviving mutation and independent rejection, are in audit_results.json. A common additive offset in both arithmetic-progression counts cancels in their difference and is not a semantic mutation of the claimed distribution; it is recorded separately rather than mislabeled as a rejected hostile.

## Separate OPEN two-step tower review

The tower is stated as c(p^r n)=c(p^(r-2)n) modulo p^(2r-2), r>=2. The target at odd depth is correctly c(pn), not c(n). Substituting n=p^a m converts the same difference to depth r+a, so the strengthened valuation checked for p-divisible n is the correct consequence of the hypothesis, not extra proved precision.

verify.py generates exact integers to coefficient index 20000, checking every recurrence division, and enumerates all admissible primes, depths and positive n. Valuations are computed from nonzero integer differences, with nonzero normalized residues retaining sharpness witnesses. These are overlapping finite comparisons, not independent samples. The recorded 15,407 comparisons do not establish the unbounded hypothesis.

The proof text correctly limits the conditional base result to modulus p². Replacing n by p^(r-2)n does not upgrade precision. Its higher correction and Cartier bounds remain missing. The discussion of split/ramified one-step conjectures correctly identifies the two-step version there as a consequence, not a stronger new conjecture. Public-search novelty is not certified by this audit.

The tower's own replay verifies exact recursive manifest membership, file hashes, normal and optimized fixtures; its hostile suite rejects changed proof bytes, an omitted fixture and an undeclared nested file. This audit reruns those checks on a manifest-only relocated copy and preserves OPEN in its results.

## Reproduction

The audit MANIFEST.json lists the closed payload set excluding itself. Inputs are frozen producer files, including producer manifests. Run python3 audit.py or python3 -O audit.py from any working directory. Generated work/ and mutants/ directories are disposable and are not publication inputs. No external packages or repository files are required. No repository writes were performed.
