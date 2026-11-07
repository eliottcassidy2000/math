---
id: THM-4904
title: "Root-change word fibers (Carmichael chamber): for the number M_K(N) of O/Z/J words with at most K root changes and value N from the rooted edge, (1) M_K(N) = [N=1] + 2 A_(K-1)(N) + B_K(N) (final-letter twin wZ/wJ, the analogue of phi(x) = phi(2x)); (2) M_K(N) = [N=1] + Fix_K(N) (mod 2), Fix_K counting dagger-fixed words with exactly K letters J (Session199's adjoint); (3) a singleton value N > 1 has root-change depth exactly K and a unique word wZ with w dagger-fixed, exactly K letters J, (K >= 1) no J at either end and no factor JZ or OJ; (4) the KNOWN height bound W = U - d >= 2^(#O/Z) (Session200); (5) singleton values 1 < N <= X number O_K(sqrt(X) (log X)^floor(K/2)) - density zero, unlike Ford's dichotomy for phi (the K = 0 chamber has Theta(sqrt X) words; that the singletons are also of order sqrt X is HYP-9494 (A))"
status: >
  PROVED: (1)-(3) and (5), given Session199's physical adjoint F(t.w) = F(t.w^dagger) (CITED; re-checked on
  7,174,452 words of length <= 14, and by the audit on 64,570,080 of length <= 16) and the inequalities
  D - E >= d >= 1, D <= U, U >= 2d (Session200; self-contained induction in the audit's Appendix A; re-checked on
  76,333,010 states). New for K >= 1. KNOWN and credited:
  - the K = 0 case of (2), (3) and (5) is Session206's converse involution T (ROOT_OBSTRUCTION sections 1-2;
    THM-4888 (B)): T is the dagger on {O,Z} prefixes, its fixed points form the K = 0 chamber, and its orbit count
    gives O(sqrt X);
  - the J-exclusion in (3) is the argument of Session200 carmichael-fibers Theorem 3, run under a root-change budget;
  - the JZ/OJ exclusion in (3) uses Session205's JZ = OJ (THM-4888 (D)1);
  - (4) is Session200 root-weighted-height section 1, which supersedes this packet's first-draft bound
    U >= 2 c0^(#J+1) 2^(#O/Z).
  The O/J twin and the profile-level Carmichael import are Session200 (singleton profiles {3, 11} to 50,001); the
  strong-class singletons {3, 5, 11, 13} are HYP-9477 (A) / THM-4884 (2b). FINITE-EXACT companions (HYP-9494,
  PROOF section 5):
  - no K = 2 word-singleton above 1 through 10^9 (pre-registered P3 CONFIRMED);
  - F_<=2 misses exactly 140 odd values below 10^9 (P1 CONFIRMED);
  - the strong C = 1 spectrum misses exactly Session197's 29 values below 10^9;
  - M_7(1897) = 1 is the only word-singleton above 1 with K >= 2 below 10^9.
  AUDIT: independent, PASS WITH CORRECTIONS (4 majors: "sharp for K = 0" retyped as HYP-9494 (A); priority of
  Session206 / Session200 credited and Lemma 4 replaced by the KNOWN W-bound; a count corrected to 4,083; HYP-9494 (C)
  restricted to K = 2 after M_7(1897) = 1; 12 minors; MISTAKE-733).
session: kind-pasteur-2026-11-07-S14
date: 2026-11-07
source: 05-knowledge/results/root_change_carmichael_chamber_20261107/PROOF.md
related:
  - HYP-9494
  - HYP-9480
  - THM-4888
  - THM-4884
  - HYP-9477
  - session200/carmichael-fibers/PROOF.md
  - session200/root-weighted-height/PROOF.md
  - session206/root-oz-noncofiniteness/ROOT_OBSTRUCTION.md
  - session199/physical-adjoint/PROOF.md
  - session198/c1-monotone/PROOF.md
  - 05-knowledge/results/session218_smith_repunit_totient_20261107/README.md
script: 05-knowledge/results/root_change_carmichael_chamber_20261107/verify_chamber.c; 05-knowledge/results/root_change_carmichael_chamber_20261107/rc_census.c; 05-knowledge/results/root_change_carmichael_chamber_20261107/literal_witnesses.py; 05-knowledge/results/root_change_carmichael_chamber_20261107/audit/chamber_audit.c
---

# THM-4904 - the Carmichael chamber of root-change words

Full proofs: [PROOF.md](../../05-knowledge/results/root_change_carmichael_chamber_20261107/PROOF.md), sections 1-4.
Audit: [AUDIT.md](../../05-knowledge/results/root_change_carmichael_chamber_20261107/audit/AUDIT.md).
Literature and repository map: [CARMICHAEL_MAP.md](../../05-knowledge/results/root_change_carmichael_chamber_20261107/CARMICHAEL_MAP.md).

## Setting

Words over {O, Z, J} act on t = (D,E,U,d) = (1,0,2,1) by the following moves. Each is a literal vertex birth, and J
is a root change.

    O = (D,0,2U,d),   Z = (D-E+U,d,2U-d,d),   J = (D-E+U,D-E,D+U,D)

The value of a word is its final D, a Hamiltonian-path count. F(s) = D - E + U is the value of the word extended by Z,
or equally by J.

M_K(N) counts the words with at most K letters J, ending in Z or J, of value N, plus the empty word for N = 1. This is
the multiplicity of the S13 root-change census. Two auxiliary counts:
- A_j(N) = #{w : #J(w) <= j, F(t.w) = N};
- B_K(N) = #{w : #J(w) = K, F(t.w) = N}.

The dagger w^dagger reverses w and exchanges O with Z, fixing J.

## Statement

1. **Twin decomposition.** M_K(N) = [N = 1] + 2 A_(K-1)(N) + B_K(N).
2. **Parity law.** M_K(N) = [N = 1] + Fix_K(N) (mod 2), where Fix_K(N) counts dagger-fixed words w with #J(w) = K and
   F(t.w) = N. The dagger-fixed words are A c A^dagger with c in {empty, J}.
3. **Chamber.** Suppose N > 1 and M_K(N) = 1. Then:
   - N is not in F_<=K-1, so N has root-change depth exactly K;
   - the unique word is wZ with w = w^dagger and #J(w) = K;
   - for K >= 1, w does not begin or end with J, and contains neither of the factors JZ and OJ;
   - for K = 2, w = B J C C^dagger J B^dagger with B nonempty and ending in Z, and C empty or beginning with O.
4. **Height bound (KNOWN, Session200).** W = U - d satisfies W(t.w) >= 2^(n(w)), where n(w) counts the letters O and Z.
   Hence F(t.w) >= 2^(n(w)) + 2.
5. **Density.** #{1 < N <= X : M_K(N) = 1} is at most the number of chamber words w with F(t.w) <= X, which is
   O_K(sqrt(X) (log X)^floor(K/2)). For K = 0 the chamber has Theta(sqrt X) words.

## Reading

The theorem transfers the shape of the classical Carmichael descent, not its numbers.

- (1) plays the role of phi(x) = phi(2x) for odd x.
- (2) and (3) play the role of the Carmichael-Klee lemma: they force a singleton into a small chamber.
- (5) is the contrast. Ford (Ramanujan J. 2 (1998), Theorem 5) proved that counterexamples to Carmichael's
  conjecture are either absent or a positive proportion of totients. Here the singletons are confined unconditionally
  to a square-root-size chamber.

For K = 0 the chamber itself has Theta(sqrt X) words. The observed singletons are:
- 500 of the 1,023 chamber words below 2^20;
- 4,083 of the 8,191 chamber words below 2^26, that is 0.498 sqrt(X).

See HYP-9494 (A). For K = 2 no singleton above 1 occurs through 10^9, which is HYP-9494 (C). For larger K, singletons
can sit at values of exact depth K: M_7(1897) = 1.

## Not claimed

The theorem does not claim:
- that K = 2 has no singletons for all N (HYP-9494 (C));
- a lower bound for the singletons;
- any statement about the numerical source profiles of Session200, or about tournament isomorphism classes, which are
  different fiber notions;
- any coverage theorem or the Chan-Pak conjecture.
