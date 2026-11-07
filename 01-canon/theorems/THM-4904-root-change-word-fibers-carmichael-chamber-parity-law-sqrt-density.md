---
id: THM-4904
title: "Root-change word fibers (Carmichael chamber): for the number M_K(N) of O/Z/J words with at most K root changes and value N from the rooted edge, (1) M_K(N) = [N=1] + 2 A_(K-1)(N) + B_K(N) (final-letter twin wZ/wJ, the analogue of phi(x) = phi(2x)); (2) M_K(N) = [N=1] + Fix_K(N) (mod 2), Fix_K counting dagger-fixed words with exactly K letters J (Session199's adjoint); (3) a singleton value N > 1 has a unique word wZ with w dagger-fixed, exactly K letters J, N outside F_<=K-1, and (K >= 1) w neither beginning nor ending with J (Session200's O/J twin); (4) growth lemma U(t.w) >= 2 c0^(#J+1) 2^(#O/Z), c0 = prod (1 - (1/4)(2/3)^i) = 0.44257; (5) singleton values number O_K(sqrt(X) (log X)^floor(K/2)), sharp for K = 0 - density zero, unlike Ford's dichotomy for phi"
status: >
  PROVED (new): (1)-(5), given Session199's physical adjoint F(t.w) = F(t.w^dagger) and Session200's physical
  inequalities D - E >= d >= 1, D <= U, U >= 2d (both CITED, both re-verified exhaustively here: 7,174,452 words of
  length <= 14, 76,333,010 states). KNOWN / CITED: the O/J twin and the Carmichael import at the level of source
  profiles (Session200 carmichael-fibers, singleton profiles {3, 11} to 50,001); strong-class singletons {3, 5, 11, 13}
  (HYP-9477 (A), THM-4884 (2b)). FINITE-EXACT companions (HYP-9494, PROOF section 5): no K = 2 word-singleton above 1
  through 10^9 (pre-registered P3 CONFIRMED); F_<=2 misses exactly 140 odd values below 10^9 (P1 CONFIRMED); the strong
  C = 1 spectrum misses exactly Session197's 29 values below 10^9. AUDIT: pending.
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
  - session199/physical-adjoint/PROOF.md
  - 05-knowledge/results/session218_smith_repunit_totient_20261107/README.md
script: 05-knowledge/results/root_change_carmichael_chamber_20261107/verify_chamber.c; 05-knowledge/results/root_change_carmichael_chamber_20261107/rc_census.c; 05-knowledge/results/root_change_carmichael_chamber_20261107/literal_witnesses.py
---

# THM-4904 - the Carmichael chamber of root-change words

Full proofs: [PROOF.md](../../05-knowledge/results/root_change_carmichael_chamber_20261107/PROOF.md), sections 1-4.
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
   - N is not in F_<=K-1;
   - the unique word is wZ with w = w^dagger and #J(w) = K;
   - for K >= 1, w does not begin or end with J;
   - for K = 2, w = B J C C^dagger J B^dagger with B nonempty.
4. **Growth lemma.** For every word w, U(t.w) >= 2 c0^(#J(w)+1) 2^(n(w)), where n(w) counts the letters O and Z and
   c0 = prod_(i>=0) (1 - (1/4)(2/3)^i) = 0.442573...
5. **Density.** #{N <= X : M_K(N) = 1} <= #{chamber words of value <= X} = O_K(sqrt(X) (log X)^floor(K/2)).

## Reading

The theorem transfers the shape of the classical Carmichael descent, not its numbers.

- (1) plays the role of phi(x) = phi(2x) for odd x.
- (2) and (3) play the role of the Carmichael-Klee lemma: they force a singleton into a small chamber.
- (5) is the contrast. Ford (Ramanujan J. 2 (1998), Theorem 5) proved that counterexamples to Carmichael's
  conjecture are either absent or a positive proportion of totients. Here the singletons are confined unconditionally
  to a square-root-size chamber.

The bound is sharp for K = 0: about 0.49 sqrt(X) singletons, about half of the 1,023 = 2^10 - 1 chamber words below
2^20. For K = 2 no singleton above 1 occurs through 10^9 (HYP-9494 (C)).

## Not claimed

The theorem does not claim:
- that K = 2 has no singletons for all N (HYP-9494 (C));
- any statement about the numerical source profiles of Session200, or about tournament isomorphism classes, which are
  different fiber notions;
- any coverage theorem or the Chan-Pak conjecture.
