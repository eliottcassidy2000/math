# Physical seed banks and the bounded-root boundary

**Status: PROVED; independent audit PASS in [the audit](../structural_audit/AUDIT.md).** No literature-priority claim. These strengthen Session225 root-descent recovery Section 6 beyond its single seed and unit intercepts. They do not settle arbitrary finite nonunit line banks.

## Inheritance and scope

Source: `05-knowledge/results/session225_root_descent_recovery_20261108/PROOF.md`, Sections 1 and 6; exact O/Z/J formulas are proved in the Session225 graft recovery. Closest mechanism: slope cone plus a missing residue progression. Hostile: the boundary h−g=1 admits every odd count arithmetically. Corrected near miss: arithmetic admissibility is not physical realizability. Least-used sidecar: the injection from Hamilton paths to rooted ordered two-path covers. Concept board: physical path-cut cone; last non-O decomposition; exact line residues; finite arbitrary seed banks. The current extension is internal to Chan–Pak; LRC(14) remains open.

For any nonempty rooted tournament, retain the source meanings D=H(A), E=#Hamilton paths ending at root, U=#ordered spanning two-path covers with second path nonempty and starting at root, d=H(A−root), with empty paths counted once. Operations are

 O(D,E,U,d)=(D,0,2U,d),
 Z(D,E,U,d)=(D−E+U,d,2U−d,d),
 J(D,E,U,d)=(D−E+U,D−E,D+U,D).

An endpoint means applying an arbitrary finite word A and then OJ to a seed. The endpoint chamber has h=D+2U, g=2U, e=h−g=D, where D,U are measured before the final OJ. Its scalar readout is N=h+g=D+4U. As in the source, an **extra closing J** realizes N as a Hamilton-path count. Strong connectivity is not claimed for arbitrary seeds.

## 1. Universal physical cone (recovered lemma)

This is the orientation-dual of the endpoint-defect inequality already used
in `adaptive_state_induction_20261103/BLOCK_SWITCH.md`, B2. We give a
self-contained path-cut proof; the inequality itself is not a new priority claim.

For every nonempty rooted tournament,

 0≤E≤d,   0≤D−E≤U−d,   1≤D≤U.

Proof. Delete the terminal root from a Hamilton path ending there to inject those paths into Hamilton paths of A−root, proving E≤d. Cut any Hamilton path immediately before its root; this injects Hamilton paths into the covers counted by U. Covers whose second path is the root singleton number exactly d; the injected Hamilton paths with that second path number exactly E. Restricting the injection to the complements proves D−E≤U−d. Every tournament has a Hamilton path: insert a new vertex into a Hamilton path of the smaller tournament, before its first out-neighbor preceded by an in-neighbor, or at an appropriate end. Thus D≥1, and D≤U follows from the cut injection. Singleton and root-sink cases with U=d and D=E are allowed; no strict-cone assumption is necessary.

Consequently **every physical seed**, not only t=(1,0,2,1), gives endpoint g/h≥2/3. Also g is even. The operation Z or J never decreases D, since its increment U−E≥0; O preserves D.

## 2. Fixed boundary intercepts form finitely many dyadic sequences

**Theorem.** Fix a finite bank S of physical rooted tournament seeds and an integer C≥1. Among every endpoint obtained from these seeds by arbitrary words followed by OJ, those with e=h−g≤C have scalar counts contained in a finite union of sequences

 N=D+4u·2^k,   k≥0,

where either 1≤D≤C and 1≤u≤2C, or (D,u)=(D_s,U_s) for a seed s with D_s≤C. In particular there are O_{S,C}(log(X+2)) distinct such counts ≤X.

Proof. Write the word before final OJ as A. If A consists only of O's, its state has D=D_s and U=2^k U_s. Otherwise split after its last non-O letter: A=P L O^k, L∈{Z,J}. Let (D_0,E_0,U_0,d_0) precede L and D_1 follow L. Monotonicity of D gives D_1≤C. The common D formula gives U_0=D_1−D_0+E_0≤C, because E_0≤D_0. If L=Z, U_1=2U_0−d_0≤2C; if L=J, U_1=D_0+U_0≤2C, since D_0≤D_1≤C. Physical positivity gives U_1≥1. The trailing O's double U_1 and preserve D_1. Final readout N=D_1+4·2^k U_1 proves the containment and logarithmic count. The argument does not bound word length: arbitrarily long state-preserving Z runs from a root-sink seed are explicitly harmless.

Thus any finite set of lines h−g=c, with arbitrary fixed integer c, contributes only O(log X) counts. Nonpositive c contribute none.

## 3. Finite-seed, finite-line residue-gap obstruction

**Theorem.** Fix any finite physical seed bank and a finite collection of lines a h−b g=c, where a,b are coprime positive integers with a<b and c is any integer, together with any finite set of boundary lines h−g=c and finitely many exceptional counts. Discard the strict-slope lines with a/b<2/3. Put M=lcm(2(a+b)) over the remaining strict-slope lines, taking M=2 for an empty set. If there is an odd residue r modulo M such that

 a r ≠ c (mod 2(a+b))

for every remaining line, the generated scalar counts omit infinitely many positive odd integers. More quantitatively, among 1≤N≤X at least X/M−O(log(X+2)) members of N≡r mod M are missing.

Proof. For a discarded line, g≥2h/3 gives c=a h−b g≤−(2b−3a)h/3, so h≤−3c/(2b−3a). Each contributes finitely many endpoint counts (or none). On every other line, substitute g=N−h and then solve for g:

 g=(aN−c)/(a+b).

Since physical endpoint g is even, aN≡c modulo 2(a+b) is necessary. The selected progression misses every such line. Section 2 bounds all the boundary-line contributions by O(log X); the discarded lines and exceptional counts add O(1). The progression is odd because M is even. Counting its members proves the assertion.

**Concrete corollary (nonunit and arbitrary seeds).** The conclusion holds if every unbounded strict-slope line has |c|<a. Take r=1: then −a<c<a, so c cannot be congruent to a modulo 2(a+b). This includes zero intercepts and nonunit intercepts such as 3h−4g=±2, and it allows arbitrary finite physical seed banks and arbitrary finite boundary intercepts. More generally it suffices that c is not congruent to a modulo 2(a+b), with no size restriction.

The residue test is exact for integer chamber candidates having even g: given N and aN≡c modulo 2(a+b), the formulas above supply integer h and even integer g on the line. It is only a necessary test for physical reachability, and does not assert positivity at small N.

## 4. Hostile: why this does not prove the unrestricted finite-line claim

The five parallel lines 2h−3g=c for c∈{0,−2,−4,−6,−8} cover every odd scalar residue by the exact even-g congruence 2N≡c modulo 10. All have slope 2/3 and each admits positive integer chamber candidates satisfying g/h≥2/3 for all sufficiently large N in its residue. Therefore the residue-gap hypothesis can genuinely fail, including for a finite collection of nonunit lines. No claim is made that these candidates are reachable; no unrestricted finite-nonunit-line impossibility follows from this proof.

The proof strengthens the old mechanism in three orthogonal directions: universal physical seeds; logarithmic control of all bounded e rather than only e=1; and exact even-g congruences handling nonunit and zero intercepts. It leaves arbitrary lines with no residue gap, non-OJ endpoints, and arbitrary nonphysical seed vectors open.

## 5. Verification

`python3 verify.py` and `python3 -O verify.py` should produce identical JSON. The verifier exhaustively obtains D,E,U,d by literal permutations for every labelled tournament through order 5 and every root, including sink roots and singletons; checks the universal cone; tests the exact line congruence against integer/even-g candidates; and exercises the finite dyadic envelope for all words through length 7 from selected physical seeds. Explicit exceptions remain enabled under optimization. These checks support, but do not replace, the all-order proofs.
