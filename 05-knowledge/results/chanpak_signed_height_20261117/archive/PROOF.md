# Exact height is a unimodular source-fiber problem

**PROVED arithmetic reduction and physical padding obstruction; PROVED exact membership for a changed-source grammar. Coverage remains OPEN.**

## Archive walk: what was actually imported

This pass read tangential sources rather than only recent Chan–Pak constructions:

- June30, `covering_min_sternbrocot_analysis_macmini_20260630.out`: the exact ray k/((n-1)k+1) suggests retaining a denominator and a determinant-one relation, not merely a limiting slope. Its historical assertions about the LRC frontier are not used.
- June30, `07-reflections/the-descent-is-a-branched-cover-the-residual-is-its-ramification.md`: the reflection itself distinguishes an exact cover computation from an unproved identification with the LRC residual. The useful methodological distinction is between a genuine pullback and a similarity of shapes; no geometric identification is imported.
- August15, `07-reflections/four-point-turnpike-xor-gram-branch-language-and-harmonic-sidecars-codex-20260815.md`: the order stratum forgets the marked origin. Here a branching reader language similarly does not choose an absolute source height.
- September25, `05-knowledge/results/duck_zeckendorf_20260925.md`, §§1–2: a two-coordinate charge and a normal form carry more information than the represented scalar. This suggested searching for a primitive integral second coordinate rather than a real approximation invariant.
- October14, `05-knowledge/results/aliquot_dyadic_frontier_20261014e.md`, §§1–2: fixing the cofactor makes the remaining candidate an exact affine solve; primality and native-source constraints stay separate. Here reader parameters can likewise be eliminated before testing physical source membership.

All paths above are relative to the repository's `math/` checkout in this workspace. None of the old conjectural cross-problem identifications is assumed. The current physical inputs are `pubsub_branch/PROOF.md` and the extreme-back-edge slot formula in `pubsub/observer/PROOF.md`.

## 1. Determinant one removes the arithmetic denominator

Let t be any legal branching-reader state (called h in the preceding branching note). Its final positive coefficients are

    a=12t-3, b=15t-4,
    N=aH+bP,                 P=P₂(source).

They satisfy5a-4b=1. Introduce the source coordinate

    J=4H+5P.

Then the transformation from (H,P) to (N,J) is unimodular, with exact inverse

    H=5N-bJ,    P=aJ-4N.

There is no lost residue or rational denominator. Equivalently,

    N=(3t-1)J+(H+P).

The old source switch (ΔH,ΔP)=(10,-8) fixes J and changes N by2. More generally, along a fixed J fiber the integer step(5,-4) changes N by1. The odd-H constraint selects every other step, exactly matching the native even count gaps. Thus the existing +2 is a **primitive source-fiber motion**, not merely an accidental cancellation of large coefficients.

For example the inherited source(187,3736) has J=19428. At t=1 it reads42779; the switched source(197,3728) has the same J and reads42781.

## 2. Exact target intervals and physically necessary congruences

For an even-order source of order at least4, necessary conditions are H positive odd, P>=4H, 8|P, and P=H mod3. The inequalities and divisibilities follow from the inherited odd-cycle cover expansion: the number of cycle components has even parity; each nontrivial cycle contributes an extra factor2. They are not sufficient for physical realizability.

For odd N and odd legal t, the inverse above satisfies all these conditions exactly when

    ceil(24N/(a+4b)) <= J <= floor((5N-1)/b),
    J=12 mod24.

Indeed P>=4H gives the lower endpoint and H>=1 gives the upper. Since a=1 mod8, 8|P means J=4 mod8; P=H mod3 means3|J. These combine to J=12 mod24 and force H odd.

The real interval has length N/[b(a+4b)]-1/b. Therefore for every odd N>=25b(a+4b) it contains such a congruence class representative. This proves that the elementary local and cone constraints alone leave target candidates at every sufficiently large height for each fixed reader t. It does **not** realize those candidates as tournament profiles. Excluding the known impossible source values H=7,21 is another guard, not the missing full classification.

This gives a precise remaining sufficient condition: find a physically attained (H,P) among the inverse profiles on one of these target-owned J intervals, for a legal t. Unlike the original branching language, this condition changes the source with the target.

## 3. Why trivial padding cannot keep the +2

Padding the original switched source with s singleton SCCs preserves ΔH=10 but multiplies ΔP by q=2^s. A general positive single-return rooted reader reads hH+RP. Its terminal cover coordinate always satisfies R>=h: cut each full Hamiltonian path immediately after its distinguished root; the prefix ends at the root and the remaining suffix is the optional second path. Concatenation recovers the path, so this is an injection into the covers counted by R.

Consequently for s>=1,

    ΔN=10h-8qR <= (10-8q)h <= -6h <= -6.

Neither orientation yields a gap of2. This obstruction allows an arbitrary change of positive single-return reader, not just the old branching invariant. To preserve a native +2, changing the source must do more than singleton padding, or the subscriber must leave the positive two-coordinate architecture.

## 4. A changed-source grammar with exact two-parameter inversion

There is nevertheless a concrete source-changing construction whose absolute membership is easy. Use the strong extreme-back-edge source E_m, m>=1, with

    H_m=2^m+1,
    P_m=4*2^m+2*(3^m-1).

Permit s>=0 singleton padding vertices and any legal branching state t. Substitution into its single-return reader yields

    N=(12t-3)H_m+(15t-4)2^s P_m.

The final tournament is strong under the same nonempty uniform single-return construction, even if the padded input is not strong. For a prescribed N and fixed(m,s), the only possible reader state is

    t=(N+3H_m+4*2^s P_m)/(12H_m+15*2^s P_m).

Check exact divisibility, positivity, and the complete branching-word inverse. There are only O(log² N) pairs(m,s) to test: every legal t>=1 gives N>=9H_m+11*2^sP_m. Thus this is a polynomial-bit-time exact target decoder for an independently changing physical source and reader, not a scan over all words or a heuristic fit.

With s=0 the formula simplifies to

    N=(72t-19)2^m+(30t-8)3^m-18t+5.

This combines the archive's fixed-core affine elimination with an actual two-base physical source family. It escapes the old fixed source's single residue progression, but the source switch companion is not retained. For instance the empty reader word t=1 gives N=159 at m=1 and N=397 at m=2. No all-target coverage or +2 for this changed-source grammar is claimed.

## Replay

`python height_select/archive/verify.py` checks the unimodular inverse, all necessary target gates on tested intervals, the old +2 fiber, and exact decoding of hundreds of independently generated changed-source outputs, including long branching words. It also recomputes the small extreme-back-edge source profiles directly. Finite tests supplement the all-order arguments above.
