# Independent audit: synchronized joint source inverse

**PASS, in the stated physical-source scope.** The synchronized inverse, raw auxiliary-coordinate bound, physical Mersenne-boundary criterion, and exact minimum-root-change formulation are supported by the proof and independent checks below. This does not prove that every admissible target is represented.

## Frozen inputs

- `PROOF.md`: SHA-256 `817ebe7d1fbf810bdbc211ae3f7b01e64d2aed8bdc6c422393c6a80c0f0253e1`.
- `verify.py`: SHA-256 `6400ffb6de25d1162c7937bc016c06a2f3d7f4fcbc81a7d41beab81bedc9a00a`.
- Both author receipts: SHA-256 `6b66dd4658ebd5d29e4fa8e0e501c7db1f14408c470e006a5b68adb4d80e9fed`.
- Physical-adjoint dependency: SHA-256 `f9c337b3cfee1cd1134ee02276890238c30ff0e84a14c26fbcd5cd78ab82dc05`.
- Carmichael-fiber dependency: SHA-256 `40facdf95adc7b97307e53a9058c30538d8ca879be90150797e96f8b4cd85e1e`.

The review requested an explicit odd-target qualification for the simplification `floor(sqrt(2n-6))`; the frozen proof incorporates it. No other mathematical repair was needed. This report approves the proof's mathematics, not a stronger coverage statement or later altered bytes.

## Proof review

1. **One common word is essential.** Sufficiency reconstructs exactly the source `t.alpha beta-dagger` from `v=t.alpha` and `(w,q)=(b,z).beta`. Necessity splits one physical source word and transfers its same suffix to both boundaries. The two bilinear equations therefore refer to one actual source rather than independently chosen witnesses. Minimum J costs must be recorded for the complete synchronized pair.

2. **The auxiliary vector remains bounded without being declared physical.** The extended positive-source cone contains z and is preserved by all three displayed nonnegative coordinate maps. The difference cone `D>=d>=E>=0, U>=d` contains b-z and is also preserved. Thus `0<=q_i<=w_i<=U_w`. The proof does not assume that z, or every transported q, is an actual source. For clarity, the general product bound on the auxiliary readout is `h<=U_v U_q<=2h`; replacing its upper endpoint by `2(h-1)` would require a positive auxiliary deletion coordinate and is not used.

3. **Both physical cover bounds hold simultaneously.** At a split, v and w are actual sources and `U_v U_w<=2(m-1)`. First crossing of sqrt(m-1), the initial split, and the no-crossing terminal split all respect the factor-two growth bound. The seed b has cover3 and m>=5 handles the terminal case. This proves a bound on both physical factors, not merely on their minimum. The same transfer preserves h and J count.

4. **The Mersenne boundary is physically realized.** Direct iteration gives `c_k=b+2(2^k-1)z=t.O^k Z`, so its physicality does not follow from a relaxed cone. The target identity is exact. From `h>m/2`, `n=m+2(2^k-1)h` gives `n>2^k m`; hence the parent's joint inverse fits the stated `C_k(n)`. The bound `m>=5,h>=3` gives `6*2^k<=n+1`, making k finite. The weighted boundary may be large; the proof bounds only its unweighted physical factors. Odd n is now explicit.

5. **Root-change budgets survive in both directions.** Reconstruction adds the two stored costs. Splitting a word preserves their sum. In the complete target decomposition, removing a terminal O and Z-run costs no J; replacing a final J branch by its balanced recipe removes a J. All-Z words use none. Thus the restricted inverse is exact, not only an upper bound on a witness cost.

## Independent hostile controls

- The complete m=15 fiber has h=9 or11, never13. More directly, (15,13) forces D=d=1,U=7 and then E=0, contradicting E congruent to U modulo2. Independent boundary choices nevertheless produce (15,13), showing why mere scalar compatibility is insufficient.
- From the paired seed, OJOJ and ZZOJ reach the same physical w=(27,9,27,9) but different auxiliary q=(15,5,15,5) and (23,7,23,7). With v=t these give (45,25) and (45,39). The author's alternative first word OJJZ gives the same first pair and was checked by its replay. Deduplicating by w alone loses genuine information.

## Reproducible checks

The author verifier was replayed serially under ordinary Python and Python -O. Both replay files are byte-identical to both frozen author receipts. Its complete forward inventory through U=10,000 has 354,227 profiles. The square-root inventories have 1,055 physical and 858 synchronized states through U=200. All 81,843 joint values and their exact minimum J costs agree; the complete C1 criterion and minimum J costs agree through10,001. Seven reconstructed literal graphs pass Hamilton-path and one-cycle checks.

The independent `independent.py` imports none of the author verifier. It separately enumerates source and synchronized inventories, reconstructs every accepted joint presentation, certifies both cone inequalities on generating rays, checks the physical Mersenne identities, and minimizes J costs in increasing-cover topological order. Through m,n<=401 it compares all 674 joint values and 183 C1 targets against a complete forward bank of 2,978 profiles through U=400. The small inventories contain102 physical and76 synchronized states through U=40. Exact minimum J costs match for both joint values and C1 targets. Both independent runs are byte-identical.

From this audit directory, run:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

The forward bounds are complete numeric-cover bounds, not guessed word-length cutoffs. Finite agreement checks the implementation; the all-order inverse depends on the reviewed proof and named physical grammar/adjoint/descent dependencies. No density, singleton-fiber prohibition, signed-response coverage or local solvability is promoted to uniform target coverage.
