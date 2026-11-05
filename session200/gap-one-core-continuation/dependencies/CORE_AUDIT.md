# Independent audit: canonical sink cores and deficit fibers

**PASS after two scope clarifications, now incorporated in the frozen proof.** The canonical sink-core projection, complete predecessor-bounded core inventory, exact joint and C1 compilers, and density-zero conclusion for the low-deficit C1 family are sound. This does not establish density zero for the whole C1 spectrum, obstruct C2/C3 repairs, or prove universal coverage.

## Frozen inputs

- `PROOF.md`: `04f7fadf15bf302c3c4f8eefa4dee392f3051a0a22d872b204a96fb51d32e3fd`.
- `README.md`: `bc1933224deba3a26dd8a7a901c3091763638275e50225894097d0238df223db`.
- `verify.py`: `56cf3c2a0254efdf2d8da67474083a478d6e5d56b14cb0db7621e26ecd032d28`.
- Both author receipts: `1f9e4a0dd1ebf346dfbda10dd545fcfbedde7a8451efafbb348a66ce250c3b3c`.
- All input and dependency hashes are recorded in `INPUTS.json`.

## Scope corrections made during review

1. The inherited positive outer-exponent descent requires the scalar `F(x)` to lie outside `B={3} union G(S)`. The condition is not merely `D_x!=U_x`. The actual source `x=t.ZJ=(5,2,6,3)` has unequal D and U but F(x)=9 in B; applying the unqualified reconstruction would give the impossible parent projection `(D,U,d)=(3,3,3)`. The final proof states the precise scalar-image hypothesis in both relevant places and includes this hostile.

2. The sink exponent and core triple are canonical for a supplied source projection, not for a joint pair alone. The words ZZZJJ and ZZJJJ yield profiles `(27,9,35,17)` and `(29,10,34,15)`, both with `(m,h)=(97,53)`, but distinct core triples. The final opening and Theorem4 now explicitly preserve this distinction. Both regressions are tested.

These changes restrict ambiguous prose to the hypotheses actually proved. They do not change the underlying core construction.

## Mathematical audit

**Endpoint/height boundary.** A nonseed zero-endpoint source must end in O and has U>=2D. Z outputs satisfy `2D'-U'=2(D-E)+d>=3`; J outputs satisfy `2D'-U'=D+U-2E>=3`, with D=1 handled separately. Thus a positive-endpoint core has `D<=u<=2D-3`. The forbidden two-height gap is genuine, not a relaxed cone inference.

**Canonical exponent and seed.** Stripping a maximal O suffix yields a positive-endpoint core unless it reaches t. For D>1 the core interval `[D,2D)` selects exactly one dyadic rescaling, giving the displayed floor-logarithm. This argument supplies exact divisibility; it does not presume an arbitrary U is divisible by the chosen power of two. When D=1, Z and J strictly increase D, so only t.O^a is possible. Its U is `2^(a+1)`, and a is one less than log2(U), not the generic floor-logarithm formula. The proof handles this exception correctly. Endpoint E and graph identity are not claimed canonical.

**Complete finite core bank.** Every nonseed core ends in Z or J. The predecessor satisfies `U_p<=D_q-1<=R-2`. Conversely both constructors are physical and positive-ended. The family t.O^j.Z attains predecessor height R-2 and output height2R-5. The inventory therefore counts actual projected triples, merging only when D,u,d agree and retaining a recipe. The cubic expression `(R^3-R)/6` is a coarse upper bound, not an exact count. Its counting proof legitimately forgets E, overcounts possible d and u, and separately covers the exceptional seed at D=1.

**Exact joint and target normalization.** Adding O preserves D,d, doubles U and costs no J, so `h=2^(a+1)u-d` is exact. Retaining the minimum cost over all realizations of each core triple is valid for G/H; minima are then taken over all triples and exponents representing the requested pair. The source recipe for a C1 target is exactly the stored core word followed by O^(a+1), Z^k, and the closing J. Its initial O erases the forgotten endpoint. The inherited terminal-run decomposition proves completeness, including k=0 and the separate value3. The double canonical descent for F(x) outside B uses the repaired scope above.

**Density statement.** Every term satisfies `h>=2^(a+1)` and `M_k>=2^k`, so n<=N bounds both exponents by floor(log2 N). Multiplying the finite triple bound by the number of exponent choices proves the stated count. For fixed R it is O(log^2 N); for `R(N)=o((N/log^2 N)^(1/3))` it is o(N). This concerns A_R, the family represented with that deficit budget. The union over unbounded deficits is the complete C1 image and is not shown sparse by this argument. The statement about inability to cover a positive fraction of odd targets with a small deficit budget is valid. No conclusion about the full C1/C2/C3 spectrum is licensed.

## Exact reproducibility

The author verifier was replayed serially in normal Python and under Python -O at cap 5,001. Both outputs match both frozen receipts byte-for-byte, and the source manifest passes. It uses 126,910 profiles through U=5,000, yielding 27,198 **cumulative core triples with deficit<=2,500**, 29,322 joint pairs and 2,471 C1 targets through 5,001. The reported 42,447 parameter choices are presentations, not distinct target values. Minimum J costs match the forward inventory. Twenty-nine scalar-nonbalanced profiles pass canonical double descent.

The independent `independent.py` imports no author code. It separately builds a predecessor-bounded core bank and compares it against a full positive-endpoint inventory through the proved output-height bound. Cumulative core counts at R=2,3,4,5,6,8,16,32,64,128 are respectively 1,1,2,2,4,5,15,41,119,334. This checks projected triples and exact minimum J costs, not word multiplicities.

At joint/target cap 1,001, the independent forward bank has 11,721 profiles through U=1,000. A core bank of 2,505 triples at deficit<=500 yields exactly 2,593 joint values and 475 C1 targets, with matching minimum J costs. Every generated arithmetic choice is checked against its literal source word. All 9,841 words through length 8 pass the canonical suffix tests, including non-minimal recipes. Sixteen small graph controls from the four displayed core rays were independently built and counted by exhaustive vertex permutations; all have the stated Hamilton-path count and exactly one Hamilton cycle. Normal and optimized independent receipts agree byte-for-byte.

For additional count controls, through 100,001 the core-budget families at R=2,4,6,16,32 contain 66,170,334,1,007,2,313 distinct targets, respectively. These are finite low-deficit-family counts, not whole-spectrum counts or evidence by themselves for the asymptotic theorem.

From this audit directory:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

The all-order conclusions rest on the reviewed proof and its named physical dependencies. Finite inventory agreement verifies the implementation. Universal coverage and the all-order singleton-fiber conjecture remain OPEN.
