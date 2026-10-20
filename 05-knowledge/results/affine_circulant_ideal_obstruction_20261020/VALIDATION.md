# Independent audit and completed-run receipt

2026-10-20. **VERIFIED.** This packet is ready for publication in
`05-knowledge/results/affine_circulant_ideal_obstruction_20261020/`.
All mathematics is scoped in `AFFINE_CIRCULANT_IDEAL_OBSTRUCTION.md`.
No historical priority claim is made.

## Proof audit

- The source is an affine torsor, not a vector space. Complement invariance follows from T^m; quotienting by complementation is exactly the odd augmentation hyperplane. Its affine-to-linear extension is unique because that hyperplane spans the whole word space.
- The commuting algebra is the regular cyclic algebra F₂[x]/(x^m−1), including repeated factors. No diagonalizability assumption is used when m is even.
- The ideal obstruction uses the all-one socle word J=M/(x+1). An ideal missing the full (x+1)-primary factor must contain J, which has odd primitive-root weight.
- The adjacent-pair hostile has least period m only for m>2. The m=1 and m=2 exceptions are explicitly retained.
- A type-surjective equivariant map is labelled-surjective because its image is a union of orbits. A type bijection is not a labelled bijection; the positive odd-m family has two-element labelled fibers.
- All inputs to h lie in the odd augmentation hyperplane, hence have unit first primary component. This fixes the output primary height. Powers-of-two heights are exactly the bad heights, witnessed by killing all other primary components.
- The shell orbit formula uses fixed words of every divisor period, not just primitive words. Their first primary block has length 2^v₂(d); the zero shell and positive shells have respectively 1 and 2^(k−1) choices there.
- The sharp coverage addition compares divisor contributions separately for every odd r|q. For e≥4, Σ_(E₀|E|e)φ(e/E)=e/E₀ and e/E₀≤2^(e−k−1); the zero-shell bound is e≤2^(e−2). Thus the maximum is A(m)/2, not merely a numerical pattern.
- The graph/tournament interpretation uses the inherited prime-order Sylow/normalizer dictionary. The all-m necklace theorem is not promoted to a composite-order graph-isomorphism theorem.

No mathematical error was found in the inherited staged argument. The initial
`verification.json` was stale: its check count and source hash predated 60
additional checks, although every data row agreed. It was replaced by fresh
normal and optimized runs. The final source adds another 39 checks for the
closed coverage law. Historical `audit.*` scratch receipts are not part of the
publication inventory.

## Exact computation

Main result: PASS, 45,248 explicit checks; every h and every input word
for m=1 through 12. This covers 22,369,620 source-map evaluations
and 11,184,810 odd-word evaluations. Every arbitrary affine
matrix/translation is independently enumerated at m≤3 inside this checker.

The supplementary coefficient-level audit independently derives
AR=RA and A(e₀)=b+R(b) from T(f)=R(f)+e₀. It covers every matrix and translation
through m=4, including 1,048,576 parameter pairs at m=4, and obtains exactly
2^m equivariant maps. Matrices failing AR=RA are discarded for every b by the
proved coefficient identity, not by a heuristic filter.

Normal and optimized main JSON outputs are byte-identical. The supplementary
normal and optimized stdout outputs are also byte-identical. All checks raise
explicit exceptions; Python optimization does not remove them.

Reproduce from this directory:

    python3 verify_affine_circulants.py --max-m 12 --output verification.json > verification.out
    python3 -O verify_affine_circulants.py --max-m 12 --output verification.optimized.json > verification.optimized.out
    cmp verification.json verification.optimized.json
    python3 independent_affine_audit.py > independent_audit.out
    python3 -O independent_affine_audit.py > independent_audit.optimized.out
    cmp independent_audit.out independent_audit.optimized.out

Main script SHA-256: `2133b46abb7c4cd4c4f3919f5d0decd8cad6d097fe4032c21c907361ca6276c0`.
Independent script SHA-256: `e6c9114daad08496053d2153afca25507dcc8ac0759eb9eac8cff71340ff0634`.

## Precise next targets

1. Determine the least Boolean degree of an equivariant type-bijection for the
   first obstructed primes p=13 and p=17. This packet proves degree at least two,
   but does not claim degree two suffices uniformly.
2. Classify equality and next-best coverage for restricted-support h, retaining
   the stable-factor gcd and the primary-height coordinate separately.
3. Separate Boolean degree from generator-free naturality and first-order
   definability: these are distinct constraints, not interchangeable no-go claims.

The relevant inherited meta-pattern is controlled forgetting with a sidecar:
complementation is safely forgotten by Δ, while primary height must remain
available to a nonlinear routing decision. No new general research-process
principle is claimed from this one packet.
