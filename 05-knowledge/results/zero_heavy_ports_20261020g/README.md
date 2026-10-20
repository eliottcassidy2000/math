# Zero-heavy cores, exact ports, and product symmetry

**2026-10-20. Current scope:** unrestricted A000568 compression and all-order
black-core connectivity remain **OPEN**. This packet separates exact
operation-preserving encodings from counts that forget their ports.

| Lane | Result | Status and source |
|---|---|---|
|Zero-heavy boundary|One-feedback-vertex class count A_1(n)=2^(n-1)-(n-1)-binom(n-2,2) for n>=2; exactly two pair-free classes at each n>=4|PROVED, [anchor quotient](zero_heavy.md)|
|Fixed anchor number|Canonical re-anchoring by bounded-delay port-word transduction; rational class-series proposal|Proof candidate under independent audit, [construction](zero_heavy.md#5-candidate-all-fixed-a-canonical-regular-language)|
|Separating anchors|The ordered join of m directed triangles needs exactly floor(3m/2) anchors to distinguish all remaining vertex signatures|PROVED, [sharp hostile](zero_heavy.md#4-why-a-uniformly-small-separating-anchor-cannot-suffice)|
|One-vertex numeric growth|Every directed path-cover count has an exact quadratic attachment response, with specified component endpoints; its labelled average depends only on the old profile|PROVED, [port tensors](profile_growth.md)|
|False smaller growth state|Pair-free six-cores with equal full profile, scores, feedback and automorphism cycle index have maximum child H values159 and147|FINITE-EXACT counterexample, [witness](profile_growth.md#2-smallest-reinforced-profile-hostile-in-the-audited-universe)|
|Square/diagonal microcosm|A product of pointed Boolean sets has literal commuting involutions whose orbit counts give the primitive Pythagorean triple; actual deletion fibres have sizes4,2,1|PROVED, [symmetry and growth](microcosm.md)|
|Capped weighted core|Same-order eleven-vertex family members have identical cap and scores but blue counts17 and0; actual chain-length equality is lost|PROVED example, [predicate boundary](microcosm.md#6-the-same-five-vertex-core-reaches-the-zero-heavy-boundary)|

The board comprises feedback anchors, canonical word ownership, homogeneous
chains, two-component endpoint tensors, blue positional involutions, and
the exceptional state in the square product. The anchor is the unresolved
zero-heavy/strong-blue boundary; the niche is exact one-vertex response;
the wildcard is the commuting-involution interpretation of (3,4,5).

The closest inherited mechanisms and hostile cases are named at the start
of each note. The inherited run-word clone dichotomy is checked again but
not claimed as new. THM-1900's incorrect deletion-of-path premise is repaired
in place, retaining the true neutrality theorem and its orphan terms;
THM-012b's three-cycle-only range is corrected from n<=5 to n<=4.

The exact checks are ordinary Python with explicit failures. Their declared
universes, independent paths, and limits are recorded beside each proof.
No finite census is used as an all-order connectivity argument.
