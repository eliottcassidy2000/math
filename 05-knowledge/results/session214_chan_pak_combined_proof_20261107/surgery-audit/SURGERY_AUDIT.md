# Session214 two-apex surgery and all-order omissions audit

**Verdict: PASS.** The restricted presentation class F has exactly the four missing positive odd counts 7,21,7477,27371 through 27371. The two new exclusions are all-order inside F. Explicit two-exit presentations attain both 7477 and 27371. No ambient exclusion or completeness of the enlarged class follows.

## Frozen author proof and computation

Staging packet: `session214-surgery`.

Initial certificate:

- PROOF.md: `1675fa2d54dc89aecf4383941215f47dfbb620a4a971a5fa3edbffddf40c89ea`
- counting.py: `825d5beb1f350208a47e0b76e1b299fdaf4b6d9e11d428776bc6ffc8532d0a32`
- verify.py: `4c6eda625d8bbe2428555c0335a97bd3c91ad7731c169bf8d7fa466c8000f4b2`
- RESULT.json: `65244bc13d57a2dcae9c14757aefa84b1d31fd06e940f97e7a8d8658f1d1b232`
- SOURCES.json: `60c982245295afa4fe7d14c6906677f45b28ef2ce3526d7ca3017ba093450831`

Composite-target addendum:

- ADDENDUM.md: `b86c0062d1d73a046181a8aa6362af4293006a4f782dc57b3ba8ac1540cece96`
- addendum_verify.py: `ea04a9dd4b99902f2353782cba8465a11c4ce2ca094f32b3dafb9fc5c0d68ebd`
- ADDENDUM_RESULT.json: `39a2a27ace8c3e73852cdb2bfc6036f9d0956e1eaf1414a2f83932e0f517d093`

The initial frozen manifest was rehashed. Author normal/optimized outputs agree. The all-order strong-floor input is the explicitly cited Busch theorem, not an inference from a finite census.

## Independent proof audit

F is hereditary when either apex disappears: the remaining unrestricted apex can retain that role. Strong-component products are forced by the acyclic dominance order. Primality of 7477 therefore reduces its occurrence to a single strong component in F. Busch gives f(18)=9375, so that component has order at most17. The author's count-preserving sink padding covers zero exits, one apex and the small boundaries. It does not assume bounded order for arbitrary padded nonstrong tournaments.

The addendum correctly distinguishes a single strong two-apex core from two nontrivial SCCs. Every nontrivial SCC must contain an apex; with two such components the restricted apex must be later, since otherwise it beats at least two chain vertices in the other component. Its sole exit is forced to its own first chain vertex. Hence it is R_m, of count 2^(m-2)+1, while the earlier component is one-apex. The displayed converse product presentation is valid.

Trial division independently confirms 7477,101,271 prime and 27371=101*271. No divisor of 27371 has the form 2^k+1 for k>=1. Thus a hypothetical 27371 presentation has only one nontrivial SCC. The cited f(20)=28125 bounds its order by19. A zero-exit restricted apex makes the remaining chain-plus-apex transitive after deleting the unrestricted apex, so that boundary is one-apex rather than an omitted two-apex case.

The one-apex bound H>=1+2^(L-2) is valid for L>=3: strongness forces the two extreme chain incidences; every even subset of interior vertices gives a distinct odd directed cycle, and the inherited positive odd-cycle expansion contributes at least two for each. It bounds one-apex cores at this target by chain length 16. The cases L<3 are already in the finite census. Hence the all-order composite reduction is exhaustive.

The new scope paragraph is important: arbitrary dominance products of F still omit prime 7477, but recover 27371 as 101*271. The latter is an F omission only, not an omission from its multiplicative closure. Both positive controls have the stated adjacency rows and exactly two restricted-apex chain exits.

## Independent counting engine

The author's three-run transfer was not imported. `surgery_cover.cpp` and `surgery_cover_extended.cpp` instead delete the second apex b and compute all subset start/end path counts of the remaining tournament X. For i!=j, let W_ij count two disjoint nonempty spanning paths of X, the first ending at i and the second starting at j. If z_i indicates b->i, direct splitting at b gives

    H(X+b)=sum_i E_i(1-z_i)+sum_i S_i z_i
             +sum_(i!=j) W_ij(1-z_i)z_j.

Thus the independently derived Boolean quadratic has constant H(X), linear coefficient S_j-E_j+sum_i W_ij, and pair coefficient -(W_ij+W_ji). A Gray traversal evaluates every cut. This is a distinct cover-convolution route to the same exact finite counts, rather than a replay of the run-state interpolation.

Initial independent cutoff:

- 1,835,012 exact-one-exit two-apex presentations with chain length<=15.
- 131,071 one-apex cuts with chain length<=16, including the zero-exit boundary.
- No 7477; all smaller positive odd values except 7,21 attained.

Extended independent cutoff:

- 8,388,612 exact-one-exit presentations with chain length<=17.
- The same 131,071 one-apex boundary cuts.
- No 7477 or27371; all smaller positive odd values except 7,21,7477 attained, hence 13,682 positive values through 27371.
- 1,795 separate direct-DP controls of the polynomial counts.
- Both supplied two-exit tournaments directly counted, giving 7477 and 27371 respectively.

Reproduce either engine with:

    g++ -O2 -std=c++17 -Wall -Wextra -Werror -pedantic surgery_cover.cpp -o /tmp/session214-surgery-audit
    /tmp/session214-surgery-audit > /tmp/session214-surgery-audit.json
    cmp surgery_cover.json /tmp/session214-surgery-audit.json

Replace surgery_cover by surgery_cover_extended for the addendum. All checks use runtime exceptions. The numerical universes are labeled parameter presentations, not isomorphism-class counts. The proofs provide the all-order cutoffs; the exact computations decide those finite universes.

## Final synthesis consequence

Swapping apex labels and taking tournament converse preserve H and a transitive chain. Therefore every two-apex presentation of either new target must give each apex at least two incoming and two outgoing chain arcs. This elementary consequence is valid and does not classify all target realizers. Chan–Pak completeness remains OPEN.
