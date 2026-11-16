# A positive readout for literal span 3, and the exact missing cap

**Update:** the physical comparison is no longer missing. [Sharp one-vertex growth](SHARP_GROWTH.md) proves the factor 39 on all span-3 prefixes and gives an explicit cap on the delayed sum. The discussion below records the original reduction.

Use the [THM-4880 frontier transfer](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/ribbon_tournaments_bounded_span_transfer_20261104.md): 119 reachable fragment profiles, nonnegative integer maps M_a for the eight letters, initial empty profile, and a row f that reads the 14 one-fragment terminal profiles. The coordinates of v(w) count actual partial configurations, and H(w)=f v(w) on every nonempty word.

## Two unconditional upper growth certificates on the frontier

Use integer rows throughout. Starting from r=1 on all 119 profiles, repeatedly replace each coordinate j with

    max(r_j, max_a ceil(10(r M_a)_j / 39)).

The exact recipe stops after 64 iterations, has smallest weight 5 and largest 78, and satisfies all 952 inequalities `10 r M_a≤39 r` entrywise. The accompanying verifier generates the raw matrices and checks those inequalities directly. Since r≥f, this proves both the uniform one-letter bound for **this positive frontier mass** and the global upper bound H(w)≤r_empty·(39/10)^|w|. This is not the desired factor 39 bound on the ratio H(wa)/H(w); raw r may count large unfinished states that f ignores. It is also not a lower-growth certificate.

A canonical alternative is to define

    g = Σ_{|v|≤3} f M_v,        G(w)=g v(w)=Σ_{|u|≤3} H(wu).

The counts of raw profiles on which the cumulative row is nonzero are 14, 67, 111, 119 for depths 0, 1, 2, 3. Thus three-letter smoothing makes the row strictly positive without changing the underlying tournament family. Its entries range from 208 to 5512, and the verifier checks the sharp maximum *over raw states and letters*, `g M_a≤(1809/214)g`.

## Exact relation to the required bounded-growth lemma

Uniform bounded literal growth, H(wa)≤C H(w) for all words and letters, is equivalent to a single comparison `G(w)≤K H(w)` on the **physical reachable orbit**.

- If the one-step bound holds, the comparison holds with K=1+8C+64C²+512C³.
- Conversely H(wa) is one of the summands of G(w), so the comparison implies the one-step bound with constant K.

Any positive raw row (including r) is comparable above and below to g by constants on the nonnegative orthant. Thus bounding `r v(w)/H(w)` along the actual orbit is also equivalent. The explicit positive envelopes remove the raw zero-row problem if used as observables, but do not show that the output is comparable to them. Certifying an upper comparison on the full raw orthant would be false, since f is zero on multi-fragment profiles.

If an orbit cap is proved and a single-letter potential grows at rate ρ on every rank-preserving segment while charging a fixed loss to each rank drop, the actual H will inherit its lower growth. The cap alone supplies neither the rank-preserving transition relation nor the lower-growth inequalities. This is the distinction between the verified restricted thirteen-unit lower-growth certificate and a new all-literal one.

## Exact sharpness in the transitive-prefix subfamily

When the old ribbon is transitive, put its vertices in actual topological order. Any label more than three below the oldest potential marked out-neighbour of the appended vertex must precede all of its out-neighbours. Every directed cycle uses the appended vertex and, afterward, an increasing old path starting at an out-neighbour. Thus such early vertices cannot lie in any directed cycle. By the odd-cycle packing formula they have no effect on H. With old labels 0,…,n−1 the possible out-neighbours have labels n−3,n−2,n−1, so only the last six old labels matter. Restriction preserves span ≤3 and leaves all possible marks intact.

The verifier enumerates all transitive orders of at most six old labels compatible with span 3 and the eight append masks. The maximum is exactly 39. One six-vertex old topological order is (3,0,4,1,2,5); the appended vertex beats labels 3 and 4. The full outgoing bitmasks are (118,100,96,55,38,64,24), with H(old)=1 and H(new)=39. This proves the proposed constant in the transitive case and makes 39 a rigorous lower bound on any all-prefix constant, not an upper bound on arbitrary ribbons.
