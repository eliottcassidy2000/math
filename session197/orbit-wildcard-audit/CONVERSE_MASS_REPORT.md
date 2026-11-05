# Independent weighted converse-mass audit

**PASS.** The universal weighted mass theorem, its positive odd fibre weights, its complete independent-bit construction, and its ordinary self-converse class parity consequence are correct. The proof handles the one-vertex case explicitly and the optional empty-tournament convention separately. No correction is required.

## Frozen source snapshot

The source files are in `session197-orbit-wildcard/`:

- `CONVERSE_MASS.md`: `07f4b20dde9f9d23b7711f65bd484d7468bdbef3cc6da0cc2888aa91d402a620`.
- `verify_converse_mass.py`: `18eb958c031f8313db66adb3b842743b7fd479f6f7fbc153c35289fcb6b32a28`.
- `converse-mass-normal.json`: `7731d4f313ae6029a0170f371fba2c2b931e7e2b885abac1c822f8b533edf4c9`.

For a self-converse tournament `T`, an anti-involution `tau`, and `G=Aut(T)`, the exact claim is

    c(T) = Fix(J_tau on Hamilton words) / |C_G(tau)|,
    sum_[T], self-converse of order n c(T) = 2^floor((n-1)^2/4).

This is a weighted enumeration identity, not an unweighted class count or a spectrum coverage theorem.

## Proof checks

1. **Existence and conjugacy.** For `n>=2`, an arc-reversing permutation has even order. Its square belongs to the odd-order automorphism group, so its order is `2m` with `m` odd and its `m`th power is an anti-involution. For two anti-involutions, `g=tau sigma` has odd order `m`, and direct multiplication gives `g^((m-1)/2) tau g^(-(m-1)/2)=sigma`. Thus the conjugating element is actually in `G`. The `n=1` identity convention is valid without an order-two assertion.

2. **Integrality and parity.** `J_tau` permutes Hamilton words as an involution. The odd Hamilton count therefore makes its fixed set nonempty and odd. The centralizer preserves this set and acts freely, since fixing a full vertex word fixes every vertex. Its order is odd and divides the fixed-set size, so the quotient is a positive odd integer. Conjugacy makes it independent of `tau`.

3. **Normalized bits.** Fix the Hamilton word `0,1,...,n-1`. Consecutive arcs are forced, and reflection of a forward-orientation bit is equality, not complementation. Among the `binom(n-1,2)` nonconsecutive pairs, exactly `floor((n-1)/2)` are reflection-fixed. Thus the number of free pair-orbits is `floor((n-1)^2/4)`. The sole fixed consecutive pair for even `n` is correctly excluded. At `n=1` and `n=2` there are no free bits. At `n=0`, the separately defined empty case has one object and exponent zero.

4. **Exact fibres.** Two `tau`-symmetric Hamilton words give the same normalized labelled tournament precisely when an automorphism carries one word to the other. Symmetry of both words forces that automorphism to commute with `tau`, by equality on every vertex. Conversely, an arbitrary normalized object in the class pulls reflection back to an anti-involution; conjugacy permits adjustment of its isomorphism to the fixed `tau`. Hence the normalized presentations of the class are exactly the centralizer-orbits of symmetric words, not all automorphism-orbits and not an additional quotient by reflection.

5. **Equivalent ordinary path-orbit count.** An invariant free `G`-orbit contains `|G|` pairs consisting of a word and its unique word-reversing anti-involution. The uniqueness follows by prescribing its value at every word entry. Summing over the single anti-involution conjugacy class of size `|G|/|C_G(tau)|` gives precisely `f_tau/|C_G(tau)|` fixed ordinary path-orbits. This agrees with the direct normalization proof.

The inherited facts used here are oddness of tournament Hamilton counts and free automorphism action on full vertex words. No classification of tournament automorphism groups or unproved numerical closure is needed.

## Independent exact controls

`check_converse_mass.py` imports no author verifier. It rebuilds all tournament isomorphism classes through order seven by every one-vertex extension and exact outdegree-preserving canonicalization. The complete class counts at orders zero through seven are

    1, 1, 1, 2, 4, 12, 56, 456.

It decides self-converseness separately by canonicalizing the converse, then compares with a separately generated normalized universe. This gives self-converse class counts

    1, 1, 1, 2, 2, 8, 12, 88.

For **all 115 self-converse classes** in this range, the checker explicitly enumerates Hamilton words, every anti-involution, the exact automorphism group, and the centralizers. It verifies the displayed conjugating power, the exact set of normalized presentations, each individual centralizer fibre, and the fixed count on ordinary Hamilton-word orbits.

In particular the controls include nontrivial centralizers, so the division by `|C_G(tau)|` is exercised rather than vacuous. At order six, the graph with rows `[10,20,33,38,13,19]` has automorphism group and centralizer of order three, nine fixed words, and fibre weight three. The ordered pair of directed triangles has automorphism order nine, centralizer order three, three fixed words, and fibre weight one.

At order seven the fibre histogram is

    weight 1: 10 classes; weight 3: 12; weight 5: 16;
    weight 7: 32; weight 9: 18.

Its weighted sum is 512, the required `2^9`. All smaller histograms match the source receipt. Normal and optimized runs of the independent checker agree byte for byte. The source verifier was also replayed normally and under `-O`; both outputs match its frozen receipt exactly.

## Reproduction and scope

Run `python check_converse_mass.py` and `python -O check_converse_mass.py` from this audit directory. The only local helper dependency is the independent `check_graph.py` for explicit exceptions and exact automorphism backtracking.

- Independent checker SHA-256: `91b979f85f1c22df2c7cea8efa32a7f26256322da534261559752bf9791fbed5`.
- Independent receipt SHA-256: `43fb5b2e762594cf68ed5bcbc7c00c3ade09e4fe1d726a9c97d5f1fce1d34bd4`.

The eleven- and fourteen-vertex example contributions 111 and 6111 are supported by the separate `SELF_CONVERSE_REPORT.md` audit. All source files remained unchanged. The theorem supplies a complete finite search universe at each order; it does not establish that every odd prime or `7p` target occurs.
