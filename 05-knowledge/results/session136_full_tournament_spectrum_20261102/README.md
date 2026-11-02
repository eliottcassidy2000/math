# Session136: structure of the full tournament spectrum

**The conjecture that H takes every positive odd value except 7 and 21 remains OPEN.** Here H(T) counts directed Hamiltonian vertex words, not path orbits. This session establishes five complementary frameworks, with independent proof review and exact executable controls. It makes no literature-priority claim.

## 1. Every local congruence is realizable

Strong one-apex tournaments realize every parity-compatible residue modulo every modulus, simultaneously for compatible CRT requirements. Every such residue contains arbitrarily large counts. Their profinite closure, and that of the strong and unrestricted spectra, is `1+2 Zhat`.

The proof uses actual positive binary append words, a modular commutator, and a dyadic reset preserving the odd part. Its explicit word bound is polynomial in the numeric modulus, not its bit length. Yet this same family misses 13 exactly. No finite congruence filter on H alone can explain an odd missing value; magnitude or additional structure must be retained.

[Proof](local_closure/PROOF.md) · [Constructor](local_closure/check.py) · [Independent audit](independent/LOCAL_CLOSURE.md).

## 2. Exact context equivalence and substitution width

If a tournament is inserted as a homogeneous block, contexts with m exterior vertices see precisely its first `floor(m/2)+1` path-cover coordinates. The first differing coordinate c requires exactly `2c-2` exterior vertices to distinguish.

An interacting module of order n forces `H >= 2^n+1`. This gives a finite exact algebra for capped H and an all-depth fixed-point algorithm for any finite substitution grammar. Define the **spectral substitution width** of h as the smallest maximum quotient order needed to construct some tournament with H=h. The certified thresholds are

    w(13)=w(19)=5,       w(31)=6.

The complete width-at-most-five spectrum through 127 contains 49 odd values. Its missing odds are

    7,21,31,35,43,53,67,71,77,83,89,101,103,113,127.

These are grammar-relative obstructions, not new globally forbidden values. Adding a graph as a seed is weaker than allowing it as a substitution operation. Whether a bounded quotient width suffices for the full spectrum remains OPEN.

[Proof and finite algebras](gluing/PROOF.md) · [Exact width-five certificate](gluing/width5.normal.json).

## 3. Insertion responses are positive network cuts

Adding one vertex with prescribed out-neighborhood S gives exactly `F(S)=H(T)+2 R(S,V\S)`, where R is a canonical nonnegative integer network built from even-cardinality paths and complementary H counts. Minimum response is a directed mincut once this network is supplied.

Strict uncrossing bounds the bottom response band. Vertex-transitive symmetry restricts its cuts to blocks of imprimitivity; primitive symmetry leaves only singleton/complement cuts. For a circulant parent of order n, the exact minimum needs only `tau(n)-1` subgroup-cut representatives, rather than all `2^n-2` nonconstant cuts. Computing each queried H is still a counting problem.

The cyclic order-eleven parent has exactly 1,062 missing odd levels in its bottom response band. Other parents can fill those values: this is a one-parent obstruction, not a hole in the full spectrum. A six-vertex hostile shows that singleton cuts need not minimize without the symmetry hypotheses.

[Positive-cut proof, controls and symmetry refinements](ear_geometry/positive_cut_geometry.md).

## 4. Forbidden-basis lifetimes and catalogue growth

A core with `h=H(T)` and `d=max_v H(T-v)` is a minimal forbidden tournament for `H<=B` exactly when `d<=B<h`. The resulting class-labelled intervals yield the conservation law

    b(m)-b(m-1)=beta(m)-a(m),

where a counts exact-H core fibres and beta counts births. Keeping births is essential; equal basis cardinalities alone cannot certify a missing interval.

Complete explicit forbidden bases necessarily have superpolynomial spikes. This is a lower bound on the number of listed tournaments, not on compressed certificates or decision algorithms. A constructive companion gives `b(B)>=2^floor(log_15(3B))` for every B>=5. The exact census includes 104 distinct strong order-seven minimal obstructions at threshold 37.

[Proof, descent graph and hostile examples](obstructions/THEORY.md) · [Exact certificates](obstructions/normal.json).

## 5. Smaller witnesses through finite-bank repacking

Every **attained** h>1 has a witness on at most

    ceil(22 log_{80267}(h)) = ceil(1.9480898075965785... ln h)

vertices. The proof repacks SCC factors using inherited finite strong-witness banks. It improves the general induced-core bound while allowing an unrelated, possibly nonstrong witness. It supplies no witness for an unproved target and does not regenerate the inherited large bank.

The minimum-order function has an exact min-plus/divisor-DAG description. Its cost-saving atoms differ from arithmetic multiplicative atoms. Separately, a two-value multiplicative mesh shows why vanishing relative gaps still do not imply full integer coverage.

[Compression and minimum-order theorem](spectrum_structure/THEORY.md) · [Relative-gap boundary](spectrum_structure/COFINALITY.md).

## Synthesis, verification and inheritance

[What changes and what remains](SYNTHESIS.md) · [Verification](VERIFICATION.md) · [Independent aggregate audit](independent/REPORT.md) · [Concept board](BOARD.md).

The inherited GF(2)-similarity lens is explicitly tested: two order-six adjacency matrices with the same single nilpotent Jordan block have H=1 and H=41. [Literal certificate](LENS_BOUNDARY.normal.json). The [interface checks](INTERFACES.normal.json) also show actual counts 2007,2013,1021 agreeing modulo 1000 with 7,13,21, and the exact contextual effect of transitive padding.

Chan–Pak record the exact-spectrum conjecture in [Section 6.5 of Computational complexity of counting coincidences](https://arxiv.org/html/2308.10214v2#S6.SS5), crediting bof and Royle's 2016 question. Rédei parity, the all-order exclusions 7/21, the OCF, the strong floor, prior substitution kernel, and prior finite witness banks are explicitly inherited. None is presented as newly discovered here.

[Recent extremal source check](literature/EXTREMAL_PATH_SOURCES.md): the cited sharp asymptotic maximum constrains witness capacity, not interior spectrum membership.
