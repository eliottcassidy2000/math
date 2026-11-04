# A strong 24-vertex realization of 4,611,841 and its complete physical ports

2026-11-04. **CONSTRUCTIVE and FINITE-EXACT:** four distinct labelled strong tournaments with 4,611,841 Hamilton paths (six prefix/suffix representations), all physical endpoint fields and directed joints, and a reusable exact compiler for tournaments which become transitive after deleting one named root. **Inherited restriction:** a clock-uniform transitive-suffix grammar has at most one ternary digit at each physical port signature. No native guarded program, continuation of the frozen 17/18-vertex lineage, or unrestricted two-branch closure is claimed.

## 1. Literal construction

Use physical vertices 0,1,...,22,r, with r encoded as -1 in JSON. Every edge between 0,...,22 points from smaller to larger. The vertices beating r are exactly

    {2,5,6,10,12,13,15,17,20,22}.

The root beats every other vertex, including 0. These rules specify every edge. `candidate.json` contains the literal 24 rows in the declared name order, along with independently compiled endpoint data. This is a new representative, not a renaming of the frozen S182 source.

It is strong: r reaches 0, 0 reaches every other old vertex, every old vertex reaches 22, and 22 reaches r. More generally an almost-transitive tournament on ordered old vertices is strong exactly when the first old vertex loses to the root and the last beats it (for at least two old vertices). Necessity follows from a source or sink; the displayed cycle routes prove sufficiency.

Regard 0 as the transitive prefix TT1 and 1,...,22 as suffix TT22. The actual suffix port bitmask is

    f=2,710,066, W={1,4,5,9,11,12,14,16,19,21}

in zero-based suffix positions. Applying the inherited S182 exact multiport law gives

    j=808,292,
    H=1+2f-j=4,611,841.

The defect is not zero. It is retained and recomputed below.

## 2. Independent large-count audit

Delete r. Every path segment on the remaining transitive vertices is forced to follow increasing order. A full Hamilton path therefore corresponds bijectively to a subset B placed before r, with its complement after r. The subset is admissible precisely when

- B is empty, or its maximum beats r;
- its complement is empty, or its minimum loses to r.

`audit.cpp` enumerates all 2^23 subsets separately for each of six representatives, using only these two actual edge tests. It obtains H=4,611,841 in all six cases. This is independent of the binary-defect formula used for discovery and of the finite-state endpoint compiler used for the port data.

## 3. Six representatives and completeness within this family

The six pairs (prefix length t, suffix bitmask f) are

    (1,2710066), (1,3035358), (1,3186353),
    (1,3568813), (2,1355033), (2,1517679).

All have total order 24. Their root defects are respectively

    808292, 1458876, 1760866, 2525786, 808292, 1458876.

For a transitive suffix with terminal vertex a port, the inherited exact scan starts g=j=e=0 and reads f from low to high. At bit b it makes

    j <- j+b(g+e), g <- 2g+b, e <- b.

Equivalently, if port positions are i<k, their contribution to j is 2^(k-i-1)+[k=i+1]. This retains every physical pair and is unchanged by prepending nonport positions.

The search in `search.cpp` is exhaustive for this narrowly specified family: t>=1, nonempty transitive suffix, arbitrary actual port subset meeting its terminal vertex. Its count is 1+2^t f-j. Since j<=f (every joint cover is among the f endpoint covers), H>=1+(2^t-1)f. Thus f<=(H-1)/(2^t-1), and t<=22 for H=4,611,841. The code conservatively includes t=23, where the bound gives no candidate. Every positive f canonically fixes the suffix length by its highest set bit; appending nonport suffix vertices would fail strongness. Exactly the six displayed pairs pass the equality test. There is no smaller-order representative **in this family**, not a lower bound for general strong tournaments. The last two pairs repeat the first two labelled adjacencies because they have the same root cut f*2^t. There are four distinct labelled adjacencies, not six nonisomorphic graphs; no isomorphism classification is claimed.

## 4. Exact reusable endpoint compiler

Let the old vertices form a transitive order, with arbitrary actual root incidences. Count ordered k-path covers with empty slots allowed, where k=1 or 2 here. Specify which slot contains r. Split that slot into the increasing segment before r and the increasing segment after r. Every other slot is one increasing segment. Assign each old vertex to one of these k+1 segments.

The only edge restrictions not forced by increasing order are the last vertex of the before segment beating r and the first vertex of the after segment losing to r. During an increasing scan it suffices to retain:

- before state: empty, last vertex beats r, or last vertex loses to r;
- after state: empty or nonempty, rejecting a first vertex which beats r.

At completion discard the bad before state. Adding a vertex to an ordinary segment changes neither state. Sum over the root's slot. This is a bijective compiler with at most six states for each endpoint constraint, not an approximation or a modular compression.

To require a specified path slot to end at an old physical vertex w, force w into its terminal segment and disallow later vertices from that segment. A start constraint forces w into its initial segment and disallows earlier vertices there. If the endpoint is r, require the corresponding before or after segment to be empty. If an ordinary slot is required to have endpoint r, reject it. Multiple constraints are imposed simultaneously. In particular this computes every joint J_wv by requiring slot one to end at w and slot two to start at v.

Every slot's physical path is uniquely reconstructed from its increasing assigned segments and r. Conversely every actual cover supplies exactly one assignment. This proves the compiler for arbitrary order, cut, physical constraints, and repeated reconstruction after births or re-marking. It does not say that only scalar H,f,j, or observer residues suffice for future updates: the full order and root cut are retained.

For each of six representatives, `verify.py` computes H, E1, S1, deletion H, E2, S2 at all 24 physical vertices and every one of 552 directed joints. At each w it independently checks the bad-join identity

    E2_w-H = sum_{v -> w} J_wv

using actual rows, so all local defects are transported. The sum is over the actual incoming neighbors, not a field-phase alias.

## 5. Current field observers and hostile controls

The inherited codec is copied into `codec.py`, with provenance in its header; it is not a new field construction. At every physical port form

    golden=(H,(S2-E1)/2,(S2+E1)/2).

Integrality, inverse linear coordinates, and F343 Frobenius closure are checked. F343 uses X^3+X^2+1 over F7. At the first representative's root:

    E1=S1=1, deletion H=1,
    E2=5,420,133, S2=7,137,627, defect=808,292,
    golden=(4,611,841,3,568,813,3,568,814),
    F49/charge=(6,4;6), Frobenius=(1,3;6),
    F343=(3,3,4) -> (4,6,4) -> (3,5,6) -> (3,3,4).

These are recalculated from this actual graph. They are not inherited from the S182 39/50 ports. In general, with n ordered old vertices and cut mask c of vertices beating r, E2_r=1+c and S2_r=1+rev_n((2^n-1) xor c). This follows directly by choosing the last before vertex or first after vertex in the two-slot cover; it is also checked by the compiler. The field views do not replace the cut.

`controls.py` recomputes the inherited hostile pair with independent subset endpoint DP: both H=462915, but actual birth cut 256 gives 1061775 and 1395279. The inherited adjacency-preserving automorphism enumerator again returns full orders 81 and 27. A triangle search after every vertex deletion proves neither hostile control is almost transitive. Thus this special compiler does not purport to compress those controls or identify their graph automorphisms with F49/F343 Frobenius.

## 6. Inherited branching restriction and the stronger companion result

This section restates the one-digit restriction already proved in the [preceding transitive-chart result](https://github.com/eliottcassidy2000/math/blob/3f2fa002801868c280fafed0c05a1d9060df6c86/05-knowledge/results/multiport_defect_grammar_20261103/TRANSITIVE_CHART.md); no novelty is claimed for it. Within this transitive-suffix family, f is the actual binary encoding of the suffix ports; it uniquely fixes the canonical suffix and j(f). A clock-uniform ternary branch with unchanged external prefix clock must have

    f'=3f, 1-j(f')=3(1-j(f))+d.

Consequently its only possible digit is

    d(f)=3j(f)-j(3f)-2.

There cannot be both a +2 and a +4 output from the same canonical transitive-suffix state while requiring the identities at every common clock. Allowing arbitrarily many suffix ports does not remove this obstruction. A useful grammar must enlarge the suffix family, change the interface/clock requirements, or relax clock-uniformity. This does not obstruct a finite count-specific repair or general ordered-component suffixes. The inherited f=7 example has d=4; it is included as a positive control in the search output. The companion [all-mask classification](../other/PROOF.md) goes further: no +2 exists anywhere in this chart, and the only +4 maps are padded 7 to padded 21, allowing compensating fixed clock shifts. Those are results of the companion packet, not inferred from this finite search.

## 7. Verification and scope

`verify.py` compares the compiler to independent full subset endpoint DP for all 127 almost-transitive root cuts on orders 1 through 7: 3,076 endpoint comparisons and 4,092 directed-joint comparisons. The six large representatives additionally check strongness, exact H, the formula, all physical defects, and field codecs. Normal Python and optimized Python have identical stdout and receipts. All guards use explicit exceptions.

Run:

    g++ -O3 -std=c++17 search.cpp -o search
    g++ -O3 -std=c++17 audit.cpp -o audit
    ./search
    ./audit
    python verify.py
    python -O verify.py
    python controls.py

No installs, startup helper, or forbidden mistakes file were used. Source ideas are the S182 `induction/PROOF.md` binary-port formula, S182 `defect/PROOF.md` physical joint/defect distinction, and S182 `phases` codec and hostile controls. This packet adds the explicit target realization, independent large-count subset audit, full physical endpoint compiler, current field data, and an explicit reminder of the inherited family-specific branching boundary, with the sharper companion classification linked above. It does not furnish a legal native program from the frozen named source or a uniform two-digit construction.
