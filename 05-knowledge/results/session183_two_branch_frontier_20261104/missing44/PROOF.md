# The missing 44 child on the actual defect-carrying lineage

2026-11-04. **FINITE-EXACT physical realization; PROVED arbitrary-prefix isolation.** Native-operation certification is a separate independent audit. No reusable two-digit closure or global spectrum theorem is claimed.

## 1. Exact physical continuation

The source is Session182 `closure/candidate4.json`, itself the actual Session181 lineage. Its eighteen physical names, root39, prefix50→51→55, fourteen suffix names, and every directed edge are frozen in `source.json`. Its H is1537279, suffix h58063, F190302, J43200, and c=h−J=14863. It is not replaced by a numerically equal representative.

Keep every name, keep the entire TT3 prefix and every prefix-to-suffix incidence, and make precisely these three suffix reversals:

- {1,4};
- {14,21};
- {33,54}.

Change root39's actual in-port set from{33,40} to{1,4,14,21,32,40,54}. This means exactly seven root-edge reversals at names1,4,14,21,32,33,54. No permanent vertex is born or deleted. `candidate44.json` contains the complete named endpoint rows and all ten changed physical pairs.

Independent suffix endpoint-table recomputation gives

    h=40085,
    F=80067+105391+227311+72911+41747+40675+65795=633897,
    J=499420,
    c=h−J=−459335.

Thus the full rooted tournament satisfies

    H=c+8F=4611841=3·1537279+4,
    E_39(2)=h+8F=5111261,
    E_39(2)−H=499420.

Its full adjacency is strongly connected, checked by direct reachability from every physical name. The root defect has increased from43200 to499420; it is not zeroed or copied from an ancestor.

`verification.json` retains all196 directed suffix coefficients J_uv, every one of the fourteen physical one-port responses F_u, every inner defect F_u−h, and all49 selected-port directed terms. It checks the bad-join identity F_u−h=sum_(v→u)J_uv separately for each physical suffix vertex. Nothing identifies automorphisms or observer phases with these directed responses. Any inherited field-phase claims require independent recomputation from these actual rows.

## 2. Exact search scope

`search.cpp` enumerates unordered triples of distinct internal suffix edges in lexicographic index order. At each triple it considers every nonempty physical port set meeting the terminal suffix SCC and every transitive prefix length t≥3 until the target is exceeded. This preserves the old TT3 names and allows transitive prefix extension. It stops at the first exact missing44 hit, indexed(a,b,c)=(1,63,89), W=13509, t=3. This was the7,569th tested triple. The complete three-flip menu was **not** exhausted; no new negative census is claimed.

The loop is finite and its clock cutoff is sound: J_W≤F_W and h≥1 imply H_t≥h+(2^t−1)F_W. In particular a target below2^23 cannot require t≥23. The implementation's limit32 is redundant. The pruning F_W≤target−h is safe already at t≥1. All endpoint tables and products used at suffix order14 fit signed64-bit integers. The seven actual selected ports, rather than an aggregate scalar substitution, realize the hit.

Session182 excluded missing44 for at most two suffix reversals, all valid physical ports, all transitive prefix extensions, and arbitrary one- or two-birth prefix extensions. This new positive witness crosses the excluded suffix-reversal boundary by exactly one additional suffix reversal. It uses no prefix birth at all.

## 3. This particular rewrite is isolated to TT3

Keep the source and output suffixes and port sets above, but replace the unchanged prefixTT3 by any nonempty tournament A. Write a=H(A), p=P_A(2), with both ordered path slots allowed empty. The inherited exact evaluator gives

    H_out−3H_in=−503924a+62991p.

Hence this rewrite has digit4 precisely when

    62991(p−8a)=4(1−a).                       (1)

In fact (1) forces A to be transitive of order3. First, a≥1 and (1) imply p≤8a. If A has n vertices, cutting each Hamilton path at each of its n+1 positions injects (Hamilton path,cut) into ordered two-path covers, so p≥(n+1)a. Thus n≤7. Also gcd(62991,4)=1, so (1) forces a≡1 mod62991. But a≤n!≤5040; therefore a=1 and p=8.

A tournament with exactly one Hamilton path is transitive: a directed triangle supplies three distinct Hamilton paths, each extendible by deterministic insertion of the remaining vertices while preserving their relative order. A transitive tournament of order n has p=2^n, so p=8 forces n=3. Conversely TT3 has(a,p)=(1,8), and the displayed construction verifies digit4.

Therefore this particular simultaneous suffix/root rewrite works at exactlyTT3 up to physical names, even allowing arbitrary prefix interiors. It is a finite return, not a repeatable prefix-independent digit map. For TT_t the offset is−503924+62991·2^t, which equals4 exactly at t=3.

## 4. Reproduction and native guard boundary

Compile with the installed C++17 compiler:

    g++ -O3 -std=c++17 search.cpp -o search
    ./search 4611841 3 44591 570906 < input.txt > hits.jsonl 2> search.log
    python build.py
    python verify.py
    python -O verify.py

The final two runs use explicit exception guards and byte-identical printed results. They verify tournament orientation, full strongness, actual root ports, unchanged transitive prefix and cross edges, exact changed old pairs, physical name retention, all directed joint terms and all inner defects.

The direct ten-edge specification is not itself a native Hamilton-complement program: legality requires an ordered repair schedule, each actual Hamiltonian-path guard, and deletion of only fresh helpers. The independent audit must supply those certificates. Endpoint strongness does not silently certify intermediate native legality. The reusable open obligation remains a family closed under both physical digit maps; filling this finite missing44 does not settle it.
