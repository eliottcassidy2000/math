# Both next digits from the actual nonzero-defect multiport source

**FINITE-EXACT physical constructions and bounded searches; PROVED arbitrary-prefix isolation. No reusable two-digit closure or global spectrum theorem is asserted.**

## 1. Actual lineage, not numerical identification

Start with Session181's literal 17-name H=512425 tournament, root39, prefix50→51, physical suffix names

    [1,2,4,6,7,8,14,21,22,31,32,33,40,54].

Its root in-ports are {4,21}, suffix fields h=18583,F=124383,J=3690, and evaluator H_s=14893+124383·2^s. The source is the full adjacency in `source.json`, frozen into the independent native audit. This is NOT the different 18-vertex zero-ready H512425 source (46073,58294,3,port40) of `zero_ready_strong_switch_20261103`, which already has numerical children1537277/1537279. Our contribution is continuation of the actual S181 physical state and its nonzero defect, not first discovery of these numerical counts.

## 2. Two actual next endpoints

Both branches retain all17 names and birth prefix55 after50→51. The new root39 beats55, both existing prefix vertices beat55, and55 beats every suffix vertex. The final prefix is TT3.

**Digit2:** Reverse suffix pairs{2,6},{31,40}; change root pair{39,4}, leaving its sole in-port21. The suffix is strong, with

    h=80637, F_21=182080, J=0,
    H=80637+8·182080=1537277=3·512425+2.

The actual endpoint is zero-ready at39. This is a physically realized reset after rewrites and a birth, not a scalar substitution or an assumption that the old3690 defect disappeared.

**Digit4:** Reverse suffix pairs{14,32},{21,33}; change the root in-port set from{4,21} to{33,40}. Thus four root pairs change. The suffix is strong, with

    h=58063,
    F_33=127669, F_40=62633, F=190302,
    J(33,40)=37224, J(40,33)=5976, J=43200,
    H=58063−43200+8·190302=1537279=3·512425+4,
    E_39(2)=1580479, E_39(2)−H=43200.

All diagonal J terms vanish by squarefreeness. The two directed cross terms are retained separately. `candidate2.json`, `candidate4.json` contain every physical row and name. `verify.py` independently recomputes the suffix fields, exact joint terms, changed old edges, tournament identities and full strongness. It passes normally and under Python optimization.

The sibling `closure_audit` packet supplies independently checked native complement guards, a strong completed repair schedule and full-tournament subset-DP counts for both branches. The order of literal native operations is the audit's recorded order, not an arbitrary order of the pair lists above. Temporary helpers are deleted; no old name is deleted. Field observer coordinates must be read from the actual-state phase audit, not inherited from a numerically equal source.

## 3. Exact arbitrary-prefix isolation of these two transformations

Let A be any nonempty prefix, and use the same old suffix and ports. Replace A by A→{55}; this preserves H(A) and doubles P_A(2). Perform the respective suffix/root rewrites above. From the all-order multiport evaluator,

    H_out,2−3H_in =35958 H(A)−8989 P_A(2),
    H_out,4−3H_in =−29816 H(A)+7455 P_A(2).

Both transformations have their desired digits simultaneously if and only if H(A)=1,P_A(2)=4. Indeed eliminating P gives50866 H(A)=50866, then P=4. A tournament with H=1 is transitive: any nontransitive tournament contains a directed triangle with three Hamilton paths, and deterministic insertion of the other vertices injects these into its Hamilton paths. A transitive tournament of order a has P(2)=2^a. Therefore **the two displayed transformations jointly work exactly at A=TT2**, up to names. They are not a shared arbitrary-prefix induction mechanism.

In particular for A=TT_s the two offsets are35958−8989·2^s and−29816+7455·2^s, respectively. Their required offsets2 and4 each force s=2. Merely lengthening the same prefix does not iterate this step.

## 4. Complete bounded suffix menus, unbounded clocks

`search.cpp` considers every zero-, one-, or two-edge reversal of the actual14-name suffix:4187 labelled edits,3946 strong suffixes. Every nonempty W and every transitive clock s≥1 is allowed. There are25 numerical+2 hits (21 at s1,3 at s2,1 at s3) and26 numerical+4 hits (22 at s1,2 at s2,2 at s3). Both selected endpoints are among the s3 records, hence preserve and extend the old TT2 prefix. No hit satisfies either corrected shared-clock target(c,F)=(44681,373149) or(44683,373149).

`search_terminal.cpp` widens the suffix menu to all4187 edits. It requires W to meet the suffix's terminal strong component, which is necessary and sufficient for the full attachment to be strong. The terminal component is the intersection of the reachability sets from all vertices. This yields29+2 hits and28+4 hits, adding4 and2 nonstrong-suffix records. It still has no corrected shared-clock target. The selected stronger-suffix witnesses above remain unchanged.

Every clock is covered, not only clocks below a numerical cap. Since J_W≤F_W and F_W≥1,

    H(M(TT_s;B,W))≥h+(2^s−1)F_W.

Thus a positive fixed target bounds s and F. The loop bound32 is redundant for targets below2^21. No candidate at the corrected uniform target is lost by the pruning: h=c+J≤c+F<1537277. Counts and intermediate products at order14 fit signed64-bit integers; each two-slot count is bounded by15·14!, and tested clock values remain within range. The search does not classify arbitrary prefix interiors, more suffix rewrites, suffix births, or all possible physical transformations.

## 5. Reproduction

Compile the included C++ source with the installed compiler. For the strong-suffix menu run `search 1537277 2 < input.txt` and `search 1537279 2 < input.txt`; for the widened menu use `search_terminal` with the same arguments. JSONL files retain every record; summaries give exact census sizes. Run `python verify.py` and `python -O verify.py`. No installs or repository startup helper are used.
