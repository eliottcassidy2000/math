# Session 86: relations as attachment objects, F64/F256, and unrestricted tournament counts

2026-10-28. **Proved scoped carrier theorems, complete finite-field tests, independent verification, and an all-order enumeration formula.** The packet distinguishes labelled counts, orbit counts, stabilizer classification and actual class generation. No single small statistic is claimed to solve all four tasks.

## Main results

1. **F64 and F256:** Frobenius invariance, nondegeneracy and one marked non-subfield K-line determine the intended trace form. Two F64 line choices nevertheless retain a threefold frame symmetry.
2. **Positive joint-port reconstruction:** declared norm-one normalization gives closed multiplicative port carriers with 9 and 17 states. Enlarging the probe vocabulary saves two fixed, nonadaptive queries in F64, with exact lower bounds and witnesses.
3. **The 81/27 hostile is stronger than expected:** all nonempty named pointwise stabilizer groups agree literally, yet transporters and attachment orbits differ. Local isotropy is not joint transport.
4. **Closed attachment objects:** weighted joint relations, element-fixed boundary tables and subgroup-marked fibers give different exact carriers for different observations. Frobenius/LCA constraints repair correlations only when their actual equations force them.
5. **A000568:** an unrestricted formula recreates every order; independent computations match orders 0–50. Canonical composition is reconstructed through 20, and a precise coarsest cycle-index quotient is proved for universal rigid-palette inflation counts. Individual prime representatives and a class-faithful atlas are not thereby recovered.

## 1. “Relations become objects” has more than one operation

The prior 7→15 construction takes the four-dimensional relation space of Hamming 7 and makes its 15 nonzero vectors into new projective objects. Repeating that operation on Hamming 15, whose relation space has dimension 11, gives 2047 nonzero objects. It does **not** immediately give 63 or 255.

A concurrent [relation-tower packet](https://github.com/eliottcassidy2000/math/blob/20a31127389aaadabf0500b6dbf31b564abc9a13/05-knowledge/results/relation_tower_frobenius_ports_20261028/THEORY.md) makes the 2047 stage concrete with a supplied [prime bank](https://github.com/eliottcassidy2000/math/blob/20a31127389aaadabf0500b6dbf31b564abc9a13/05-knowledge/results/relation_tower_frobenius_ports_20261028/prime_bank.json) for Hamming[2047,2036,3]. This session independently checked its 4,085 recursive Lucas certificates, all 2047 distinct prime columns and all 698,027 integer square-total triples. The largest bank prime is 10,028,024,506,499. The 2^2036 valid subsets follow from rank-nullity; they were not exhaustively enumerated. This verifies the provided finite construction, not arbitrary-rank prime supply.

The F64/F256 experiments instead study binary dimensions 6 and 8, with quadratic extensions E/K. Their nonzero vectors factor as

\[
63=7\cdot9,\qquad255=15\cdot17,
\]

where q-1 counts nonzero scalars on a K-line and q+1 counts the lines. These are precise related constructions, not one unlabeled numerical sequence. Turning a relation into a new object requires stating its operations, coefficient domain, ports and observation.

## 2. Frobenius plus one marked line: exact tests in F64 and F256

[FIELD_RESULTS.md](FIELD_RESULTS.md) works with E=F_(q²), K=F_q, q=2^m and absolute Frobenius F(x)=x². The inherited all-m spread theorem identifies the trace form

\[
B_1(x,y)=\operatorname{Tr}_{E/\mathbb F_2}(x y^q)
\]

as the unique nondegenerate alternating binary form invariant under F and making every K-line isotropic. This session independently derives that classification and closes the requested **one-mark** tests at m=3,4.

|Field|Invariant alternating forms|Nondegenerate survivors|Minimum marked K-lines|
|---|---:|---:|---:|
|F64|8|2|1|
|F256|16|8|1|

Every K-line other than K itself identifies B1 among those nondegenerate survivors. K distinguishes none. A line is not one bit: its isotropy supplies multiple equations. Their independent ranks are 2 on six F64 lines and 1 on the two exceptional lines; every F256 non-K line supplies rank 3. In the exceptional F64 case, nondegeneracy excludes two remaining nonzero singular forms.

The actual finite-field models are verified irreducible and primitive. All invariant forms, restrictions and relevant linear stabilizers are exhausted with exact arithmetic. The all-m spread theorem is inherited/proved; the all-m **one-line minimum** is not inferred from these two cases.

### A unique form is not a unique frame

The joint stabilizer of Frobenius, the full spread and B1 is exactly the cyclic Galois group, of orders 6 and 8. A marked non-K line kills it in F256. In F64 six line choices kill it, but two short-orbit lines retain a subgroup of order 3.

Even with B1 fixed, relative port phase remains observable. In the declared F64 coordinates, F² sends 3 to 17 within 3K, while

\[
B_1(14,3)=0,\qquad B_1(14,17)=1.
\]

Applying Frobenius simultaneously to both ports preserves the form; moving only one port is a different alignment. Likewise in F16, individual Frobenius-orbit labels fail to distinguish port pairs(2,3) and(2,5), whose pairings are 0 and 1. Form recovery does not invent an unspecified relative transport.

## 3. A genuinely closed, normalized port carrier

[NORM_PORTS.md](NORM_PORTS.md) supplies a positive reconstruction mechanism. Retain the actual field/norm structure and declare that each port is the unique norm-one representative of its marked K-line:

\[
\nu(x)=x/\sqrt{N_{E/K}(x)},\qquad
U=\{u:N_{E/K}(u)=1\}.
\]

In characteristic two the square root in K is unique. Thus nu is invariant under K* scaling and induces a canonical group isomorphism E*/K*→U. The joint table is determined by

\[
B_1(u,v)=\operatorname{Tr}_{E/\mathbb F_2}(u/v).
\]

This does not recover arbitrary unknown scalars. It uses an explicit normalization rule and retained field data. If the line is already known, its normalized representative is computable without partner queries.

U has 9 states in F64 and 17 in F256, and is closed under multiplication, inversion and Frobenius. Its full partner-response rows are distinct, so no quotient identifying two U states preserves every named partner query. This is a sharp carrier statement for that grammar.

**Addition changes the closure requirement.** U spans E over F2. Allowing arbitrary port addition therefore enlarges the closed carrier to all 64 or 256 vectors. The 9/17-state carrier cannot silently stand in for an additive attachment language.

### Unknown-port query costs

Now suppose the normalized port u is unknown and answers are B1(u,w). The exact worst-case costs depend on both the probe alphabet and adaptivity:

|Field|Normalized ports|Fixed probes in U|Adaptive probes in U|Arbitrary linear probes, fixed sets suffice|
|---|---:|---:|---:|---:|
|F16 control|5|4|4|3|
|F64|9|6|5|4|
|F256|17|5|5|5|

The normalized pairing graphs are respectively K5, three disjoint K3 graphs, and the undirected Paley graph 17. Paley 17 is not the directed Paley 7 tournament from the previous session.

For F64 the three K3 blocks are orthogonal binary planes. A fixed partner set needs two probes per block, so six is sharp; there are 27 minimum sets. Adaptivity reduces the worst case to five: probe once in each block, then use at most two follow-ups. Four cannot suffice because after two adversarial zero answers at least five candidates remain. Allowing probe vectors outside U gives four probes, which meet the information bound. All 651 two-dimensional kernels were checked; exactly 36 avoid every difference of distinct U states. An explicit four-probe set is{3,12,17,33}. F256 has 714 minimum five-partner sets.

The six-to-five saving uses adaptivity; the further five-to-four saving admits richer probe objects. Neither creates information. If a compound probe must itself be implemented using several paid primitive measurements, that implementation cost must also be counted.

The same concurrent packet also studies fixed named **normal Frobenius-orbit** ports modulo simultaneous phase. That is a different candidate set and attachment language from the norm-one U ports here. Its phase/profile carrier and query counts must not be substituted into this table merely because both use Frobenius.

## 4. The 81/27 control: identical isotropy, missing transporters

[81_27_TRANSPORTERS.md](81_27_TRANSPORTERS.md) fixes the actual labelled pair

\[
T=C_3[A_0,A_0,A_0],\qquad T'=C_3[A_0,A_0,A_1],
\]

where A0 is a singleton source above a directed triangle and A1 is its literal converse on the same labels.

The pair has 67 identical labelled cyclic triples, H=462915, and identical Hamilton/path-cover profiles on all 4096 induced subsets. Nevertheless its automorphism groups have orders 81 and 27.

The stronger statement is literal: if G=Aut(T) and H=Aut(T'), then H is the within-block subgroup of index 3, and every element of G outside H moves all vertices. Hence

\[
G_D=H_D\quad\text{for every nonempty named subset }D,
\]

where the subscript denotes **pointwise** stabilizer. This does not apply to setwise stabilizers.

Yet the transporter taking vertex 0 to vertex 4 has 27 elements in G and none in H. Ordered distinct-port orbit counts are 2 versus 6 for one port,12 versus 36 for two, and 72 versus 216 for three. One-marker attachments can even have the same stabilizer at every named port while having different numbers of attachment types.

An explicit orientation-sensitive red-marker attachment at port 8 changes Hamilton counts to 1,395,279 versus 1,061,775. The source notes keep their converse and labelling conventions separate; that named-port statistic is not transferred to a differently labelled/conversed implementation.

There is a literal connection with the F64 experiment: its norm-one pairing graph is the same three disjoint three-point blocks. Choosing cyclic orientations for those blocks and their quotient, then adjoining one source/sink marker to each block, realizes the two twelve-vertex constructions above. The form supplies the blocks; the orientations and marker attachments are additional declared data, not consequences of the form.

This is why joint port data cannot be replaced by marginal stabilizer numbers or even all nonempty pointwise stabilizer groups. The missing information consists of maps **between** ports and the coherent global action.

## 5. What the smallest closed attachment carrier means

[CLOSED_ATTACHMENT_CARRIERS.md](CLOSED_ATTACHMENT_CARRIERS.md) specifies the operation language before minimizing. At finite named boundaries:

- For labelled weighted counts, retain the joint boundary-response table. Glue by compatible join and sum the products of witness weights over the shared assignment.
- If arbitrary named singleton probes are allowed, that table is the coarsest exact response carrier. Without those probes, a coarser quotient may suffice.
- With a fixed finite group and equivariant attachments, element-fixed boundary tables compose and suffice for final orbit counts via Burnside.
- Full subgroup fixed-fiber tables classify the actual finite fiber actions and stabilizers. They are stronger than necessary for orbit counts alone.

The odd-order C3² hostile makes the last distinction explicit: a regular nine-point orbit plus three fixed points and a union of four three-point coset orbits have the same element-fixed character and four orbits, but different global fixed points and stabilizer distributions.

These are closed relations for their declared operations. Exposing a forgotten port, changing the interpretation of complement, introducing new named probes or dropping a square/positivity/injectivity guard can enlarge the language and invalidate a smaller carrier. Native square-sum, graceful-tree and tournament guards are not assumed to be affine merely because the field model is linear.

For a genuinely affine model, the previous LCA quotient gives an exact criterion. Let W be all required rows and W_local the span retained by local scopes. Frobenius rows repair the join exactly when the supplied rows span the lost quotient, with compatible right-hand sides. A consistent tree arrangement does not create a row that was never retained.

Finally, invariance of a **family or distribution** is not pointwise invariance of every configuration. Three independent binary decorations have a cyclic-invariant uniform law on all eight choices. Equations requiring each decoration to be fixed under cyclic rotation reduce that to two choices; those equations are extra constraints, not a free consequence of the symmetric law.

## 6. Full A000568 recreation and the composition quotient

[A000568_COMPOSITION.md](A000568_COMPOSITION.md) gives an unrestricted all-order formula with no tournament-class list. For an odd partition lambda=(l₁,...,lᵣ), put

\[
e(\lambda)=\sum_i(l_i-1)/2+\sum_{i<j}\gcd(l_i,l_j),
\qquad z_\lambda=\prod_d d^{m_d}m_d!.
\]

Then

\[
a_n=\sum_{\lambda\vdash n,\ \text{all parts odd}}
\frac{2^{e(\lambda)}}{z_\lambda}.
\]

Even-cycle permutations reverse an antipodal edge and fix no tournament. Independent pair-orbit/Burnside calculations verify all 31583 odd-partition cases through 50; all 51 counts agree with OEIS. This is an independent recreation of a classical formula, not a new counting record.

The canonical grammar is also explicit: every tournament is an ordered sequence of strongly connected components, and a nontrivial strong component has a prime quotient with decorated children. The full aggregate prime cycle index is recovered and recomposed through 20, matching 129 coefficients and the ordinary counts. **Recovering prime data from Burnside and composing it back is a structural round trip**, not an independent discovery of the same numbers.

### A precise minimal quotient for ordinary composition

For prime-root inventories of arities 3 through K, equality of aggregate cycle indices is equivalent to equality of ordinary inflation counts for every palette of distinct nonempty rigid transitive children. A single palette with sizes

\[
1,(K+1),(K+1)^2,\ldots,(K+1)^{K-1}
\]

already separates unequal indices: independent color multiplicities become carry-free digits of the total size. Multiplicity tags, when present, are retained rather than silently deduplicated. The proof gives a coarsest **ordinary-count observer**, not a minimum bit encoding or a practical tiny palette.

For actual class generation, the rooted metagraph still needs aligned parent identifications and attachment masks. Its root-switch matrix satisfies Q²=NQ and, for N>=1,

\[
a_N=\operatorname{tr}(Q)/N=\operatorname{rank}_{\mathbb Q}Q
=\operatorname{rank}_{\mathbb F_2}Q.
\]

Counts alone cannot regenerate its rows, individual prime representatives or transporters. The empty order a₀=1 is handled separately. An unrestricted prime-class atlas and a smaller class-faithful attachment grammar therefore remain separate goals.

### A separate constructive bank, with its audit boundary retained

The concurrent [root-switch implementation and saved data](https://github.com/eliottcassidy2000/math/blob/20a31127389aaadabf0500b6dbf31b564abc9a13/05-knowledge/results/relation_tower_frobenius_ports_20261028/THEORY.md) reports 191,536 order-nine representatives from 1,716,608 pointed states, with trace 1,723,824. Our bounded check confirmed the saved bank contains 191,536 valid, literal-distinct tournaments and matches its recorded hash; the manifest arithmetic is consistent with 9*191536 and the labelled mass 2^36. We did **not** independently certify isomorphism uniqueness/coverage, the rooted-state census, those automorphism orders or matrix ranks. The source's independent native-metagraph histogram comparison stops at 8; the order-nine prime-class count and Hamilton-path-orbit histogram were not computed there. This is a constructive complement to the independently verified count formula, not a reason to erase the different verification scopes.

## 7. What unifies these results

The constructive hierarchy is now sharper:

1. **Choose the object language:** field-normalized ports, arbitrary labelled attachments, equivariant witnesses or canonical tournament inflations.
2. **Retain the composition data:** a group element, joint table, fixed-fiber response or cycle index, as required by that language.
3. **Recover what constraints actually determine:** the trace form, normalized joint pairing or missing affine row.
4. **Keep residual transport explicit:** the 81/27 pair and relative Frobenius phase show exactly where local summaries stop.

This is progress toward closed minimal carriers with proofs and hostile controls. It is not a universal finite summary, a proof of the graceful-tree conjecture or square-sum Hamiltonicity, or a claim that Frobenius alone reveals unnamed attachments.

## Reproduction

Five public Python checkers:

- `recreate_a000568.py` (standalone; embedded OEIS values are checked only after construction);
- `check_hostiles.py`;
- `check_attachment_carriers.py`;
- `check_field_forms.py`;
- `norm_ports.py` (imports the co-located `check_field_forms.py`).

They use only Python's standard library and the stated co-located helper. Run normally and with `python3 -O`; explicit checks remain active. JSON receipts are generated outputs, not required inputs. Independent implementations audit the field forms and stabilizers, port-probe minima, attachment fixed-point laws, the literal 81/27 permutation groups, and the enumeration/composition calculations.
