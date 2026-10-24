# Independent audit: native affine ternary phase codes

24 October 2026. **PASS after local prose corrections, in the stated scope.** The audited producer is `inputs/producer/AFFINE_PHASE_CODES.md`. Its all-size conclusions follow from the proofs reconstructed below, not from extrapolating finite checks. `independent_check.py` imports no producer code and uses standard-library integer arithmetic only. All checks remain active with Python optimization.

## 1. Scope, inheritance, and correction

The starting mechanism is the frozen literal C3 kernel carrier in `inputs/base/RELATIONAL_KERNEL_HOLONOMY.md`, together with the endpoint-marked gluing and workspace sidecars. The independent audit of that carrier is retained in `inputs/base-audit/AUDIT.md`. The closest hostile examples are duplicate resource keys, incompatible named observer successors, and higher-arity guards whose binary projections lose their constraint.

This extension explicitly supplies the predicate A k=b on canonical phases extracted from actual suffix owners by k=owner(q)-q modulo3, not on raw owner residues with an unchanged right side. It does not derive arbitrary linear equations from pairwise owner equalities. In particular the three-coordinate relation x+y+z=0 has full binary projections. The L=2 helper family likewise has 27 states but every two-coordinate projection is all nine pairs, so its binary-projection closure has 243 states. No conjunction of binary restrictions on those same five variables can recover the 27-state relation. This statement does not exclude an auxiliary-variable construction under separately authenticated rules.

The final producer retains the following necessary gates:

- Fixed distinct prefix keys are normalized capacity units. Three phases are three full-record lifts of one unit, never three units.
- The entire other retained guard must be constant on the family, or explicitly kept. Fixed source/target class incidence does not establish arbitrary owner-sensitive electrical constraints.
- Packet edges join accepted arrivals atomically. The result is not a microscopic routing theorem for an older metagraph.
- The source/target matrices and endpoint maps authenticate individual records. Equal target classes do not give a common named partner or one global observer.
- Normalized-resource non-reacquisition and per-slot full-record non-return are different predicates.
- The all-size designed graph is not the full budget graph once L>=3, and its degree and move budget grow.

A final precision pass also corrected wording that described canonical-phase equations as equations on raw terminal owner values. All existing producer and independent computations already extracted owner(q)-q correctly; their accepted states, spectra, and route conclusions do not change. Both packets now explicitly retain the raw-owner right-side shift q A1 and separate it from the converse-head phase shift2A1. An added independent control covers135 two-slot one-row systems for r=1,...,5, including72 nonzero raw-coordinate translations. At r=1, row(1,0), and canonical right side0, the correct raw right side is2: leaving it at0 incorrectly selects first phase1 instead of0. At r=2, the raw shift is1 while the converse-head phase shift is2, so the two translations cannot be silently conflated.

One sentence in the developing §6 suggested that extra higher-cost mode bridges could enter the full budget-(L+1) graph. That is false at this exact budget: mode-changing weights are L+1+j with j>=0. The producer corrected the sentence before freeze. The missing edges at budget L+1 are exactly same-mode changes in between two and floor((L+1)/2) pairs. Higher-cost bridges enter only at larger budgets. No theorem, executable result, or stated spectrum changed.

## 2. Native record and resource authentication

Let q=8r. The source has r ordered rigid eight-owner blocks followed by a directed C3. The actual prefix paths are P0=(3,0,6,7,1,5,4,2) and P1=(3,1,6,7,0,5,4,2); the seed matrix has upper-triangle code 29558785. A pinned key chooses a prefix and retains exactly the three paths with terminal sequence q+k,q+(k+1 mod3),q+(k+2 mod3).

The phase extraction k=owner(q)-q modulo3 is literal. If o is the vector of raw owners at clock q, then o=q·1+k, so A k=b is equivalent to A o=b+q A1 modulo3. Keeping the same right side on raw owner residues is generally wrong. This raw-coordinate translation is distinct from the converse-head phase shift2A1 discussed below. Raw owner arithmetic uses the declared suffix numbering q,q+1,q+2; under another owner embedding the phase chart must be transported explicitly, since arbitrary owner labels are not automatically F3 coordinates. Pairwise difference rows or r divisible by3 can conceal the translation because q A1 then vanishes. Distinct phases alter all three suffix clocks and no prefix clock. Suffix rotations are source automorphisms and actual partner transports. The map g=(4,0,2,3,1,5,6,7), or its inverse, transports prefix partners while preserving both actual block endpoints. It therefore preserves the exceptional retained junction edges in the ordinal construction. This is the necessary marked gluing datum.

The inherited all-size proof identifies the ordered strong components, the rigid prefix, the terminal C3 automorphism group, and the four-way source/K/converse reservoirs. The independent finite checker reconstructs the entire path-coordinate upper triangle, rather than the producer's shortened word, and confirms:

- complete prefix-pinned menus and literal Hamiltonian paths;
- suffix source and partner transports;
- both marked endpoints under every single-prefix-key partner change;
- four words per normalized key and disjoint reservoirs across keys;
- all literal converse K squares and exact three-clock support;
- nonexistence of a same-source suffix reflection;
- seed rigidity by all 8! permutations;
- the named-partner obstruction from prefix observers beginning 3,0 and 3,1.

These checks cover r=1,...,5, up to source order43, beyond the producer's r<=4 native range. Actual fixed keys may be reconstructed from the sidecars at every affine state; the code coordinate never manufactures a vacancy or releases a unit.

## 3. Affine state space, low-support components, and routing

After capacity and any constant guard pass, consistency of A k=b is essential. An inconsistent system has no accepted state and zero graph components. Otherwise a solution x0 gives the exact affine state set F=x0+C, C=ker A. Completeness of the three-record menus proves both directions of the correspondence. With d=dim C there are 3^d actual assemblies.

For x,y in F, the difference v=y-x has support exactly the changed slots; the actual clock cost is 3 wt(v). Thus the based endpoint mismatch enumerator is W_C(z t^3). It is a Hamming enumerator, not a smaller-budget shortest-path distribution. Empty fibres have no based enumerator; a singleton has W=1.

For nonnegative integer B let S_B be the nonzero kernel vectors of weight at most B and H_B their F3 span. Every packet path stays in an H_B coset because its increments lie in S_B. Conversely choose an independent subset of S_B spanning H_B. A representation of y-x in that basis is realized by applying each nonzero coefficient times its original vector. Each partial sum remains in C, so every arrival remains in the exact full affine guard and uses the same fixed resources. Coefficient2 has the same support as coefficient1. This proves exact components x+H_B, component count 3^(d-rho_B), size 3^rho_B, and at most rho_B events per connected endpoint pair.

For an owner-clock budget T the correct slot budget is floor(T/3). A full connectedness threshold is the least B with rho_B=d, equivalently the least maximum support among kernel bases. An endpoint threshold only requires its difference to lie in H_B. Neither is determined by dimension or minimum distance alone.

The fixed-B candidate enumeration visits exactly sum_(j=1 to min(B,n)) binom(n,j)2^j nonzero ambient vectors. Retaining those in C and selecting an independent ORIGINAL subset is polynomial for fixed B. Elimination may be used to test membership and recover coefficients; it may not replace the move vectors by arbitrary reduced rows. The audit includes the explicit hostile rows

    u=(1,1,1,0,0), v=(1,0,0,1,1).

Both have weight3, while their reduced second row u-v=(0,1,1,2,2) has weight4. Using that reduced row as a budget3 move would be wrong. This distinction is implemented correctly in the frozen producer.

## 4. Full-budget Fourier spectrum

On the additive group C, a character chi is an eigenfunction of the translated Cayley graph with Laplacian eigenvalue

    |S_B| - sum_(g in S_B) chi(g).

Translation of the affine origin only multiplies character functions by constants and does not change the spectrum. Since S_B=-S_B, the counts of omega and omega^2 evaluations coincide. If each is c, the eigenvalue is 3c. All eigenvalues are therefore nonnegative integers divisible by3.

It vanishes exactly for characters trivial on H_B. The annihilator of H_B in the dual group of C has cardinality 3^(d-rho_B), matching the connected-component count. The characters are those of C, not all ambient characters counted repeatedly: equivalently ambient characters must be quotiented by C-perp. The implementation parametrizes characters using a true basis of C and obtains each exactly once.

The independent checker calculates all such spectra for all ternary subspaces through length4 and all budgets0,...,n. It compares the first four spectral moments against graph-derived counts. In particular a regular graph with degree d and N vertices has traces N, Nd, N(d^2+d), and N(d^3+3d^2-t), where t counts ordered generator pairs g,h with -g-h in S_B. This calculation uses literal adjacency independently of the Fourier formula. Whenever the code has at most27 states, every cyclotomic character column is additionally multiplied by the entire literal matrix using integer pairs for 1,omega,omega^2. This includes the full 27-state example and certifies its exact eigenspaces, not numerical eigenvalues.

## 5. Exact no-return criterion and hostile separation

For a per-slot full-record no-return route from x to y, a coordinate with x_i=y_i cannot move: any departure would require returning to that initial record. If x_i!=y_i, the only possible sequences are the direct change or a visit to the unique third phase followed by the target. Writing delta_i=y_i-x_i, the third phase is x_i-delta_i. The two increment sequences are delta_i or (-delta_i,-delta_i), because -2 delta_i=delta_i in F3.

Assign stages0,1,2 to initial, third, and final phases at every differing coordinate. Retain precisely the stage vectors satisfying the original full affine equations and fix all equal-endpoint coordinates. A legal no-return packet is precisely an increasing edge between two retained stage vectors, with weight the number of changed coordinates. Thus no-return routes correspond in both directions to directed paths in this finite DAG. Bottleneck shortest paths give the minimum budget; shortest event paths after a budget cutoff give minimum event count. There are at most3^|supp(delta)| vertices and at most2|supp(delta)| nontrivial events in any path. No polynomial-time claim follows.

The direct endpoint packet always exists, with budget |supp(delta)| and minimum possible total clock work 3|supp(delta)|, absent a smaller imposed budget. Under a fixed smaller budget the DAG must retain its stage restrictions; unrestricted span membership alone is not enough.

For v=(1,0,0,1,1), w=(0,1,1,0,1), and C=span{v,w}, the only nonzero weight<=3 vectors are ±v,±w. Budget3 therefore connects all nine states, whereas budget2 has no edges. The target v-w=(1,2,2,1,0) has fifth coordinate equal to its initial zero. Every nonzero budget3 first step changes that fifth coordinate. Any route to the target must eventually return there, violating strict per-slot no-return. The ordinary route 0,v,v-w visibly has fifth phases0,1,0. A direct weight4 packet reaches the target without a return, proving sharp no-return budget4. The normalized resources remain continuously occupied on both routes.

The checker independently compares a literal visited-phase-mask breadth-first search with the stage-DAG criterion, for every zero-to-codeword endpoint, every code through length4, and every budget0,...,n: 13,071 endpoint-budget tests. Translation makes zero-based tests sufficient for arbitrary affine origins. This expands the producer's 811 no-return tests through length3. The five-resource hostile matrix is additionally checked on all243 ambient phase vectors and reconstructed as actual records in S_3.

## 6. Pair/helper family and exact designed gap

For L>=1 impose b_j-a_j=m in F3 on the canonical phases for all j and keep all three helper phases. In raw owner coordinates the same equation is o_bj-o_aj-o_m=-q modulo3, not generally0. The parameters (a_1,...,a_L,m) give dimension L+1 and 3^(L+1) states on 2L+1 distinct resources. At fixed m, each independent pair change costs2. Any change in m changes the helper and at least one endpoint of every pair, so it costs at least L+1; choosing delta a_j=0 or -delta attains that cost for each pair independently.

At m=0, each pair contributes 1+2z^2. At either nonzero m, the helper contributes z and a pair contributes 2z+z^2. Hence

    W_L(z)=(1+2z^2)^L+2z(2z+z^2)^L.

For L>=2, budgets below2 give singleton components, budgets2 through L give exactly the three fixed-mode components, and budget L+1 connects them. L=1 is correctly exceptional: its first bridge already has cost2.

The designed graph retains only the 2L single-pair moves and all 2^(L+1) minimum mode bridges. Thus d_L=2L+2^(L+1). It equals the full budget-(L+1) graph only at L=1,2. For L>=3 the latter also has same-mode multi-pair moves. At L=3 the designed/full degrees are22/34; at L=4 they are40/64. Degree and normalized gap must never be transferred between those graphs without recalculation.

For a character indexed by (r,t), write s=wt(r) and h=t+sum r_j modulo3. The bridge sum is

    (-1)^s 2^(L-s) times (2 if h=0, otherwise -1),

because 1+omega^(-delta r_j)=-omega^(delta r_j) at nonzero r_j. Together with within-mode adjacency 2L-3s this proves the producer's entire spectral multiset.

To minimize the positive Laplacian eigenvalues, at odd s the favorable branch is 3s+2^(L+1)-2^(L-s), minimized at s=1. At positive even s it is 3s+2^(L+1)-2^(L-s+1), minimized at s=2 and larger by3. The other branches are larger. The s=0 nontrivial characters have eigenvalue3·2^L, equal to the s=1 minimum only at L=1. Consequently the gap is 3+3·2^(L-1), with multiplicity4L for L>=2 and6 for L=1.

The sharp normalized bound follows from the exact identity

    40 gap_L - 27 d_L = 6(2^L-9L+20).

The parenthesis is positive at L=1,2,3, zero at4, and strictly increasing thereafter because its increment is2^L-9>0 for L>=4. Thus the normalized gap is at least27/40, exactly at L=4. This is an all-size proof, not an extrapolation from the checker's L<=256 arithmetic controls. It retains the growing slot budget L+1 and exponential degree.

At L=2 the full budget3 graph has 27 states, degree12, and exact spectrum 0^1,9^8,12^6,15^12. The independent literal-matrix character test confirms this, separately from the producer's split-polynomial and trace checks. Designed spectra and weight enumerators are independently enumerated through L=6, beyond the producer's L<=5 range.

## 7. Coordinates, converse, and endpoint data

The mathematical substitution k=D P z+c transforms the predicate into A D P z=b-Ac. An affine monomial map preserves every Hamming distance; a general invertible linear map does not. The audit tests all1,296 permutations, nonzero coordinate scales, and shifts at length3, against the full guard and every ordered pair of ambient vectors. The invertible shear (u,v)->(u,u+v) gives a weight1-to-weight2 hostile example.

This is coordinate equivalence, not automatically a symmetry of native owners or of the unchanged guard. Literal rotations have actual source and partner maps; slot permutations must carry keys and records with them. A directed C3 reflection is not a same-source automorphism. Reflections or independent incidence signs cannot be inferred from F3 notation.

Converse transposes the source and reverses the literal path. The old suffix clock q+t becomes2-t. At converse clock2 the original k is recovered unchanged. At converse clock0 the extracted phase is ell=k+2, so A ell=b+2 A1. The changed right side is necessary unless that shift happens to vanish. The K square is checked on literal matrices; there is no guessed minus sign. The varying part lies at the front of the reversed path, rather than being silently identified with the original terminal suffix.

## 8. Finite universes, reproducibility, and limits

Subspaces are independently enumerated by adjoining ambient vectors and closing under addition, not by the producer's RREF catalogue. For lengths0,1,2,3,4, the counts are1,2,6,28,212. All affine cosets are then enumerated, giving1,4,22,184,2452 fibres respectively. This includes the zero-dimensional and zero-slot boundaries. The audit checks13,071 affine budget graphs,1,195 kernel Fourier graphs, and13,071 literal-history/DAG endpoint-budget tests.

A separate exhaustive rectangular-system pass covers every A over F3 with0<=n<=3 and0<=p<=2, including zero rows and repeated/dependent rows, and every right side:7,504 systems,904 inconsistent. This supplements the code/coset census rather than assuming the latter represents inconsistent inputs.

The package includes byte-frozen producer inputs and all needed proof sidecars. During final packaging, historical CUBES.md snapshots were renamed CUBES.txt without changing their bytes, so historical relative Markdown links are not presented as live package links. Original repository paths and source digests remain retained; only local archive paths and dependent packaging metadata changed. All affected normal/optimized and manifest-only relocated replays were repeated, with mathematical receipt fields unchanged. Normal and optimized independent receipts agree exactly. Producer normal and optimized replay is recorded separately. A manifest-only temporary-copy replay checks both independent modes, producer replay, and strict byte-exact closure from outside the original workspace. The root manifest excludes only itself and bytecode caches; unexpected files fail closure. Paths and byte hashes are retained in `INPUT_PINS.json`, and replay commands are in `README.md`.

These conclusions are **PROVED within the declared model** or **FINITE-EXACT for the stated universes**. They do not prove pairwise-equality universality, a global observer, arbitrary reflected native holonomy, phase-sensitive electrical neutrality, unrestricted Hamiltonian-path connectivity, microscopic metagraph routing, or fixed-locality all-size mixing.
