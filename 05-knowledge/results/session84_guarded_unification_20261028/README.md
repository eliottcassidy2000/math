# Session84: guarded observations, divisor balance, and Berggren addresses

2026-10-28. **Proved applications, independently checked computations, and a proposed organizing framework.** Standard context-congruence/rank theory and inherited number-theoretic results are credited separately. This packet proposes no new universal conservation law, no cross-domain isomorphism, and no theorem that the repository's discoveries share a dynamical attractor.

## The larger framework

The strongest common structure is an **object equipped with named operations, legality guards, observations, and a scale**. The central question is:

> Which information may be forgotten now without changing the answers to permitted later experiments?

The prior [observation/composition reflection](../../../07-reflections/why-the-advances-work-observation-composition-and-non-backtracking-macmini-20261026.md) already points toward this framework. It is a source of questions, not a theorem overriding corrected canon. The [THM-840 operation-kernel theorem](../../../01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md) supplies an earlier exact mechanism. This session makes the common framework precise in several arithmetic and tree settings.

Three different reductions must remain distinct:

1. **Histories to objects.** Divisor multiplication merges commuting histories; the native Berggren tree has unique addresses. A divisor diamond is not a tree merely because both have recursive constructions.
2. **Objects to observations.** Divisor counts, weighted sums, residues and signs discard different information. Equality of one observation does not transport every other predicate.
3. **Observed objects under future operations.** Some observations are closed under the chosen updates; others need extra coordinates. Changing the permitted operation language can change the exact minimal quotient.

This is a proposed organizing framework for the research, grounded by the theorems below. It is not a claim of priority for contextual equivalence or a replacement for domain-specific proofs.

## 1. F=U+S: recover the exact object first

[FUS_DIVISOR_FRONTIER.md](FUS_DIVISOR_FRONTIER.md) uses the original counts on D(N)={d|N:1<d<N}, N>=2:

\[
F=\tau(N)-2,\quad U=\omega(N)-[N\text{ prime}],\quad
S=2^{\omega(N)}-1-[N\text{ squarefree}].
\]

The prime set is contained in the squarefree set, so F=U+S is a count equality, **not a disjoint partition**. S here is not the divisor sum sigma. The balanced shapes are exactly p,p³,p²qr, with distinct primes in the mixed shape. Their proper-divisor term counts including1 are1,3,11; that list is not the prime-label set{2,3,11}.

Let δ=F-U-S. The exact update laws prove that δ is monotone on nonsquarefree integers under divisibility, with just one nontrivial equality family p² | p²q. Consequently balanced composites form a divisibility antichain: every proper multiple has positive defect.

The positive-defect region has four minimal prime-labelled divisibility profiles:

\[
p^4,\quad p^3q,\quad p^2q^2,\quad p^2qrs.
\]

Balance is therefore a thin equality locus, not a conserved quantity. It is not even a separator every path must visit:4→12→24 has defects-1,-1,+1 and jumps over zero.

## 2. Count balance and square divisor sums meet, but are not interchangeable

At66=2*3*11, sigma is144, but(F,U,S)=(6,3,6), so δ=-3. Promoting any one prime exponent gives132,198 or726, all balanced p²qr shapes. Their divisor sums336,468,1596 are nonsquares.

The elementary reason is that for v_p(N)=1,

\[
\frac{\sigma(pN)}{\sigma(N)}=\frac{p^2+p+1}{p+1}
\]

is a reduced nonsquare rational. Thus a square-sigma source necessarily loses that property under this promotion. This is not global incompatibility:364=2²*7*13 is the least mixed balanced square-sigma source, with sigma784=28².

There is a stronger exact result, **using the inherited classical Ljunggren prime-power theorem**: on square-sigma sources, multiplication by a prime preserves square sigma iff the prime is3 and does not already divide the source. Hence the prime-multiplication graph on all positive square-sigma sources consists precisely of pairs{N,3N}, for3∤N, and isolated sources with v_3>=2. The pair{1,3} is included in that auxiliary graph; removing1 also isolates3. This graph is not the aliquot graph.

The finite tests do not prove Ljunggren's theorem. Its precise inherited role is exposed in the note.

## 3. Complete abundance of the balanced family:650 is singled out

[BALANCE_AND_EXCESS.md](BALANCE_AND_EXCESS.md) gives a complete all-size abundance classification among balanced sources. The abundant set is the disjoint union

- 12r for prime r>=5;
- 20r for prime r>=7;
- 28r for prime r>=11;
- 18r for prime r>=5;
- 6p² for prime p>=5;
- and the seven exceptions350,490,550,572,650,748,836.

**No balanced source is perfect; every other balanced source is deficient.** In particular,

\[
F(N)=U(N)+S(N)\ \text{and}\ s(N)=N+2\quad\Longleftrightarrow\quad N=650.
\]

This recovers an exact link with20's alternate representation. The650 uniqueness is already present in the incoming overlap/dyadic proofs and is credited there; the full abundance classification above supplies its broader context. Concurrent [factor-pair overlap work](../aliquot_divisor_overlap_20261028/THEORY.md) proves that s(N)=1+a+b with ab=N is equivalent to either tau(N)=4 or s(N)=N+2. Thus the balanced sources admitting that compression are exactly prime cubes and650. The literal proper-divisor counts remain distinct:20 has five terms;650 has eleven. This does not reclassify20 as a literal three-term source.

The proof uses exhaustive prime cases, not extrapolation from a census. The checker independently compares literal divisor counts through100000, finding2971 abundant balanced sources and the unique excess-two intersection650.

## 4. The exact {2,3,11} and Berggren-root bridge

[BRACKETS_AND_66.md](BRACKETS_AND_66.md) recovers the original bracket meaning. In B_j=((2j-1)²,(2j+1)²], the primes whose doubles remain in their own bracket are exactly2,3,11. More generally, the same prime set is selected precisely when

\[
25/13<\lambda\le25/11.
\]

This is an exact parameter plateau, not an iterative attractor. Its prime product is66.

[BERGGREN_SEMIPRIME_BRIDGE.md](BERGGREN_SEMIPRIME_BRIDGE.md) supplies the second explicit route. Every even distinct-semiprime square-sigma source has

\[
N=2q,\quad q=12k^2-1\text{ prime},\quad
(q-2,2\sqrt{\sigma(N)},q+4)=3\Phi(2k,1).
\]

Here Phi(m,n)=(m²-n²,2mn,m²+n²), and Phi(2k,1) runs along the D-ray from the primitive root(3,4,5). At k=1 this gives q=11 and22→14. The coprime sigma-square multiplier3 then gives66→78. The bracket construction and this guarded root construction therefore agree exactly on66, with their operations and assumptions stated. They are not conjugate dynamical systems.

All lifted targets t=72k²+6 have **no class-II aliquot predecessor**, since the required even-target prime t-3=3(24k²+1) is composite. They need not have singleton fibers: at k=2, target294 receives both282 and293².

## 5. A complete decorated Berggren bijection, not a prime-preserving tree

For a primitive triple(a odd,b even,c), odd distinct-prime square-sigma sources have exactly two charts:

\[
(p,q)=(d(c-a)-1,d(c+a)-1),\quad d\ge1,
\]

or

\[
(p,q)=(d(c-b)-1,d(c+b)-1),\quad d\ge2\text{ even}.
\]

Both endpoints must be prime. Retaining the primitive tree address, scale and choice of half-gap leg gives a bijection onto all such odd sources. Dropping any of those coordinates changes the statement.

The exact child formulas retain one prime endpoint but require retesting the other. An invalid node can have a valid descendant: the odd-gap root at scale1 gives(1,7), while its U-child gives(7,17), hence119→25. A failed prime guard cannot justify pruning that entire subtree.

For each fixed factor gap q-p=2r, the candidates reduce to factor pairs(c-z)(c+z)=r². Hence a bounded gap gives only finitely many square-sigma semiprimes. There are no twin-prime examples; the least odd-prime gap is10 at7,17. These primes need not be consecutive, and this predicate differs from the older twin-prime/Berggren dictionary.

A complete source census through1000000 contains209867 distinct-prime products and901 square-sigma sources, split540 odd-gap,285 even-gap,76 even D-ray. Independent divisor sieving and a separate triple-matrix traversal verify the results; the tree/Euclid comparison covers15919 primitive triples through hypotenuse100003. No quadratic-prime infinitude is claimed.

## 6. Two kinds of memory that must not be confused

[BRANCH_MEMORY.md](BRANCH_MEMORY.md) starts from the audited Berggren branch-character theorem. For every modulus M, the root and the word A D^(M-1) U D have the same endpoint residue and opposite determinant characters. Thus no fixed endpoint modulus determines the character.

Nevertheless **one extra bit is exactly sufficient**: A flips it, U and D preserve it. Every reachable endpoint residue occurs with both signs, so the joint reachable observer is exactly E_M×{+1,-1}. Endpoint dynamics were already closed; the sign merely did not descend through that quotient. The full branch matrix modulo3 also determines the sign. Projective pole data and a normalized pole vector are not interchangeable.

This corrects an overbroad reading of “growing memory”: failure of a chosen residue observation does not prove that the whole branch word must be stored. The bit is finite-state even though the tree is infinite.

## 7. A precise unifying theorem and sharp arithmetic applications

[GUARDED_OBSERVERS.md](GUARDED_OBSERVERS.md) proves the standard context-congruence construction with explicit legality failures. Two objects may share a sufficient deterministic state exactly when every permitted continuation gives the same output and is defined in the same cases. The quotient is the coarsest exact carrier for that specified operation language.

The arithmetic applications sharply separate three meanings of “small”:

1. **Linear arithmetic dimension.** On coprime rectangles, s(ab)=sigma(a)sigma(b)-ab has sharp rational separation rank2. On coprime rectangles with nonsquarefree products the count defect kernel has sharp rank4, with feature vector(tau,2^omega,omega,1). Explicit minors certify both lower bounds.
2. **Complete guarded state.** Under arbitrary named coprime multiplications, the full-domain scalar FUS defect has exact quotient(rad(n),tau(n)). Prime-support guards require unbounded information even though the arithmetic update has few rational coordinates. For aliquot s, future observations distinguish every integer.
3. **Operation-language dependence.** Allowing named existing-prime promotions makes the FUS defect quotient distinguish every source N>=2 too. The pair972=2²3⁵ and288=2⁵3² shares(rad,tau)=(6,18) and current defect11; a promotion by2 gives defects17 and14. The richer language observes the previously hidden exponent assignment.

Rational coordinates can hold unbounded integers. A rank bound is not a bit bound, state-count bound or universal compression theorem. This distinction is the direct connection with the earlier tournament/shared-port/Vitali work: the allowed later experiment determines which correlations and names must remain observable. No arithmetic–tournament isomorphism is asserted.

## 8. What is established and what remains a research direction

The organizing framework is now backed by explicit quotients, lower bounds, operation laws, complete restricted classifications and hostile examples. The next questions are domain-specific: further prime-guarded intersections, richer tournament contexts, and exact carriers for their own operations. Neither recurring numerals nor a metaphorical “attractor” proves those extensions.

Open problems are not silently closed: the work does not prove quadratic-prime infinitude, the prime-untouchable conjecture, general aliquot termination, or a universal theorem equating arithmetic and tournament structures. The square-preserving prime-step result has a classical dependency; the unifying contextual theorem is standard theory with concrete applications, not an invented general theorem.

## Reproduce

The six public checkers are standalone Python standard library programs:

- `check_brackets_66.py`
- `check_branch_memory.py`
- `check_fus.py`
- `check_observers.py`
- `check_berggren_bridge.py`
- `check_balance_excess.py`

Run each normally and with `python3 -O`; explicit checks remain enabled. Some scripts emit deterministic JSON, which is generated output rather than an input dependency. Independent audits used different implementations for the principal counts, matrices, divisor formulas and guarded finite systems. All universal claims rely on the displayed proofs or explicitly cited theorems, not on finite tests alone.
