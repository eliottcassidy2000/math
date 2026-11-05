# A finite positive fixed-port interface for both birth directions

## 1. Scope and results

Fix a nonempty tournament and k distinct physical vertices as ports. Retain exactly those ports forever. At a right birth the newborn beats an arbitrary selected subset of the ports and loses to every other old vertex. At a left birth it loses to the selected ports and beats every other old vertex. No newborn is exposed, no original port is forgotten, and no other old vertex is re-exposed.

There is an exact finite nonnegative integral cover interface for arbitrary mixed words in this language. Its numbers of coordinates for k=0,1,2,3,4 are respectively

`1,8,45,253,1487`.

These are transparent interface sizes, not minimal dimensions. In particular, fixed two-port mixed growth has a 45-coordinate model. This does not contradict the infinite continuation rank proved for the larger language that permits port replacement.

For periodic words, a further all-order obstruction holds:

**Mixed periodic theorem.** Suppose a fixed nonempty word W is repeated from any actual initial graph with fixed ports. If its Hamilton counts eventually satisfy a nonconstant recurrence `h_(n+1)=3h_n+c`, then

- one birth direction occurs exactly once in W;
- `c=2 modulo4`.

In particular **no fixed-port mixed periodic word implements `3H+4`, at any finite number of ports, even eventually**. Port-changing words, state-dependent words and graph rebuilding remain outside this theorem. A mixed periodic `3H+2` example is given in Section 6.

## 2. Covers and exact two-sided closure

For an unordered directed path cover of the entire old graph, let c be its number of nonempty components. Let P be the number of port occurrences among its path endpoints, with a singleton port counted twice. Keep precisely the covers satisfying

`c-P<=1`.

Encode each component by its initial and terminal port labels, or `*` for a nonport endpoint; record the multiset of these pairs. A port cannot occur in two components. If it occurs at both ends of one component, that component is its singleton. Other unmentioned ports are internal vertices of their paths. There are at most2k port endpoint occurrences, hence at most2k+1 components, so finitely many shapes. Repeated `(*,*)` components are allowed.

The coordinate of a shape is the number of actual unordered covers having that shape. The Hamilton count is the sum of the one-component coordinates, which all satisfy the restriction.

Insertion selects an optional predecessor component and an optional successor component, joins them through the newborn, and may instead leave the newborn singleton. They cannot be the same component, which would create a directed cycle.

- Right birth: the predecessor may end at any nonport or unselected port; the successor must start at a selected port.
- Left birth: the predecessor must end at a selected port; the successor may start at any nonport or unselected port.

A missing predecessor or successor is always allowed. Choices among repeated shape components have their actual multiplicity. Discard a result failing c-P<=1. This is a nonnegative integral matrix for every mask and side.

For exactness, delete the newborn from a counted output cover. If it is internal, deletion increases c by one and restores at least one port endpoint: its successor is a port for a right birth, or its predecessor is a port for a left birth. Therefore c-P does not increase. If it is a path endpoint, c does not increase and P can only increase. If singleton, c decreases by one. Thus the recovered old cover is always retained. Deletion uniquely recovers the old cover and the two attachments, proving the exact transfer at every graph order.

One explicit enumeration of the shapes chooses disjoint port-containing components, then adds copies of `(*,*)`. If there are p components with two different port endpoints, a singleton-port components, and b components with exactly one port endpoint, their multiplicity is

`k! 2^b / (p! a! b! (k-2p-a-b)!)`.

There are p+a+2 permitted numbers of additional `(*,*)` components. Sum over nonnegative p,a,b with2p+a+b<=k, then subtract one to exclude the empty cover. This gives the displayed coordinate counts.

Capping every coordinate at X+1 commutes with these nonnegative matrices. If the Hamilton readout is at most X it is exact; if it exceeds X, induced-count monotonicity permits permanent rejection. Thus capped-state closure also gives exact arbitrary-depth target reachability in this fixed-port language.

## 3. Directional triangularity and opposite-side commutation

Let U=2k-P be the number of used path-cover incidences at ports. For right transfers use the potential

`Phi_R=U+(number of nonport starts)`.

Port incidences only increase, and an old nonport start cannot be consumed by a right birth. The potential therefore does not decrease. Its range is0,...,2k+1: since c<=P+1, the number of nonport starts is at most the number of port ends plus one.

Equal-potential transitions are exactly appending the newborn to one of the nonport ends, without another attachment. They keep the shape unchanged. Hence every right transfer is triangular in this potential, with diagonal d equal to the number of nonport ends, in0,...,k+1, independent of its mask. Left transfers have the analogous potential `Phi_L=U+(number of nonport ends)` and diagonal equal to the number of nonport starts.

A right-only word of length r consequently has eigenvalues among `d^r`, d=0,...,k+1. Its characteristic polynomial, on this one fixed finite interface, supplies an annihilator independent of the initial graph. Possible Jordan blocks give polynomial factors in n, not additional exponential bases. The analogous statement holds for a left-only word.

A right and a left birth commute up to interchanging their newborn labels, while fixing every old vertex and every port. The left newborn always beats the right newborn in either chronological order, and their incidences to the old graph are unchanged. This proves physical commutation and equality of their cover transfers. Within each side, the order of masks is preserved.

## 4. A mixed period separates into two one-sided words

Let W contain r right births and l left births. Let R and L be its two ordered side-subwords. The graph after n periods is isomorphic, fixing the old graph and ports, to the graph obtained by R^n followed by L^n.

The two period matrices commute. Commuting finite complex matrices can be simultaneously triangularized: choose an eigenvector in an eigenspace of one matrix invariant under the other, then induct on the quotient. Their product therefore has eigenvalues among

`d^r e^l`, with d,e in0,...,k+1.

Alternatively, the uniform directional annihilators give a tensor exponential-polynomial representation for the two independent repetition counts, with the same diagonal bases. If a direction is absent, use only the other one-sided spectrum, without a0^0 convention.

A nonconstant eventual affine ternary recurrence has the form

`h_n=A*3^n+beta`, with A>0 and beta=-c/2.

Thus3 must be a period eigenvalue. If both r,l are at least two, no product d^r e^l is three. If one side is absent, its nonempty side must have length one. Therefore at least one direction occurs exactly once in W.

## 5. The leading ternary amplitude is even

First, if only one birth occurs per period, the repeated newborns form one transitive clone module. For a fixed rooted source, let a_j count directed root skeletons with j root occurrences and no consecutive roots. The clone formula is

`H_m=sum_j a_j j! S(m,j)`.

The j-1 nonempty internal fragments between roots can be permuted freely. Hence `(j-1)!` divides a_j, in particular a_3 is even. A pure leading ternary mode has coefficient a_3, so its amplitude is even.

It remains to treat a genuinely mixed word. By reversal assume its left subword has length l=1 and r>=1. After n periods, the n left-born vertices form a transitive module over the source containing the initial graph and R^n. Every left-born vertex beats every right-born vertex and has the same incidences to all original vertices. Let a_j(n) be that source's root-skeleton coefficients. Then

`h_n=sum_(j=1)^(k+1) a_j(n) j! S(n,j)`.

The upper limit is UNIFORM in n: the left root loses only to its selected subset of the k fixed ports, and every internal fragment must end at a distinct one of those ports. Thus j-1<=k. Each a_j(n) is nonnegative, integral, and a_3(n) is even.

For fixed q, the count with q left clones and R^n is a right-only continuation sequence from a fixed actual source. It obeys the uniform right-word recurrence, whose roots are d^r. The triangular Stirling system at q=1,...,j expresses a_j(n) as a fixed rational combination of those counts. Hence each a_j(n) is an exponential-polynomial sequence with bases d^r, uniformly in n.

Now suppose `h_n=A3^n+beta` eventually. Positivity of every term in the run sum gives:

- For each of the finitely many j>=4, `a_j(n)=O((3/j)^n)`, so it is eventually zero.
- a_3(n) is bounded. Its only bounded nonzero recurrence base is1, so it is eventually a constant even integer a.
- `a_2(n)=O((3/2)^n)`. Every nonunit base d^r is at least2, so a_2(n) is eventually polynomial in n. Its2^n contribution is o(3^n).
- If r>=2, `a_1(n)=O(3^n)` and its first possible nonunit base is at least4. It too is eventually polynomial and contributes o(3^n). Hence A=a is even.
- If r=1, a_1(n) is itself the ordinary clone sequence for the one right-born vertex, with the one left root fixed. Its growth bound excludes run indices above three. Its3^n coefficient b is even by the same free-fragment action. Hence A=a+b is even.

This proves even amplitude in every case. Since all Hamilton counts are odd, beta must be an odd integer. Therefore `c=-2beta=2 modulo4`, proving the theorem. The eventual formulation is included throughout; finitely many initial terms are immaterial.

## 6. A sharp mixed +2 example

Use the actual eight-vertex seed

`T=[2,12,145,149,99,143,175,19]`, ports `[0,1,4]`.

Its root7 transitive clones have count60*3^m-1. Repeat the fixed-port two-birth word

`L(mask0), R(mask7)`.

Every left newborn is a universal source relative to all previous vertices. After commuting opposite-side births, the output after n periods is exactly a transitive prefix of n vertices followed by the root7 clone graph with n+1 clones. All ports are still the same original vertices. Consequently

`h_n=60*3^(n+1)-1`, and `h_(n+1)=3h_n+2`.

Thus mixed periods of length greater than one can return +2. The obstruction is specifically to the stated additive congruence, not to every nonzero ternary affine law.

## 7. Verification and limits

The supplied verifier checks all formal directional potentials and diagonals through k=4, opposite-side transfer commutation through k=3, and physical insertion against literal unordered-cover enumeration for every labelled graph through order four and every choice of up to three ordered ports. It also checks the mixed +2 graph identities and direct Hamilton counts. Normal and Python-optimized receipts agree.

The all-k closure, uniform spectra, bounded-run argument and free-fragment action are mathematical proofs; finite tests only validate their implementations and explicit witnesses. The underlying cover dimension is not claimed minimal. Nothing here rules out an adaptive nonperiodic fixed-port algorithm, a port-changing periodic algorithm, or a count-preserving graph replacement.
