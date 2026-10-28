# Session85: a bit from constraints, arithmetic codes, and higher-order relations

2026-10-28. **Proofs, complete finite realizations, independent checks, and explicit limitations.** The organizing tools are established mathematics: binary codes/matroids, Boolean coordinate rings, affine gluing and integral signed incidence. This session supplies concrete arithmetic realizations, recovery certificates and typed bridges between those tools. It does not claim that matching numbers alone unify all mathematics.

## The working synthesis

A useful step beyond “retain what later operations can observe” is:

> Write the relations that make some coordinates redundant, and produce a certificate that a missing coordinate follows from the retained ones.

A bit can then cost no **new measurement**, while still costing a valid constraint, a guard check, a symmetry-breaking convention or computation. We distinguish six mechanisms:

1. a dual parity equation reconstructs an erased bit;
2. nonlinear legality removes an otherwise possible alternative;
3. dependent checks save repeated work;
4. separator constraints avoid paying twice for shared information;
5. a gauge choice removes representational signs while retaining reconstruction data;
6. a finite typed correspondence transports a recovery rule to another domain.

These are constructive mechanisms, not information created from nothing. Integer free parameters and integrality obstructions are different resources and are not automatically bits.

## 1. A real one-bit gain, with a sharp hostile example

[CONSTRAINT_CERTIFICATES.md](CONSTRAINT_CERTIFICATES.md) starts from the inherited atom bank(4,2,5,3,7,13,11). Square-class parity gives16 formal zero-syndrome subsets, but only12 are coprime legal selections. Atoms4 and2 may not be selected together.

Erase the four choices for{4,2,5,13}, and retain that3 is selected while7 and11 are absent. Linear parity alone leaves two completions: all four erased choices absent, or all four present. Coprimality rules out the latter. The missing bit is recovered without another measurement; the remaining source is3, with sigma4.

The same erased positions in another context remain ambiguous between70 and364. Both have square sigma and obey all coprimality guards. Their distinction cannot be inferred from those observations alone. A known F=U+S predicate would select364, but that is additional information.

The recovery has a short Boolean-algebra certificate. If the parity equations imply x0+x1=0 and legality gives x0*x1=0, then

\[
x_0=x_0(x_0+x_1)+x_0x_1\pmod{x_0^2-x_0}
\]

forces x0=0. This motivates a precise code-plus-conflict-complex model. The Boolean coordinate-ring theorem makes such inference expressively complete on finite domains; it does not promise efficient elimination or solve unbounded arithmetic truth.

## 2. The circuit is the ambiguity; the dual check supplies recovery

[ARITHMETIC_FREE_BITS.md](ARITHMETIC_FREE_BITS.md) treats the seven-prime bank

\[
\{2,5,7,19,29,59,89\}.
\]

The parity vectors of p+1 are exactly the seven nonzero vectors of F2³. Square-sigma subsets form Hamming[7,4,3], with enumerator

\[
1+7z^3+7z^4+z^7.
\]

There are **14 minimal circuits**: seven triples and seven four-point complements. The all-seven word is the fifteenth nonzero codeword but is not minimal. A circuit is a smallest hidden toggle; a reconstruction equation belongs to the dual code. Conflating these two roles gives the wrong local test.

The note constructs explicit prime banks for:

- Hamming[15,11,3], realizing every simple binary matroid of rank at most4 by restriction;
- a different pure-square-sigma extended Hamming[16,11,4] bank, with five independent checks and dual RM(1,4).

All98,432 subsets of the three banks are checked by exact arithmetic, not only formal signatures. General prime supply remains conditional; the finite tables require no prime-producing conjecture. For instance, two distinct prime loop columns cannot be supplied: p+1 is square for a prime only at p=3.

### Where7,15 and5 are genuinely related

The seven-coordinate Hamming code has dimension4, hence15 nonzero vectors. Taking those vectors as projective points gives PG(3,2); choosing a basis gives the15 columns of the next Hamming matrix. This is a vector-space construction, not a canonical numerical map between prime products.

The affine lift has16 points v in F2⁴ and check columns(1,v). The constant coordinate plus four linear coordinates explains the five-dimensional dual RM(1,4). A different occurrence of5 comes from a precise recovery theorem below. These counts have declared meanings; they are not identified merely because they recur.

### Error correction versus legal multiplication

One arbitrary bit toggle corrects a Hamming syndrome. It may require removing a selected prime. If only coprime **additions** from the bank are allowed, up to r primes can be necessary for rank r; the bound is sharp exactly when the unused columns form a basis. This is a genuine guard-cost distinction, with exhaustive distributions for r=3,4.

Coset marginals need equally careful typing. Hamming7 cosets agree through order3, Hamming15 through order7. Their separate arithmetic-sign laws also agree. But the joint observation of the sign together with selected bits can distinguish one order earlier: guaranteed joint equality is only through orders2 and6. Separate equal marginals are not joint equality.

## 3. An actual finite arithmetic–Paley correspondence

[FANO_PALEY_INFORMATION.md](FANO_PALEY_INFORMATION.md) recovers the repository's earlier Fano/Hamming/Paley results and makes the arithmetic transport explicit. Fix a labelled Fano plane, a point–line polarity and a baseline Paley7 tournament. Label each Fano line by its nonzero normal a, match that normal to the corresponding prime-bank column, and reverse that triangle when the prime is selected.

For all128 selections,

\[
\sigma(N)\text{ square}\quad\Longleftrightarrow\quad T_N\text{ is Paley7}.
\]

The arithmetic operation is the squarefree product

\[
N\star M=NM/\gcd(N,M)^2,
\]

which corresponds exactly to symmetric difference and commuting named-line flips. This is not ordinary multiplication, aliquot iteration or an unrestricted tournament-composition law. The unit selection maps to the baseline.

The compatible16 masks have the explicit form

\[
c(a)=S+a\cdot v,\qquad S\in\mathbb F_2,\quad v\in\mathbb F_2^3.
\]

The three coordinates of v select one of eight disjoint partner Fano planes. S selects global reversal and is the selection-parity bit, so mu(N)=lambda(N)=(-1)^S. Relative to Fano-preserving GL(3,2), S is the chirality bit; it is not an unrestricted tournament isomorphism invariant. All sixteen Paley tournaments are isomorphic after arbitrary relabelling.

Any three coordinate normals forming a Fano line sum to zero, so their three c-values XOR to S. Geometrically these are three concurrent **primal** Fano lines, not three arbitrary tournament arcs. The exact access rule is:

- no two local bits suffice;
- seven of the35 triples suffice;
- some four-bit sets still fail;
- **every five of the seven bits suffice**.

This is an explicit7/5 recovery theorem. The same bank's fifteen nonzero words give fourteen minimal three/four-triangle flip patterns plus global reversal.

The baseline is a real choice. The Fano automorphism group has two Paley orbits of eight and no fixed Paley tournament, so a fully Fano-symmetric distinguished baseline cannot be selected from the unmarked plane alone. The polarity also matters: a point action and its induced line-normal action are contragredient. The prime-column Frobenius must not silently be identified with the same vertex permutation under an arbitrary dot-product polarity.

The finite census independently checks30 labelled Fano planes, eight disjoint partners per fixed plane, two reverse tournaments per partner union, H=189 and automorphism order21. The quotient189/21=9 counts Hamiltonian-path orbits; it does not identify21 with a forbidden Hamiltonian-path count.

## 4. The old fifteen supersingular primes are a decisive control

[FIFTEEN_PRIME_GRAPHIC_CONTROL.md](FIFTEEN_PRIME_GRAPHIC_CONTROL.md) uses the actual classical supersingular/Ogg list, with its source cited there. Its square-sigma subset matrix also has15 columns, rank4 and2048 kernel words. Nevertheless it is **not permutation-equivalent to Hamming15**: prime3 is a loop, and repeated signatures supply ten weight-two words.

There is an explicit ordinary **five-vertex multigraph** whose labelled-edge cycle space is exactly this arithmetic code. The graph keeps its loop and parallel edges; simplifying it changes the problem.

Six primes from the user's seven-prime bank—2,5,7,19,29,59—form a K4 restriction in that graph. The remaining prime89 supplies the missing class10, completing the Fano plane and destroying graphicness. It is the least prime in that square class, since10k²-1 gives9,39,89 first.

Thus the contrast is precise:

- the classical fifteen-prime square-class model is graphic;
- six of the chosen seven primes form K4;
- adding89 creates the non-graphic Fano restriction.

The missing89 column has an exact higher-order replacement: a four-vertex hyperedge. Its boundary equals each of the three perfect matchings of K4, explaining the three new Fano triples; the other four triples are K4 triangles. These same three matchings index the earlier four-charge energy average, although equal binary boundaries do not imply equal energies.

This is a real bridge between the previous15 and the present7, not a proposed Monster action on graph edges. The old code also has the two-coordinate recovery check x13+x41=0, unlike Hamming15's dual minimum weight8. Equal length and dimension do not imply equal information structure.

## 5. LCA gluing: pay for independent constraints, not named guards

[LINEAR_GLUING_INFORMATION.md](LINEAR_GLUING_INFORMATION.md) recovers the earlier tournament-LCA and correlated-port work, then gives an exact repair formula for affine constraints. If W is the global parity-check space and B ranges over observed scopes, let

\[
W_{\rm local}=\sum_B(W\cap V_B),\qquad
 d=\dim W-\dim W_{\rm local}.
\]

Joining the projected local relations enlarges the true affine set by exactly a factor2^d. The quotient W/W_local identifies which constraints were lost. A basis repairs them. This is a count of missing linear constraints/full syndrome information, not a lower bound for an arbitrary nonlinear membership test.

For covering bags with running intersection, separator subtraction gives an exact rank/entropy formula and an algorithm carrying affine supports plus extension counts. At an LCA a guard is charged only for rank beyond the child constraints; a redundant guard costs zero new binary dimension. These assumptions cannot be discarded for arbitrary nonlinear tournament choices.

For Hamming7, every three-coordinate projection is full. Thus even all seven Fano-line projections discard all three checks. The complementary-four checks recover the code, but the seven named checks have rank3, not7. Their simultaneous success probability under uniform independent input bits is1/8, not1/128.

The earlier K5/Petersen result concerns a different resource: ten pair scopes on five labels need five literal laminar hierarchies. Algebraically dependent parity checks and literal scope covers are not the same optimization problem.

## 6. Going beyond graph cycles without abandoning graphs

[SIGNED_INCIDENCE_RELATIONS.md](SIGNED_INCIDENCE_RELATIONS.md) proves a sharp boundary: existentially hiding auxiliary **edges** in an ordinary graph cycle space still leaves a graphic cycle code. Auxiliary edges alone therefore cannot turn it into the Fano code. A Tanner graph succeeds because its nodes carry equality and parity relations rather than the uniform even-incidence law. The graph plus its local semantics is the actual object.

For integer label reconstruction, signed incidence supplies another kind of information. Balanced signed components have a free integer parameter; unbalanced components give a Z/2 integrality obstruction. Vertex switching removes n-c representational sign bits after a forest gauge choice, but the transformed offsets and vertex reconstruction must travel with it. Positivity and square-value guards are not automatically preserved.

Petersen ear suppression illustrates the danger concretely. Eliminating the middle vertex of a two-edge sum path produces a **difference** equation. The corresponding mixed signed K3,3 retains the source's Z/2 obstruction; replacing it by bare all-sum K3,3 loses that obstruction and invents a free parameter. Hamiltonian lifting and label reconstruction need different information even when they use the same marked ears.

Higher-order incidence goes beyond one parity bit: the k-uniform simplex-facet matrix J-I has cokernel Z/k. Fano point–line incidence has Smith factors(1,1,1,1,2,2,6). That integer matrix is not the3-by7 binary parity-check matrix. Coefficient ring, matrix role and admissibility conditions must all be specified.

## What is established, and what is not

This is a constructive synthesis of **constraint spaces, legal selection complexes, integer relation modules and explicit transport maps**. It supplies recoverable bits, exact costs, no-go examples and algorithms rather than only a warning about information loss.

The general tools are classical. Inherited Fano/Paley, graph-labeling and supersingular results are credited. The explicit new realizations and assembled bridges do not prove prime-producing infinitude, graceful-tree conjectures, general square-sum Hamiltonicity, unrestricted tournament enumeration, or a universal equivalence of arithmetic and geometry.

Berggren's streaming determinant bit, this packet's Paley chirality bit, and an incidence integrality obstruction are distinct objects. Each has an exact mechanism; sharing the word“bit” does not supply a map between them.

## Reproduction

The six published checkers use only the Python standard library:

- `check_constraint_certificates.py`
- `check_lca_gluing.py`
- `check_codes.py`
- `check_signed_incidence.py`
- `check_paley.py`
- `check_fifteen_control.py`

Run each normally and under `python3 -O`. Explicit checks remain enabled. Generated JSON is output, never a required input artifact. Independent audits use separate algorithms for arithmetic subsets, affine recovery, signed Smith forms, graph projections and order-seven tournament/plane censuses. The notes state their finite universes and the proofs supporting all-size conclusions.
