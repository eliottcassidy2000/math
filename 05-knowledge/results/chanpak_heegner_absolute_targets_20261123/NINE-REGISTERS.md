# Nine vertices, 36 directed-edge choices, and the five Heegner registers

## 1. Exact equality of register capacity, not of dynamics

For \(p=3,5,11,17,41\) let \(D=4p-1\). This gives exactly the five indicated discriminant magnitudes \(11,19,43,67,163\). The element \(\omega=(1+\sqrt{-D})/2\) in the corresponding imaginary quadratic ring satisfies \(\omega^2-\omega+p=0\), so its field norm is \(p\). In particular \(e^{\pi\sqrt{163}}=e^{\pi\sqrt{4\cdot41-1}}\); the algebraic relation here is the norm/discriminant formula, not a deduction about the near-integrality of that exponential. The quotient register \(\mathbb F_p^\times/\{\pm1\}\) has dimensions \(1,2,5,8,20\) when each quotient position holds one bit. Their total is \(36=\binom92\). Every orientation of \(K_9\) is a bit vector on these 36 coordinates after choosing a labelling and any packing into pairs. The equality of capacities alone does not yield a native modular transition.

**Exact obstruction.** Multiplication by two on the unsigned registers has cycle lengths \(1,2,5,4,4,10,10\). Indeed the first powers of two equal to either sign of one occur at exponents \(1,2,5,4,10\), respectively, and all five return with sign minus. Each multiplier on the quotient group acts by translation, so each orbit at a given prime has that same length. The combined unsigned action has order 20.

An equivariant bijection between these positions and the pairs of \(K_9\), intertwining one clock tick with one fixed vertex permutation, is impossible. The action on pairs of \(K_9\) is faithful, so the vertex permutation would have order 20. The only cycle type on nine vertices with that order is \(5+4\). Its pair orbits have lengths \(5,5,4,2,20\): two distances within the 5-cycle, two within the 4-cycle, and one cross orbit. This disagrees with the register cycle structure. Even an arbitrary bijection of *all* register bit states with *all* orientations cannot globally conjugate the actions. The unsigned clock has fixed states, while the \(5+4\) native permutation fixes no tournament: after two steps it swaps the endpoints of an opposite pair within the 4-cycle.

If an actual arc transports an unsigned bit, choose its reference orientation and transport it around the orbit. This works exactly when the first return to the unordered edge does not interchange endpoints. More generally, for a vertex permutation with cycle lengths \(a_i\), within each cycle the distance-\(d<a_i/2\) pairs give an untwisted orbit of length \(a_i\). The diameters of each even cycle give one **twisted** orbit of length \(a_i/2\). Between different cycles of lengths \(a,b\), there are \(\gcd(a,b)\) untwisted orbits of length \(\operatorname{lcm}(a,b)\). These exhaust all orbits.

The [partition replay](nine-vertex-registers-check.py) gives exact minimum vertex counts for three coordinatewise encodings; unused edge coordinates need not be preserved:

| Coordinates to place on distinct edges | Minimum vertex count | Realizing cycle type |
|---|---:|---|
| Unsigned, on unordered pairs only | 13 | \(5+4+2+2\) |
| Unsigned, with untwisted physical arc bits | 15 | \(5+4+2+2+1+1\) |
| Signed, retaining the choice between each \(x\) and \(-x\) | 72 | \(2+4+8+8+10+20+20\) |

For the first row use length-one orbit from a 2-cycle, length-two diameters of the 4-cycle, length-four orbits within the 4-cycle and between it and one 2-cycle, one within-5 orbit, and two 5–2 cross orbits. For the second use the fixed pair for length one and a 2-cycle–fixed cross orbit for length two, keeping the untwisted length-four, five, and ten choices. The third row follows without search: every register orbit requires order-minus on first return. The only physical return of that type is a diameter orbit, one per even vertex cycle. Therefore the required cycle lengths force all \(72\) vertices, and the diameters attain it. The minima are for a single fixed relabelling and a coordinatewise assignment, not for arbitrary multi-swap gates or general nonlinear encodings.

## 2. A 36-edge grammar from the affine plane of order three

On \(\mathbb F_3^2\), every unordered pair belongs to one affine line. There are four directions, each having three parallel lines, and each line contributes three distinct pairs. Orient the three edges of each line as a directed triangle, independently in either direction. These \(2^{12}\) labelled tournaments are regular, since each of the four lines through a vertex supplies exactly one outneighbor. All have order nine and exhaust this specified grammar, not all regular tournaments. Every line flip preserves scores.

The [exact replay](affine-nine-check.cpp) gives this Hamilton-height histogram (labelled orientations):

| Height | 3159 | 3243 | 3255 | 3275 | 3281 | 3299 | 3303 | 3333 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Multiplicity | 64 | 288 | 576 | 864 | 864 | 864 | 288 | 288 |

In particular no line flip is a unit correction. More strongly, the replay checks all \(4096\cdot36\) individual arc reversals **including exits from the regular grammar**, and none changes the height by \(\pm2\). Hence this especially symmetric 36-edge chart does not even have an outgoing unit arc at width one. This is a finite exact statement about a fixed family.

## 3. Why fixed nine-port charts cannot alone select the exact integer

**Proposition.** For any fixed \(r\), the heights at most \(X\) of all tournaments formed by blowing up vertices of any order-\(\le r\) core into nonempty transitive blocks number \(O_r((1+\log X)^r)\). Every such height has a realization in the same chart on \(O_r(1+\log X)\) vertices.

Indeed ordinal concatenation of strong components multiplies Hamilton counts. A singleton component, replaced by a transitive block of arbitrary width, contributes the factor 1; keep its width one. Every vertex \(v\) of a nontrivial strong component occurs in a directed triangle: the outneighborhood and inneighborhood of \(v\) are nonempty, and some arc from the first to the second exists, since otherwise the outneighborhood could not leave itself. For the width \(m\) at \(v\), take one representative from each of the other two triangle modules. This induced tournament has \(2^m+1\) paths: a single block visit in each of the three positions gives three, and two nonempty visits give \(2^m-2\) ordered partitions of the transitive block. The count of an induced subtournament is at most the ambient count. For a single inserted vertex, every old Hamilton path can be extended by an endpoint or by an adjacent crossing, and deletion of the inserted vertex recovers the original path, giving an injection; iterate. Thus \(2^m+1\le X\). There are finitely many cores and at most \(1+\lfloor\log_2 X\rfloor\) choices per relevant width.

Taking \(r=9\), even the full bank of labelled nine-vertex cores misses a density-one set of odd target heights if the only unbounded operations are these transitive inflations. An exact-integer selector therefore needs recursive composition, nontransitive fillers, or genuinely growing charts. This sparsity says nothing negative about those mechanisms. It separates the capacity identity \(1+2+5+8+20=36\) from the dynamical and absolute-height conditions a constructive solution would actually need.

## 4. Prime-field capacity and a genuine exact lift on order nine

The bit-packing picture in §1 stores one binary marker at each quotient slot. If the slots instead independently carry values of their respective prime fields, the actual joint state space is
\[
 |\mathcal S|=3^1 5^2 11^5 17^8 41^{20}=2^{163.38\ldots}.
\]
A labelling with a fixed directed Hamilton spine leaves \(\binom{n-1}{2}\) freely orientable pairs. Even allowing all anchored orders \(n\le19\) gives fewer than \(2^{154}\) states (the largest exponent is 153, and the other exponents strictly decrease). A lossless encoding of *arbitrary independent joint states* needs order \(\ge20\), or \(\ge21\) if order must be odd. This is only a cardinality lower bound; dynamically reachable subsets might be much smaller.

There is nevertheless a positive, bounded absolute lift. Every order-nine Hamilton count satisfies \(1\le H\le22680=9!/2^4\): pair the first eight positions of a permutation into four disjoint adjacent pairs. In each orbit under swaps of those pairs, just one of sixteen orderings can direct all four within-pair arcs correctly. Since \(3\cdot5\cdot11\cdot17\cdot41=115005\), residues at these five primes separate all actual order-nine heights. Even the product \(5\cdot11\cdot17\cdot41=38335\) suffices with that bound.

Three of the primes in fact suffice. The replay enumerates the \(2^6\) and \(2^{10}\) labelled orders four and five, giving respective maxima 5 and 15. Splitting a nine-vertex Hamilton path after the fourth vertex gives \(H\le\binom94\cdot5\cdot15=9450\), allowing the joining-edge condition to fail in the upper count. Both the height and its desired target are odd, so congruence at \(11,17,41\) forces their difference to be divisible by \(2\cdot11\cdot17\cdot41=15334\). Two numbers in \([1,9450]\) with that difference must coincide. This lifts a modular witness to an exact count **if** one has already constructed an order-nine tournament with those residues. It does not show every in-range residue triple has a tournament, and still less does it lift to arbitrary targets without a proved height envelope.
