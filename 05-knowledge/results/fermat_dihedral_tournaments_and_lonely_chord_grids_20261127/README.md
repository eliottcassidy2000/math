# Dihedral tournaments, chord lengths, and what a Fermat runner grid forgets

This extends the [Fermat polygon chord split](../fermat_polygon_chord_energies_and_period_splits_20261127/README.md), the [Fermat circulant census](../fermat_monster_hamiltonian_counts_prime_classes_20261102/README.md), and the repo's [Farey-prefix runner work](../../../07-reflections/the-cut-prefix-is-a-continued-fraction-times-an-euler-product.md). The finite grid is a tool for constructing witnesses or certifying that the grid itself fails; no statement below settles the unrestricted [Lonely Runner problem](https://arxiv.org/abs/2411.06903) or Chan–Pak.

## 1. A polygon can act dihedrally on paths even when its reflections reverse the tournament

Take any odd \(p\), and choose a set \(S\) containing exactly one of \(v,-v\) for each nonzero residue mod \(p\). Orient \(x\to y\) when \(y-x\in S\). As an ordinary symmetry of the polygon, a translation preserves the tournament and each reflection \(x\mapsto a-x\) takes it to its converse. A reflection cannot preserve all arrows of a nontrivial tournament, since its transposed pair would have to support both directions. On **ordered directed paths**, however, reflecting the labels and then reversing the word preserves directedness. This gives the dihedral action on paths.

If \(H\) is the Hamiltonian-path count and \(F\) is the number of paths fixed by the signed reflection \(x\mapsto-x\) followed by word reversal, the number of orbits of size \(p\) is \(F\) and the number of orbits of size \(2p\) is
\[
(H/p-F)/2.
\]
A rotation fixes no ordered path except for the identity, and a stabilizer contains at most one reflection. In a size-\(p\) orbit, oddness of \(p\) makes all reflections conjugate and each specified reflection fixes exactly one path. Explicitly, a fixed path runs through one representative of each antipodal pair on either side of \(0\); once its half starting at \(0\) is chosen, the other half is forced. This is a signed-transversal path problem, not a path problem on an ordinary unoriented quotient graph.

### At a Fermat prime, this symmetry is exactly dihedral

Let \(p\) now be a Fermat prime. The multiplier action on the orientations \(S\) is free: a nontrivial subgroup of the cyclic 2-group of units contains \(-1\), which cannot preserve \(S\). One can go beyond affine rigidity without a permutation-group classification. The eigenvalues of the circulant adjacency matrix on nonconstant Fourier frequencies \(a\) are \(\sum_{s\in S}\zeta^{as}\). Equality at two nonzero frequencies \(a,b\) would make the difference of their two 0–1 exponent polynomials a multiple of \(1+X+\cdots+X^{p-1}\). The constant coefficient is zero, hence \(aS=bS\), contrary to multiplier freeness unless \(a=b\). The constant eigenvalue \((p-1)/2\) differs from the others by the strict triangle inequality for \(p\ge5\); at 3 the three eigenvalues are visibly distinct. Therefore the spectrum is simple.

Every ordinary automorphism commutes with adjacency, hence with the polygon rotation (both are diagonal in the same simple Fourier eigenbasis). The permutations commuting with a \(p\)-cycle are exactly its powers. Every arrow-reversing symmetry is consequently a rotation followed by reflection. Also an isomorphism between two such labelled circulants must normalize their common translation group, so must be affine. This yields exactly \(2^{(p-1)/2}/(p-1)\) ordinary classes, all self-converse. At \(p=17\) there are 16, and the replay gives their Hamiltonian counts and both types of dihedral path orbit. The total count varies from 13,478,929,757 to 13,689,269,499. Thus the dyadic group controls the symmetry and the path-orbit sizes, not the exact path total or the Chan–Pak spectrum.

## 2. Squared chord length is an exact lonely-runner coordinate, but must be read by a minimum

For nonzero relative speeds \(v_1,\ldots,v_k\), let \(L(t)=\min_i\|v_i t\|\), where the norm is distance to the nearest integer. On the unit circle,
\[
\min_i|1-e^{2\pi i v_it}|^2=4\sin^2(\pi L(t)).
\]
The right side is increasing for \(0\le L(t)\le1/2\). This makes the runner problem an exact shortest-chord problem: at the conjectural time every squared chord must be at least \(4\sin^2(\pi/(k+1))\). It does not turn a sum of squared chords into a minimum.

Consider a prime grid \(t=a/p,\ a\ne0\bmod p\), with all \(v_i\) nonzero mod \(p\). For threshold \(\delta\), let \(B=\{b\ne0:\|b/p\|<\delta\}\), and \(C(a)=\sum_i{\bf1}_B(av_i)\) count bad chords. With quadratic character \(\chi\), the expected bad count on the square/nonsquare half-grid is exactly
\[
\frac{k|B|\ \pm\left(\sum_i\chi(v_i)\right)\left(\sum_{b\in B}\chi(b)\right)}{p-1}.
\]
For each speed, multiplication maps the candidate half-grid to the half-grid with sign \(\chi(v_i)\) times the candidate sign, so the formula counts the size of its intersection with \(B\). If either mean is strictly less than one, that half-grid contains a winning time. At Fermat primes \(p\ge17\) these half-grids are the full \(9/4\) orbits of the earlier note. This is a rigorous selection criterion, not a claim that there must always be sufficient imbalance. Pair overlaps have the equally exact formula \(\sum_a C(a)(C(a)-1)=\sum_{i\ne j}|B\cap (v_i/v_j)B|\); in particular at mean one any overlap forces an uncovered time.

## 3. The sharp masking example: the tournament orientation grid can hide every runner statistic of the chords

Take \(k=(p-1)/2\) and \(\delta=1/(k+1)=2/(p+1)\). Then \(B=\{\pm1\}\). A grid witness exists exactly when the unoriented speed residues do not hit every class \(\{\pm v\}\): the bad class in time \(a/p\) is the single class \(a^{-1}\{\pm1\}\). This includes all cases with a duplicate unoriented speed class.

If the speeds form a transversal, exactly one speed is bad at each nonzero grid time; the loneliness is always \(1/p\). In fact the **entire multiset** of squared chords is the same for every such speed choice and every nonzero grid time:
\[
\{4\sin^2(\pi j/p):1\le j\le(p-1)/2\}.
\]
Its sum is \(p\), because the sum of the cosines over a signed transversal is half the sum over the full nonzero residues. Thus even retaining all symmetric chord moments cannot distinguish these runner schedules on this grid. And exactly these transversals are the data defining the circulant tournaments above. The tournament keeps the orientation of each antipodal pair; measuring chord lengths discards it.

At 17 this single grid signature hides all 16 distinct tournament classes and very different real-line speed lifts. For example, the integer speed sets \(1,\ldots,8\) and \(9,\ldots,16\) are related by reflecting the **residue** orientation, so they yield isomorphic tournaments. Yet their continuous loneliness maxima are respectively **\(1/9\)** and **\(9/25\)**. Reversal of the modular labels is not negation of actual integer velocities once positive representatives are taken. Both cases have exactly the same chord multiset at every nonzero time of denominator 17. This gives a finite, explicit obstruction to deriving the real problem from a single Fermat-prime grid or from an unmarked tournament invariant.

The continuous maxima in the replay are exact rational computations. Each distance \(\|v_it\|\) is piecewise affine with slopes \(\pm v_i\). A maximum of the lower envelope occurs at a break of one distance or an intersection of two affine pieces, so a candidate denominator is at most \(2\max v_i\). The replay tests all these rational times directly. It also computes Hamiltonian counts via subset/endpoint recursion and reflection-fixed paths via the signed halves above.

The useful cross-problem pattern is precise: a Fermat prime makes the multiplier and reflection bookkeeping rigid in group theory, and it makes the square/nonsquare period splitting constructive in polygon geometry. It also supplies explicit finite-grid runner witnesses in an imbalanced regime. But the same symmetry creates a maximally balanced, completely blind grid for a whole class of inputs. Crossing to another denominator or keeping the integer lifts is essential there.

Run `node replay.js` for all 16 representatives, exact runner maxima, dihedral orbit counts, and square-class cover checks.
