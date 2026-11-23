# Sharp continuation memory and physical residue reachability

**PROVED / FINITE-EXACT where marked.** Companion to [reentry depth](REENTRY-DEPTH.md) and the [fixed-order unit toggles](PHYSICAL-TOGGLES.md). Here the physical order can grow, the exterior and two marked exit arcs do not change, and the resonant state is the full modular register, not merely the first two reads.

## Sharp memory at every odd prime

Fix an odd prime \(p\) and order parity. Write \(d=(p-1)/2\), \(v=(Q_1,\ldots,Q_d)\bmod p\). These coordinates suffice for every future uniform read by [the finite modular controller](REENTRY-DEPTH.md#a-finite-modular-controller-and-what-the-two-triples-actually-distinguish). In fact every vector in \(\mathbb F_p^d\) occurs for **strong** fillers of either parity, with independent choices at every prime of a finite set. These guarded moves extend the transitive chain at its marked apex, so they are not uniform grafts of the entire previous filler. See [the explicit apex-word construction](ALL-PRIME-WORDS.md), which does not assume that the two basic moves below act transitively.

For every \(p\ge11\), any prescribed center and gap residues can therefore be synthesized inside the original exterior through \(C=21Q_1+840Q_2\), \(G=2Q_1\), including zero coordinates and any independently prescribed higher register. Choices can be simultaneous at arbitrary finite sets of such primes and at arbitrary powers of those primes. In fact the apex-word construction prescribes all \(d\) displayed reads mod \(p^e\); the remaining derivative data needed for full future prediction modulo a prime power are a separate question. The exceptional first-prime degeneracies at 3, 5 and 7 still come from the physical coefficients, not from unreachable registers.

If two realized vectors first differ at register \(j\), after \(j-1\) successive cyclic sandwiches their ordinary Hamilton residues differ by \((j-1)!\) times this difference. In the sandwich recurrence only the all-up shift can carry register \(j\) to register 1 in \(j-1\) steps, and \((j-1)!\ne0\bmod p\). Thus any state encoding predicting arbitrary future Hamilton residues must distinguish every one of the \(p^d\) vectors. The same lower bound applies to future physical fork exit readouts because their difference is \(4Q_1\), with 4 invertible modulo \(p\). This is a lower bound for *any* encoding, not just linear coordinates. Over a fixed finite odd squarefree modulus, take the product of these counts by simultaneous realizability and the same probe at any prime of disagreement.

## The two shortest physical moves are often enough, but not always

\(P\) adds two singleton ordinal sources to the filler, and \(W\) adjoins the cyclic sandwich \(M\to b\to c\to M\). Both preserve parity and the two chosen exterior reversal edges. For \(\sigma=(-1)^{|M|}\), use boundary coordinates \(v_0=0,\ v_{d+1}=\sigma v_d\). The exact base-prime actions are
\[
(Pv)_s=s^2 v_s,\qquad (Wv)_s=s^2v_s+s(v_{s+1}-v_{s-1}).
\]
On even order, directed breadth-first search from the actual transitive two-vertex filler \((s^2)_{s=1}^d\) gives the following **exact** maximum minimum word lengths for all nonzero states. In particular the \(p=13\) full nonzero orbit fits into at most 54 vertices.

| \(p\) | Reachable nonzero states | Maximum shortest word length | Maximum resulting order |
|---:|---:|---:|---:|
| 3 | 2 | 1 | 4 |
| 5 | 24 | 7 | 16 |
| 7 | 342 | 12 | 26 |
| 11 | 161050 | 21 | 44 |
| 13 | 4826808 | 26 | 54 |

The C++ replay also checks the following test twins. Begin with the transitive tournament on \(p+1\) vertices and its last-to-first arc reversal. Their difference in the register is \(e_1\): the single reversed arc belongs to every cycle, so the odd-cycle contribution to the slot polynomial is \(t[(t+1)^{p-1}-(t-1)^{p-1}]\), which is 1 at \(t=1\) and zero at \(2\le t\le d\) modulo \(p\). The identical words (applied left to right) \(\varnothing\), WWWWPW, WPPWWPP, WPPPWPPPPWPPPWWWPWWP, WPWPWPWWWPPPWWPWWPW at these five primes send the difference to \(e_d\). They are explicit examples where no earlier register sees the hidden difference.

At \(p=7\), the odd sandwich matrix instead satisfies
\[
W^2v=(6v_1+5v_2+2v_3)(1,3,1);
\]
on the surviving line \(W\) acts by 4. This collapse is caused by physical parity, not by a change in the curvature labels. More broadly, continuant determinants can be computed by setting \(a_{-1}=0,a_0=1\) and
\[
a_i=(i^2+\mathbf 1_{i=d}\sigma i)a_{i-1}+i(i-1)a_{i-2}.
\]
This gives \(\det W=a_d\). In the exact scan below 1000, even \(W\) is singular at 23 alone; the odd singular primes are \(3,7,151,307,569\).

Singularity at 23 is a genuine **obstruction to these two moves**. For any fixed seed, words without \(W\) give at most \(d\) vectors because \(P^d=1\). Every word using \(W\) lies in \(\bigcup_{j=0}^{d-1}P^j\operatorname{im} W\) by its last sandwich, a union of at most \(d\) proper hyperplanes, covering at most \(d p^{d-1}\) states. These two sets cannot cover the nonzero space when \(p=23\). This does not restrict the alternate all-prime construction.

Yet \(P,W\) always generate the **full linear matrix algebra**: the \(s^2\) on \(1\le s\le d\) are distinct mod \(p\), so polynomials in \(P\) give all coordinate projections. The off-diagonal entries of \(W\) next to the diagonal are nonzero, so projections before and after \(W\) give each adjacent matrix unit; multiplication gives the rest. In particular there is no common invariant proper nonzero linear subspace, even at 23, but this says only that the *span* of each nonzero orbit is full, not that the actual physical orbit is full.

Replay `c++ -std=c++20 -O3 fork-orbit.cpp -o /tmp/fork-orbit && /tmp/fork-orbit`. All statements concern residues, not equality to a prescribed integer count, and no exterior exit alone is treated as a new center.

[Short-word certificates](SHORT-ORBITS.md) bypass enumeration and extend pure \(P,W\) reachability to all other odd primes through 31 except 23. At 23 a four-vertex detour (singleton, odd sandwich, singleton) instead restores the full nonzero orbit. [Prime-square jets](PRIME-SQUARE-JETS.md) track the extra derivatives and certify full lifts of nonzero base vectors at eight primes.
