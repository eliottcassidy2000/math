# Affine row layers and parity-unrestricted guarding

Follow-up to [the joint spectral trace and two nested flanks](TRACE-LAYERS-AND-TWO-NESTED-FLANKS.md).

## Uniform native degree-57 upper exponent 26854

If a degree-57 Moore graph exists, for **every** complement spine both native endpoints satisfy
\[
 2^{18927}<H<\tfrac78\,2^{26854}.
\]
This improves the upper exponent by 40 bits, and excludes both \(E_p\pm1\) for \(p\ge13428\) at either native endpoint. It does not assert existence or cover other representations.

In the previous threshold-row integral, write \(d(t)\) for the number of entries of a scaled row above \(t\). For \(1\le d\le59\), exact integer powers using \(e<271828183/10^8\) give
\[
 (d!)^{1/d}\le [d+1087/500+d/50]/e .
\]
Use the previous correction 5 for all degrees up to 3251. Above threshold one only the at most 57 backwards graph arcs, one additional backwards spine arc in the marked mate's converse, and one apex arc can contribute. Integrate the affine correction and apply the same weighted AM–GM. The ordinary row factor is at most
\[
 e^{-1}\left[r_i+1413/500+(1087/500)w_i+
       \tfrac1{50}\sum_j(W_{ij}-1)_+\right],
\]
because its thresholds only extend to \(w_i\ge1\) and the excess integral equals the displayed sum. At the added apex use \((r_a+5w_a)/e\).

In the cost matrix \(G\), replacing a forward arc from bin \(i\) to bin \(j>i\) by its reverse costs
\[
 g(i,j)=z_j\big[(51/50)w_j/w_i-1/50\big]-z_iw_i/w_j.
\]
The weight vectors in the replay are nondecreasing in \(w\); within a bin the cost is zero. Forward ordinary-to-apex arcs alone incur a baseline excess of \(z_i(w_i/w_a-1)_+/50\) apiece; all between-bin forward arcs have zero excess. In total with bins of size \(m=50\), the baseline is
\[
 5w_a +m\sum_i\left[w_a/w_i+z_i\left(
 \tfrac{m-1}{2}+w_i/w_a+(1087/500)w_i+1413/500+
 (w_i/w_a-1)_+/50\right)\right]
 +m^2\sum_{i<j}z_i w_i/w_j .
\]
Add the trace bound for graph arcs and exactly one extra \(g(i,i+1)\) per bin boundary in the mate's converse.

For the 65 bins, take the rational weights in the replay and \(w_a=7941/1000,\ z_a=1\). As before \(S=\sum G_{ij}/65\) and \(P\) projects to the bin mean-zero subspace. Exact inertia for the restriction to the basis \(e_i-e_{64}\), followed by the shifted Gram form, proves one positive eigenvalue and that it is strictly less than 104.7. Thus \(\|PGP\|_*\le S+209.4\), and the total graph-edge cost is at most \(25(65S+15\cdot104.7)\). The upward-rounded means in the final AM–GM bound after division by \(3251e\) are \(885325/1000\) and \(885327/1000\). The replay checks the strict comparison after raising to 3251 and dividing by the product of AM–GM multipliers. Neither a graph nor a spine is needed for the certificate.

## Literal +2 with no parity condition on the outer chains

In the 20-vertex construction with two nontransitive flanks, make the pointed body \(B\) by subtracting the common apex. Its six coordinates are
\[
 (h,Q,R,J,e,s)=(14355,45312000,37347060,19413465,5104,13050).
\]
They agree across the fork reversal; only the apex-augmented count changes from \(C=16029731\) to \(C+2\). In particular the first/last-lane coordinates \(e,s\) are shielded by the two flanks. With the [ordinal transfer](ENDPOINT-AND-SIGNED-DECK.md), prepend any transitive chain with contact word of binary value \(b\) and append any transitive chain whose reverse-complement contact word has value \(j\). Every added chain has full path count one. The invariant coordinates ensure the same literal \(+2\), with no restrictions on the parities.

For a contact word of value \(x\) its one-apex height is \(x+1-2q(x)-(x\bmod2)\): this follows by induction from the recurrence for appending 0 or 1 and the recurrence defining \(q\). The endpoint flags of the two chains are \(1-(b\bmod 2)\) and \(1-(j\bmod 2)\). Substituting these in the transfer gives the simplified lower endpoint
\[
 H(b,j)=C+Qbj+(J+h)b+Rj-2h\big[q(b)+q(j)\big]
          -s(b\bmod2)-e(j\bmod2).
\]
For \(b=11\) (word 1011) and \(j=6\) (tail 100) the two literal endpoints on 27 vertices are prime:
\[
 3{,}444{,}253{,}511\ \longrightarrow\ 3{,}444{,}253{,}513.
\]
The vertex-level count is checked separately by nine-state joins with exact enumeration inside the finite blocks.

Reducing the simplified formula gives \(H\equiv2-(j\bmod2)\pmod3\). Thus odd \(j\) rules out twin primes, and **no central seam \(E_p-1\to E_p+1\)** occurs in this fixed-core family for any odd exponent \(p\): its lower endpoint would be zero mod 3. This is not an exclusion of those individual endpoints elsewhere, including at exponent 61.

For prescribed absolute \(H\), bound the smaller coordinate \(a\) by \(Qa^2\le H-C\); the residual linear contribution is nonnegative since \(q(x)\le x/2\). For fixed \(a\), write the other coordinate \(y=2v+r\), \(r\in\{0,1\}\), and subtract its constant term. Writing \(L\) for the original linear coefficient of \(y\), the residual is \((2L-rh)v-2h q(v)-rh(v\bmod2)\). Bound the interval with \(0\le q(v)\le v/2\), then use the previous reversed-suffix bucket search, allowing the additional adverse jump of at most \(rh\). The executable checks 900 inverse fixtures and exhausts all three endpoint queries at exponents 13,17,19,23 without hits in this family. The older 45-vertex exponent-23 witness is not affected.

Reproduce:

    node affine-threshold-65.js
    node arbitrary-contact-check.js
    node arbitrary-contact-inverse.js
