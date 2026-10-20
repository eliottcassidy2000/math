# Independent audit: finite pure-blue atoms and rigid classification

Status: PASS for the deductions from the explicitly cited published inputs, the complete finite certificate check through order 31, and the exact arithmetic cutoff. The unknown odd pure-blue strong alphabet is not claimed classified.

## Source and proof checks

The primary publisher abstract of Benhocine–Wojda, *On the existence of specified cycles in a tournament* (1983), DOI 10.1002/jgt.3190070413, was independently opened. It states the required order-at-least-seven oriented-cycle containment, with its reversal definition. Choosing one reversed edge gives the Hamiltonian bypass used in the proof. No exceptional-digraph theorem is needed.

The free action of Aut(T) on directed Hamiltonian words makes H(T)=|Aut(T)| exactly one path orbit. Camion gives a closing path, so every path in this orbit closes. Consecutive cuts of the same cycle force its full successor rotation to be an automorphism; the order is odd and the tournament is a regular circulant. The published bypass theorem excludes all orders at least seven. The complete order-five multiplier orbit and the explicit nonclosing path (0,2,3,4,1) in connection set {1,2} remove order five. Thus the rigid strong alphabet is exactly singleton and C3, without an order bound.

Consequently all rigid tournament classes are ordered singleton/C3 SCC words. Together with the separately audited even-purity and center laws, every even-order pure-blue count is the stated Narayana-type count. This does not identify odd purity with rigidity.

The cycle-successor injection admits the stated inverse pairing: exactly one of an n-cycle automorphism and its inverse traces a directed Hamiltonian cycle. For pure-blue tournaments this proves 2 HC equals the number of full-cycle automorphisms, at most |Aut|-1.

Cuckler's primary published abstract and Dixon's primary PDF were independently checked for the uniform factorial Hamilton-cycle lower bound and the automorphism-group upper bound. They give qualitative finiteness. Thomassen's 1985 primary indexed publisher abstract was independently checked for the minimum-outdegree factorial bound; its full proof was not independently reconstructed. With this cited input, the exact integer comparisons exclude every odd order at least 133, plus 121,127,129. The quantitative maximum is 131, not a cutoff inferred from an unspecified asymptotic error. The actual exclusions use positive k; no k=0 convention is needed.

The rational class OGF (1+B(z))/(1-z^2-z^6) follows uniquely from ordered SCC decomposition, with B the finite odd polynomial of pure-blue strong classes, including the singleton. Its eventual recurrence is valid for n greater than deg B, and therefore for n>=132 under the explicit bound. This is rationality with an incompletely known numerator, not a completed atom census.

The all-order carousel obstruction is correct: the displayed cycle uses only steps 1,2,m and sends an existing arc to a nonarc under its successor when m>=3.

## Independent finite verification

The standalone verifier `independent_cycle_certificates.py` imports no producer code. It checks:

- Every mask occurs exactly once in every odd order 3 through 31, for 65,534 skew-circulant presentations total.
- All 65,528 negative certificates at orders 7 through 31 have a valid directed Hamiltonian cycle, a genuine cyclic cut, a successor that fails to preserve a specified arc, and a positional reversal that fails to reverse a specified arc.
- A separate bottom-up subset dynamic program reproduces all carousel HC controls through order 15.
- Exact integer arithmetic verifies the cutoff bases at orders 133,135,137 and the three smaller exclusions. The six-order induction multiplies the lower/upper ratio by 4(k+1)/27>1.

`independent_rigidity_cutoff.py` additionally checks every order-five vertex permutation, finding H=15, Aut=5, HC=2 and five nonclosing Hamiltonian paths for regular C5, and verifies that the four order-five skew connection sets form one multiplier orbit.

The audit checks saved witnesses and complete coverage, not the producer's DFS runtime or search-node count. An earlier incomplete guarded run is not used as mathematical evidence. The completed order-31 certificate file supplies all required masks.

Together with the all-order circulant reduction, the complete pure-blue strong alphabet through 31 is C3,C5, with the singleton handled separately. The remaining possible new atom orders are odd orders 33 through 131, excluding 121,127,129. They are not resolved by this packet.

## Reproduction

Python 3 standard library only. Put the main JSON and all extended gzip JSON certificate files together, then run:

    python3 independent_cycle_certificates.py cycle_certificates.json
    python3 independent_rigidity_cutoff.py

The output JSON files record coverage, independently calculated controls, and exact integer comparison values.
