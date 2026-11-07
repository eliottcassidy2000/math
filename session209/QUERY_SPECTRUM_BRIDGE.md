# Hamilton counts as exact residual query complexity on paired-black leaves

Status: PROVED corollary of the two independently audited Session209
packets. This is a theorem about a specified marked observation protocol,
not computational complexity for unrestricted tournament reconstruction.

Let T be a labelled tournament representing a fundamental paired-black
quotient leaf. Observe its COMPLETE one-port Hamilton-count response table.
Define R(T) as the minimum worst-case number of further queries needed to
recover the labelled core exactly, over all labelled cores compatible with
that same observed one-port table. Each new query adds fresh ports to
the original core and reports only its Hamiltonian count modulo4.

Then

    R(T) = v3(H(T))
         = log2(size of the complete one-port response fibre of T). (1)

Proof. By the inherited Session14/15 classification, T is either a
nonpalindromic singleton/C3 SCC word or A -> S4 -> A^op for a rigid word A.
If r is its total number of C3 SCCs, H(T) is3^r in the first case and
5*3^r in the second. Thus r=v3(H(T)); in the second case r is even.
The cyclic-port packet's Theorem3 says that the complete one-port fibre,
among ALL labelled tournaments on the same vertex set, is precisely the
2^r independent reversals of those triangles. Each fresh modulo4 response
is binary, since tournament H is odd. Any worst-case decision tree resolving
the fibre therefore needs at least r queries. The r localized two-port
experiments in that theorem attain the bound. This includes r=0. QED.

Consequently the leaf-spectrum packet's normal limits apply without a new
asymptotic argument to R. Under uniform quotient leaves, R/n tends in
probability to0.19425400402459; under uniform labelled tournaments
conditioned on having a leaf image, it tends to0.13357525121498. In either
case its variance is v(a)n+O(1) and its centered, rescaled distribution tends
to a standard Gaussian, with the exact rho,mu,v formulas in that packet.
The statistic R is invariant under choosing a labelled representative, so
it is well-defined also in the first sampling scheme.

This transfers an exact path-count valuation to residual binary-query
complexity and one-port ambiguity, while preserving the sampling measure.
It forgets the locations of the triangle SCCs, which the complete one-port
table supplies. That table may require2^n queries or supplied entries, and
its evaluation and processing costs are not included in R. The theorem
does not make reconstruction linear-time. The fibre members are different
labelled presentations of the same unmarked tournament, not different
quotient nodes. The lower bound does not apply to unrestricted integer
answers, and it does not assert optimal total cost before preprocessing.

At order14, every H27 leaf has R=3 and a fibre of size8. The numeric packet
counts28 such leaves among105 H27 quotient nodes; (1) makes no claim for
the77 nonleaf nodes. The minimum-query identity is exact, while the
distributional constants concern growing order and the stated ensembles.

Dependencies: [leaf spectrum](leaf-spectrum/PROOF.md), especially Sections2,
5 and6; [complete response fibre and optimal decoder](cyclic-port-decoder/PROOF.md),
Theorems2 and3; [numeric H27 stratum](numeric-14-15-27/PROOF.md), Section4.
