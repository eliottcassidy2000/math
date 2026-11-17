# Invert the new reader for a source that is not yet known

The [internal reader switch](PROOF.md) normally starts from a physical (h,p). For the opposite direction, fix the desired lower endpoint N and the position of the cut j≥1. Retain the putative high block A and low block r from the reader's address, with q=2^j, ρ=rev_j(r), A>0 even and ρ even. Put G=A+ρ.

The one-arc +2 condition is qp=hG+2. Substitute it into the absolute count formula. With

    L = 2qA² + rG + q(1−U(A)−U(r)),

the exact inverse is

    h = (qN−4qA−2r)/L,
    p = (hG+2)/q.

Thus once the target chooses an address, there is no source-coordinate freedom left. This is not a generic numerical fit: both divisions are exact, and the returned profile still has to occur in a tournament. Reject immediately unless v₂(G)=1, because qp=hG+2, h is odd, and p is even. In that case the primal necessary divisibility gate is the unique odd class

    h ≡ −(G/2)^(-1) mod q.

From p≥2h and ρ≤q−2, a legal candidate has A≥q. Also U(v)≤v−1 for v>0, U(0)=0, so L>0 under A≥q. Positivity and the lower bound on L make the address set finite for each target. This inverse complements the quadratic rounding inverse for a *known* source: it eliminates both unknown counts, exposing the still-open task of choosing a productive address and realizing its exact derived profile.

For the 65-vertex crossover in the packet, N=158860293255395167777, j=25, A=634979302, r=0. The inverse returns exactly h=197,p=3728. Run `python invert_source.py` for this case and further exact controls. It does not call the known-source decoder to obtain h or p.
