# Coverage is not monotone; deletion-robust coverage does pad

23 October2026. **REFUTED:** full k-subset charge coverage alone implies
full (k+1)-subset coverage whenever k<h/2. **PROVED:** the forced-padding
equivalence below.

Let G=(Z/2)^2 and take the four distinct labels 00,01,10,11, each once.
Singleton charges cover G. Every pair of distinct labels has nonzero sum,
and all three nonzero sums occur, so pair charges omit00. Here k=1<4/2.
The missing coordinate is the requirement that the added slot be distinct
from every slot already used. Unrestricted repetition would erase it.

This has the smallest possible number h of supplied labels for positive k.
For h<=2 there is no positive integer k<h/2. For h=3 the only possible k
is1; singleton coverage forces |G|<=3. The trivial group is immediate; in
Z2 two label types among three labels give both pair sums; in Z3 covering
all charges requires labels0,1,2 and their pairs cover all charges. The
checker exhausts all multisets of length<=4 in every abelian group of
order<=4, finding precisely the displayed hostile. No larger group can
give a smaller-h singleton-coverage counterexample.

## The exact K4 structure

Make the four label slots the vertices of K4, named0,1,2,3 by their two-bit
values. Colour an unordered distinct-slot edge{a,b} by a+b in G, equivalently
their bitwise XOR. Its three nonzero colours are the perfect matchings

    colour1: {0,1}, {2,3};
    colour2: {0,2}, {1,3};
    colour3: {0,3}, {1,2}.

The missing zero colour would lie on the diagonal{a,a}: it reuses the same
resource, which is not an edge of this simple graph. Every charge occurs
at occupancy1 and also at occupancy3, since the four-label sum is zero
and the triple omitting a has sum a. Occupancy2 loses exactly zero.

This is an exact map from distinct-slot pairs to pair charges, not an
arbitrarily oriented tournament. It retains the disjoint perfect-matching
partition and identifies the failed padding step as resource reuse. If
another genuine slot with the same charge were supplied, its pair with
the first slot would be legal: equal charge is not itself equal resource.

## The stronger property needed by the native bank

For an indexed label bank B=(ell_1,...,ell_h), let S_k(B) be the set of
charges of k distinct slots. Fix a set F of r forced slots. Then exactly

    {charges of (k+r)-subsets containing F}
       = sum_(i in F) ell_i + S_k(B\F).

Thus full coverage by (k+r)-subsets **containing those exact slots** is
equivalent to full k-subset coverage after deleting F. This is a translation
identity in any abelian group, and includes multisets because slots retain
their distinct indices.

In particular, if every one-slot deletion of the53-source bank retains
complete three-subset charge coverage, then for **each prescribed source**
there is complete four-subset coverage containing that source. To solve
charge q, subtract its label, solve the residual using three distinct slots
in the deleted bank, and adjoin the prescribed slot. This is stronger than
unforced four-subset coverage. The minimal counterexample above does not
have this deletion-robust hypothesis.

The identity supplies only the charge step. Actual path choices, shared
pins, resource compatibility, target margins and physical per-owner guards
must still be authenticated by the native bank theorem. No abstract padding
argument creates a missing helper record or proves bounded connectivity.

Run:

    python 05-knowledge/results/compatible_charged_helpers_20261023/controller/padding_probe.py
    python -O 05-knowledge/results/compatible_charged_helpers_20261023/controller/padding_probe.py

The finite control enumerates the minimum-hostile universe and checks the
forced-one-slot identity independently against direct subset enumeration.
