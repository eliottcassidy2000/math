# Put exact three-odd windows on both sides of the target

**Proved:** a deterministic algorithm builds two literal width-three
windows, one below and one above every sufficiently large prescribed
N. Every member is within O(N^(2/5)) of N. No graph search is needed.
It does not imply a window containing N.

Use the [independent-address bank](INVERSE.md) and write

    H_j = P(x,y)+k_j+C(x,y),
    P(x,y) = 432xy+129x+258y,
    C(x,y) = 45(p(x)+p(y))+15(x mod 2)+30(y mod 2),
    k_0,k_1,k_2=135,137,139.

Let L=sqrt(N/432), and g be the largest power of 2 at most L^(1/5).
Insist both x and y be −1 mod g. This is not an externally imposed
congruence on N: it leaves the low log₂(g) bits of each address
all one, which add nothing to their checksums. Then

    p(x)=p(floor(x/g)) ≤ x/g, similarly for y.

Choose the smallest positive c in gZ−1 with P(c,c)+135≥N.
Then c=L+O(g+1) and the overshoot Δ=P(c,c)+135−N is O(Lg).
Set x=c+tg, y=c−tg for t≥0 integer. The large terms simplify exactly:

    P(x,y)=P(c,c)−D(t), D(t)=432g²t²+129gt.

The checksum is nonnegative and at most 45(2c/g+1).
Two target-dependent comparisons suffice:

* **Above:** take the largest t with D(t)≤Δ. Even the lowest
  branch is at least N.
* **Below:** take the smallest t with
  D(t)≥Δ+4+45(2c/g+1). Even the highest branch is at most N.

For both choices tg=O(sqrt(Lg))=o(L), so both addresses are positive
for all sufficiently large targets. The overshoot from integer
rounding is bounded by a single neighboring difference of D,
namely O(g sqrt(Lg)+g²). Adding the checksum bound gives

    O(L^(1/2) g^(3/2) + g² + L/g) = O(L^(4/5)) = O(N^(2/5)).

Algorithmically all thresholds can be tested with integers, using
binary search and the integer upper bound 45(floor(2c/g)+1)
(use that p(x)+p(y) is integral). There are O(log N) comparison/arithmetic
steps on O(log N)-bit integers; the resulting tournaments have
O(log N) vertices. Both physical triples contain all residues mod 3,
so one may choose a member congruent to N without worsening the
error estimate. This construction does not pick from a catalogue;
it programs the two addresses from N using a square-root-like inverse.

The use of trailing ones is the main change from a generic product
approximation: it bounds the actual checksum defect, not just the
size of the linear readout. The two-sided guarantee improves the
previous pair-line coarse placement exponent of about 0.8963 for
the present construction. It leaves a large gap compared with unit
repairs; a complete proof needs a further physical selector that
reduces this prescribed residual without losing the exact gaps.

## Replay

    python verify_inverse.py
    python verify_bracket.py

The latter replays over a thousand bracketing instances, including
512-bit targets; the former verifies the graph-count formula with
literal orientations. These are checks of the algorithm, not evidence
used instead of the error bound.
