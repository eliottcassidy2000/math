# Put exact three-odd windows on both sides of the target

**Proved:** a deterministic algorithm builds two literal width-three
windows, one below and one above every sufficiently large prescribed
N. Every member is within O(N^(3/8)) of N. No graph search is needed.
It does not imply a window containing N.

Use the [independent-address bank](INVERSE.md) and write

    H_j = P(x,y)+k_j+C(x,y),
    P(x,y) = 432xy+129x+258y,
    C(x,y) = 45(p(x)+p(y))+15(x mod 2)+30(y mod 2),
    k_0,k_1,k_2=135,137,139.

The key is to keep a common low-bit phase, but **charge its large
contribution to the baseline** rather than treat it as error or round
the baseline onto its grid.

## Exact low-bit transfer

Take g=2^s, b between 0 and g−1, and write its s binary digits as
b_0,…,b_(s−1) from left to right. Set A_b=sum_(j:b_j=0) 2^j and let
z_b be 1 if the first digit is zero, 0 otherwise. For s=0 use A=z=0.
Then

    p(gu+b)=p(u)+A_b u+p(b)+z_b(u mod 2).

In the checksum recurrence, appending a zero in position j charges
the current prefix 2^j u plus the previous prefix bits of b.
The parity of that prefix is independent of u except at j=0.
Appending ones charges nothing. This proves the identity.

## Choose a center before opening the grid

Let L=sqrt(N/432), and let g be the largest power of 2 at most L^(1/4).
Choose the smallest positive **integer** c with P(c,c)+135≥N.
Thus c=L+O(1), not L+O(g). Write c=gu+b and choose x,y to be
c+tg,c−tg for a nonnegative integer t. This holds their low digits
equal to b and keeps their high halves' sum at 2u. The large terms
simplify exactly:

    P(x,y)=P(c,c)−D(t), D(t)=432g²t²+129gt.
    p(x)+p(y)=2A_b u+2p(b) + p(u+t)+p(u−t)
                          + z_b[(u+t mod 2)+(u−t mod 2)].

For g≥2 both address parities are fixed, so include 45(c mod 2)
in the baseline. For g=1 handle them as a bounded error instead.
Define B=P(c,c)+90(A_b u+p(b))+the fixed parity. The surplus
Δ=B+135−N is nonnegative and O(L): A_b<g (unless s=0, when
A_b=0), p(b)≤b, and c is only O(1) above L. The remaining
variable correction is nonnegative and bounded by
E=90(u+z_b)+45[g=1], because p(u+t)+p(u−t)≤2u.

Now make two target-dependent choices:

* **Above:** take the largest t with D(t)≤Δ. Even the lowest
  branch is at least N.
* **Below:** take the smallest t with D(t)≥Δ+4+E. Even the
  highest branch is at most N.

For either choice tg=O(sqrt L+g)=o(L), so both addresses and their
high halves are nonnegative, and the addresses are positive at
all sufficiently large targets. The error in the integer square
step is at most one neighboring difference of D. Adding the
variable correction gives

    O(g sqrt L + g² + L/g) = O(L^(3/4)) = O(N^(3/8)).

Algorithmically the construction uses O(log N) comparisons and
integer arithmetic operations on O(log N)-bit numbers, using binary
search. The two explicit tournaments for each center choice use
O(log N) vertices. Both exact triples contain all residues mod 3,
so a member congruent to N can be chosen without worsening the
error. Unlike a generic approximation, the three literal repairs
are attached to the **same computed addresses**.

The low-phase method is independently useful for the
[bare one-apex baseline](bounded/TARGET_PLACEMENT.md). Here the
new feature is the simultaneous one-sided guarantee and a retained
unit repair. The remaining step is substantial: spacing 2 within
a fixed window is not a physical way to correct a general residual
of size O(N^(3/8)).

## Replay

    python verify_inverse.py
    python verify_bracket.py

The second replay checks the exact low-digit transfer and over a
thousand bracketing instances, including 512-bit targets. The first
rebuilds literal graphs and verifies the count formula. Neither
experiment is a substitute for the uniform error proof.
