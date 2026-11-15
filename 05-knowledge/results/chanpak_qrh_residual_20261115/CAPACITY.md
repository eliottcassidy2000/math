# The complementary lane is sharp in capacity but cannot be the whole selector

Let X be a dyadic target scale, fix g as a power of two comparable to
X^(1/6), and use the three-choice physical window with 432xy as its
dominant term. The complementary placement algorithm with this g
works uniformly for X≤N≤2X. Its addresses obey

    x+y ≡ −1 (mod g),   x+y ≍ X^(1/2),   |x−y| ≤ C X^(1/3).

Here C can be a fixed constant large enough to contain the
construction, independent of X. Restrict to that lane, let all allowed
address pairs and center choices run, and count distinct output
integers (not parameter pairs).

There are O(X^(1/2)/g)=O(X^(1/3)) choices for the sum, and for each
sum O(X^(1/3)) choices for the difference. Thus the lane has at most
O(X^(2/3)) distinct outputs. Conversely the uniform constructive
placement gives, for every integer N in [X,2X], a physical output in
the lane at distance at most C'X^(1/3). Each distinct output can serve
only O(X^(1/3)) prescribed targets. Therefore there are Ω(X^(2/3))
distinct outputs in the lane, in [X−O(X^(1/3)),2X+O(X^(1/3))].
This proves Θ(X^(2/3)) at scale X for the stated lane. The same conclusion
holds if one keeps only one center choice, or keeps the one choice
with the target's mod-3 residue on each branch.

Consequences:
* The order X^(1/3) covering radius is optimal up to constants for
  this fixed-width lane: a substantially smaller radius could not
  reach all X targets with only O(X^(2/3)) outputs.
* Exact modular control for each fixed modulus does not contradict
  this lower bound. A constant-bit suffix reprogramming changes the
  allowed difference width only by a constant, not its exponent; the
  lane still has density zero. A full exact selector must vary another
  parameter on the scale of X, or break the lane's geometric
  restrictions rather than just improve its modular reach.

This explains why separating the resonant main term from the
checksum error is valuable (we obtain the right scale for this
lane), but not sufficient. An analogy with a spectral or averaged
number-theoretic bound does not supply the missing extra states.
