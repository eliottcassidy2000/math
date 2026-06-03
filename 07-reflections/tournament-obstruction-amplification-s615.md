# Tournament Obstruction and Arithmetic Amplification, S615

The user's hunch is right in the structural direction, but it needs one careful
separation.

The confirmed common object is a retained carrier whose coarse visible quotient
forgets the important obstruction or amplification channel.  For tournaments,
the carrier is the OCF conflict graph of directed odd cycles, and the striking
unexplained face is the permanent gap `H=7,21`.  For LRC, the carrier is the
`C=2n-1` shell/carry ledger, and the new HYP-2177 result says the doubling
sporadic is controlled by `gcd(n-2,2n-1)=gcd(3,2n-1)`.  For Collatz, the carrier
is the two-block cycle equation around `2^E-3^k`.  For unit distance, the
carrier in the disproof of the grid conjecture is not the visible plane grid;
it is the high-degree arithmetic construction whose projection produces many
unit distances.

The numbers line up suggestively:

- `7 = Phi_3(2)` is the first permanent tournament `H` gap.
- `21 = 3*7` is the second permanent tournament `H` gap and also the LRC
  modulus `2n-1` at `n=11`, where S613 saw a survivor-burden jump.
- `27 = 3^3` is the LRC `n=14` modulus, where the carry-depth problem lives.
- HYP-2177 says the doubling sporadic is controlled by divisibility by `3` in
  `2n-1`.

That is enough to justify a new program: fully understand the low-order
tournament obstruction, then use it as a schema for LRC shell/carry
conservativity and Collatz two-block residuals.

The careful separation is the `1.014` exponent.  Sawin's unit-distance paper
does give more than `n^1.014` unit-distance pairs infinitely often.  But I did
not find a public or in-repo tournament theorem with exponent `1.014`.  The
classical maximum-Hamiltonian-path story is different: `max H(T)` is
Theta-like around `n!/2^(n-1)` with polynomial slack in Alon's upper bound and
regular-tournament constant-factor lower behavior.  Strong tournament minimum
Hamiltonian-path counts have base `5^(1/3)`, also unrelated numerically.

So the right slogan is:

```text
not yet "same 1.014 exponent",
but "same carrier-amplification problem; 1.014 is the unit-distance instance."
```

This is still potent.  It tells us where to look for the tournament analogue:
not in raw `H(T)` over all tournaments, but in constrained carrier tournaments
attached to arithmetic data.  For unit distance, build a tournament whose
vertices are split primes, norm-one generators, projection coordinates, or
proof obligations.  For LRC, build it from shell strata, carry vectors, owner
labels, and proof obligations.  For Collatz, build it from two-block residue
states and Baker-gap obligations.  Then ask whether any of those carrier
tournaments has a thin amplification exponent matching or explaining `0.014`.

S615's route tournament ranked the program accordingly.  The best route is the
side-channel carrier route: OCF forbidden packets, LRC gcd shells, unit-distance
arithmetic towers, and Collatz two-block gaps.  The weakest route is raw
numeric equality without a fixed carrier.

Sources used this session:

- Sawin, "An explicit lower bound for the unit distance problem",
  https://arxiv.org/abs/2605.20579.
- Alon et al., "Remarks on the disproof of the unit distance conjecture",
  https://arxiv.org/abs/2605.20695.
- Alon, "The maximum number of Hamiltonian paths in tournaments",
  https://doi.org/10.1007/BF02128667.
- Busch, "A Note on the Number of Hamiltonian Paths in Strong Tournaments",
  https://doi.org/10.37236/1141.
