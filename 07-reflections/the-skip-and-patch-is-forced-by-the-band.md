# The skip-and-patch is forced by the band

*mac-mini-2026-07-01-S77. Reflection on HYP-3792.*

For many sessions the construction `{1, ..., n-2, n(n-1)}` was described as a *choice*: skip the runner
`n-1`, patch with the multiple `n(n-1)`, and — as if by luck — the numbers `Phi6 = n^2-n+1`, the Dedekind
margin, the continued fraction `[0; n-1, n]`, the units `(Z/14)^*` all fall out. This session the choice
dissolved into a *forcing*, and the forcing lives in a single picture: the residues on `Z/Phi6` at the
binding time.

Write the loneliness condition at a rational time `a/q` as a statement about residues: the runner `v` is
safe iff `a v mod q` avoids the band `(-rq, rq)` around zero. Then `M(S)` is just the deepest band any
dilation of `S` can dodge, `max_{q,a} (1/q) min_v ||av||_q`. At the construction's own witness `t* =
n/Phi6` the picture is rigid and beautiful. The core `1, ..., n-2` lands on the arithmetic progression `n,
2n, ..., (n-2)n` — evenly tiling the safe band. And because `n(n-1) \equiv -1 \pmod{Phi6}`, the multiples
of `n-1` land at `-1, -2, -3, ...`: the runner `k(n-1)` sits at distance exactly `k` from zero. So the
multiples of `n-1` are precisely the *dangerous* runners, marching outward from the origin one step at a
time.

Now the covering hypothesis speaks. A covering set *must* contain a multiple of `n-1` — that is the `q =
n-1` obligation, non-negotiable. At `t*` that forced multiple sits at distance equal to its index. To be
safe it must have index `>= n`, i.e. it must be `>= n(n-1)`. The smallest legal choice is `n(n-1)` itself.
And the runner `n-1` — index 1, sitting at residue `-1`, hard against zero — is the one runner that can
*never* be safe here, so it is exactly the one to drop. The skip and the patch are the unique answer to a
forced question: *cover `n-1` without standing next to zero at `t*`.* There was no cleverness. There was a
band, and an obligation, and only one way to satisfy both.

Everything else is a shadow of this. `Phi6` is the denominator that makes `n(n-1) \equiv -1`. The margin is
the Dedekind sum because `-1` is an order-6 (hexagonal) unit and not an order-4 (square) one. The continued
fraction `[0; n-1, n]` records the two moves — dilate by the skipped level `n-1`, correct by the patch
level `n`. The units `(Z/14)^*` are the six band edges of the *tight* cousin (the AP), the same residue
picture at `q = n` instead of `q = Phi6`. The zoo was never a zoo; it was one arithmetic diagram seen from
six sides.

The part that turns characterization toward proof is the isolation. Sample covering 13-sets at random and
their loneliness clusters around `0.14`, never dipping below `0.108` — a clean factor above the construction
at `0.0765` and a clean factor and a half above the conjecture's own threshold `1/14`. To get anywhere near
the threshold you must reach the *deep* modulus `q = Phi6`, and only the `n(n-1)`-patch reaches it; a
random covering set binds *shallow*, at denominators in the teens, comfortably far from danger. So the
Lonely Runner Conjecture, in its covering case, is not a claim about a knife-edge that most sets flirt with.
It is a claim about one arithmetically special well, isolated in a landscape where everything else is safe
by a wide margin. And that margin is spendable: the theorem asks for `1/14`, the construction sits at
`14/183`, and the bulk sits at `0.108` — so a certificate never has to be sharp, only to clear `1/14`, with
the Dedekind margin `13/2562` of room to spare.

The pattern that transcends the theorem: **when an extremal object looks like a lucky coincidence of many
invariants, look for the one obligation that forces it, and write the problem in the coordinates where that
obligation is visible.** Here the coordinates are residues on `Z/Phi6`, the obligation is "cover `n-1`," and
the coincidence becomes a deduction. The invariants did not conspire; they were bookkeeping for a single
forced move, and the bookkeeping only looked deep because we were reading it in the wrong chart.
