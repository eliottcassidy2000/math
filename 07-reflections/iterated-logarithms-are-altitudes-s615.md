# Iterated logarithms are altitudes (S615)

Tao can eyeball whether a bound should carry one log, two, or four. Watching the Collatz epoch experiment, I think I
finally see why it isn't a parlor trick — it's reading the *altitude* of the dynamics.

Picture the logarithm tower over a scale `N`: the value `N`, then its bit-length `log N`, then the bit-length of
*that*, `log log N`, and so on — each floor is "the size of the address" of the floor below. A process that drives
`N` to a constant does so by contracting *some* floor of this tower at a constant rate. The bound it earns is one
floor *above* the floor it contracts:

- Collatz contracts the **value** geometrically (×√3/2 per step). The value lives on floor 0, so the cost — the
  number of steps — lives on floor 1: `log n` steps. One log.
- Coarse-grain into epochs of (bit-length)-many steps. Now the **bit-length** contracts geometrically (×0.79 per
  epoch). Bit-length is floor 1, so the epoch count is floor 2: `log log n`. Two logs. (Measured: R²=0.9987.)

Same orbit, same arithmetic — the log-depth is just *which altimeter you read*. Step-count and epoch-count are the
floor-1 and floor-2 readings of one descent. Tao's four logs in "almost all orbits attain almost bounded values" are
four altimeters stacked: value, epoch, the dyadic range of scales you union-bound over, and the slack of an
arbitrary target function. Each "for almost all" is one floor up. The constant out front — `1/log(1/ρ)` — is just how
fast that floor contracts.

What made it click is that the engine is *altitude-agnostic*. The thing I formalized (`altitude_descent`) knows
nothing about Collatz or floors or logarithms: it's the bare fact that `a_{i+1} ≤ ρ aᵢ + C` forces
`aᵢ ≤ ρⁱ a₀ + C/(1−ρ)`. Feed it the value and you get one log; feed it the bit-length and you get two; feed it a
contraction whose altitude *grows* with N and you get `log*`. The iterated-log is not in the lemma — it's in the
choice of what to feed it. That is the abstraction: **a log is a unit of altitude, and you are free to climb.**

This rhymes with the repo's two faces (HYP-2175). The LRC "almost all configs are very lonely" result is the same
altimeter trick on the additive side: the gap improvement is a union bound over dyadic *frequency* scales, so it,
too, should be loglog-deep. Resonance is what happens when you *can't* climb — when the dynamics refuses to contract
at any finite altitude and the tower never terminates. A cycle is an orbit stuck on floor 0 forever. The conjectures
— Collatz, Lonely Runner — are both the assertion that the only place you get stuck is the ground floor: the integer
1, the lonely time.
