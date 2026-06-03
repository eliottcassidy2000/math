# The window lemma ports, and stops exactly at the ramified shell (S625)

The instruction was to take the grid-disproof's machinery — modulus-one elements in CM fields, class
field towers, bounded root discriminant, the geometry-of-numbers window — and port it into a uniform
LRC bound, turning the constructive route `M ≥ 2/(2n-1)` from a per-`n` fact into an attack. I did the
port. It works cleanly, it gives a real theorem, and it stops at a wall that is itself the most
interesting thing I found.

The port is the **prime-shell window lemma** (THM-418). Evaluate the runners at a prime shell `p`:
the multipliers are `(ℤ/p)^* = Gal(ℚ(ζ_p)/ℚ)`, the rotations of the clock. Each runner can be pushed
into the width-two danger band by exactly two of them (`a = ±v_i^{-1}`), so the band blocks at most
`2(n-1)` rotations out of `p-1`. The moment `p-1 > 2(n-1)`, i.e. `p ≥ 2n`, a rotation survives and
every runner lands at distance `≥ 2`, giving `M ≥ 2/p`. That single inequality — *more rotations than
the band can block* — is the disproof's window lemma with all the analysis stripped away: in the plane
you need a high-degree field to manufacture enough modulus-one vectors; on the clock the prime shell
hands you `p-1` of them for free. The proof is three lines.

And three lines is exactly the problem. `p ≥ 2n` forces `2/p < 2/(2n-1)`; the bound is `≈ 1/(2n)`,
half of the conjectured `1/n`, trivial-strength. The factor of two is not slack in the argument — it
is a wall. To reach the optimum you must use the shell `2n-1` itself, and there the count is
critically tight: `p-1 = 2(n-1)`, the window is exactly full, and counting guarantees nothing. The
extremal configs *do* have a surviving rotation (THM-415 is attained), but you cannot see it by
counting. You have to look at the arithmetic of the shell.

Which arithmetic depends on whether `2n-1` is prime. When it is prime, the shell is "unramified" — one
clean cyclotomic field, the window critically full, and the last rotation is there but invisible to
the pigeonhole. When `2n-1` is a prime power — `27 = 3³` at `n = 14`, the first even frontier — the
shell *ramifies*. The non-unit runners (multiples of three) fall out of the band automatically, the
real structure lives in the tower of subfields `ℚ(ζ_3) ⊂ ℚ(ζ_9) ⊂ ℚ(ζ_27)`, and the surviving
rotations have to be hunted level by level. I checked: for the AP-minus-two extremal at `n = 14`,
exactly two rotations survive on the shell 27, `13` and `14`. Two, out of eighteen units. The window
is nearly closed, and what keeps it open is the three-adic tower, not a count.

So the port reproduces, in miniature and on the nose, the entire shape of the grid disproof. A single
unramified prime gets you a constant factor and no further — that is the easy half of the disproof and
it is THM-418. The true exponent requires an infinite class field tower with bounded root discriminant
— and on the LRC side the true value `2/(2n-1)` requires the cyclotomic tower of the ramified shell,
with a uniform control on how the runners can concentrate across its levels playing the role of the
bounded discriminant. The grid disproof's `ε ≈ 6·10^{-38}` and our missing factor of two are the same
gap: the part that counting cannot reach and only the tower can.

I did not close it. But I now know precisely what closing it is: a tower-descent lemma saying `n-1`
units cannot cover every `±`-class up and down the cyclotomic tower of `2n-1` at once, with the
covering controlled uniformly in `n`. That is HYP-2240, and it is the same theorem the disproof proves
in the plane, asked on the clock. THM-407 told us a year of sessions ago that `n = 14` is hard because
`2n-1` is the first even case where doubling stops being shell-transitive — that *is* ramification, and
I had been staring at a class field tower without the name. The unit-distance disproof gave it the
name.
