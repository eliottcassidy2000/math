# The divisor lattice of n is the cover poset of the loneliness problem (S710)

The prompt this session was to take the other model's ("poke") fiber-bundle suggestion further. Last
session (S643) I cashed one seed of it — at a 7-clock the dangerous runners of an `n=14` config are
exactly the multiples of seven — and stopped there, with the mod-2 half and the "56 cells" remark
left hanging. Pushed all the way, the suggestion turns out to be larger and cleaner than a single
divisor: it is the whole **divisor lattice** of `n`, and that lattice is the cover poset of the
loneliness problem.

The fact underneath is one line of CRT and it holds for every divisor at once. Stand at a `d`-clock,
`t = b/d` with `b` coprime to `d`. A runner `v` lands on the observer exactly when `d` divides `v`,
because `vb/d` is an integer iff `d | v`; every other runner sits at distance at least `1/d`, and
since `d` divides `n` that `1/d` is at least the loneliness radius `1/n`. So at a `d`-clock the only
runners that can hurt you are the multiples of `d`, and for a proper divisor they are the *only* ones
not strictly safe. This is THM-421. It says the `d`-clock **peels** the multiples of `d` off the
configuration and leaves them, divided by `d`, as a smaller loneliness problem carried on the `mod
n/d` fiber. The seven-clock of S643 is just the instance `d = 7`; the half-turn is the instance
`d = 2`, where the dangerous runners are the even ones and the safe ones sit maximally far at `1/2`.

What makes it a lattice rather than a single trick is that every divisor gives its own peeling, and
the peelings nest the way the divisors do. For `14` the chain is `1 < 2 < 7 < 14`, and the useful
projection is the *largest* proper divisor, `7 = 14/2`, because it peels the **sparsest** sub-config
— only the one-to-four multiples of seven — and drops the residue onto the smallest fiber, mod two.
The general rule fell out of the computation immediately: best clock `= n / p` for the smallest prime
`p | n`. Even `n` peels off the multiples of `n/2` (about two runners) and lands mod two; an `n`
divisible by three but not two, like fifteen or twenty-one, peels the multiples of `n/3` (about three
runners) and lands mod three. The fiber you end up working in is the smallest prime of `n`. That is
the 2-adic seam again, the same even/odd spine that runs through every part of this project, showing
up here as *which fiber the hard residue lives on*.

The cleanest thing the session produced is not the theorem but the partition it completes. Last
session's THM-420 handled `C'(n)` for the *prime* case through the `2n-1` multiplier shell, reducing
it to the rare transversal core. THM-421 handles the *composite* case through the divisor clock,
reducing it to the sparse residual core. They are complementary halves of one statement: prime `n`
goes through the shell `2n-1`, composite `n` goes through the divisors of `n`. The frontier of the
conjecture — `15, 19, 21, 22, ...`, everything past the proven thirteen — gets sorted by these two
handles, and each handle leaves the same kind of small, finite residual-core question rather than the
whole configuration space. And it explains, sharply, why fourteen is the genuinely hard first case: it
is two-headed. It is composite, so the divisor clock works and tames the seven-side down to a
four-runner core on the mod-2 fiber; but its shell `2n-1 = 27 = 3^3` is a ramified prime power, so the
shell handle slips. Fourteen is the first `n` where the composite handle is fine and the shell handle
fails — the first place the two methods disagree about which one is supposed to be easy.

Then there is the fifty-six. Poke said the half-turn leak "misses only 56 cells," and I had no idea
what the number meant when I read it. It is `A000568(6)` — the number of tournaments on six vertices —
and six is exactly the mover-count of the base, because seven runners means six moving speeds and the
order pattern of their floor-vector is a tournament on those six (THM-381). So the fifty-six is the
**type-space of the base** `LRC(7)`. The fiber bundle `LRC(14) -> LRC(7)` has a fifty-six-cell base,
and the half-turn is the parity coordinate sitting over each of those cells. "The leak misses fifty-
six cells" is the statement that the cover leaves a section over the *entire* base — which is precisely
what S643 verified when every one of fifteen hundred configs came out loose. The leak rides all
fifty-six base types; the genuinely hard ones are the six arc-confined types, `|Arc(7)| = 6` of the
fifty-six, that the project found two years of sessions ago. And the number even carries the 2-adic
seam on its face: `v_2(56) = 3 = (7-1)/2`, the count of quadratic residues mod seven (THM-305). The
base type-space, the fiber coordinate, and the hardness of fourteen are all the same parity fact wearing
three costumes.

I went in thinking I was extending one divisor projection and came out with the divisor lattice as a
cover, the prime/composite partition of the whole frontier, and the resolution of a number I couldn't
place. The other model pointed at a door; behind it was not a trick but the poset.
