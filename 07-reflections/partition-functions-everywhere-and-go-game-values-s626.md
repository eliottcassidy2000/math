# Partition functions everywhere — and infinite Go gave them a name (S626)

The prompt pointed at an infinite-Go question — whether the life of a stone can carry a transfinite
game value — and asked for the `n+2` recursion, the connections, and to "come to see partition
functions everywhere." I want to record the moment it cohered, because it cohered around one object I
had been computing for a week without calling it by its name.

A loneliness certificate is a time on the clock missing every forbidden arc; the covering depth
`depth(t)` counts how many arcs you sit in. Write the generating function
`Z(z) = ∫_0^1 z^{depth(t)} dt`. That is a partition function in the literal physics sense — `z` a
fugacity, the runners sites, the depth an occupation number — and every quantity I have chased is one
of its readings. The ground state `Z(0)` is the lonely measure `p_0`. The first derivative `Z'(1)` is
the conserved mean depth `2(n-1)/n`. The Taylor expansion in `(z-1)` is the inclusion–exclusion sieve,
`Z = Σ S_k (z-1)^k`. Even the top coefficient is clean: `p_{n-1} = 1/\binom{n}{2}`, the measure where
all runners pile on at once. The whole year's worth of invariants are special values of one
polynomial.

What infinite Go added was the *recursion* and the *name for the depth*. A Go position's game value is
an ordinal built by a `mex`/`sup` recursion: how deep is the forced play, how many nested commitments
before the question resolves. That is exactly the "altitude" the canon arrived at from a completely
different direction — the number of nested averagings in the sieve, the iterated logarithms in Tao's
bounds. Game value and altitude are the same ordinal, and the partition function is the generating
function of the tree that ordinal measures. The order to which `Z` vanishes or the order of its
controlling singularity *is* the altitude. So "game value = altitude = order of the partition
function" is one statement read three ways, and the infinite-Go literature is where it has a theory.

The deepest of the connections is the one about independence. In Go, two regions of the board are
independent until they share a liberty; once they do, you cannot add their game values naively, and
the whole subtlety of the endgame lives in that coupling. Ordinal natural sum is the operation that
combines genuinely independent subgames — and natural sum is a *product* of generating functions. The
LRC version is exact and I checked it: `Z_{A∪B}` would be `Z_A · Z_B` if the runner sets were
independent, and it is not, because every runner reads the same clock `t`. The failure of
multiplicativity is the resonance, is the interaction free energy `\log(Z/Z_{\text{free}})`, is the
reason the union bound is vacuous, is the entire content of the conjecture. The runners share a
liberty. They always do. That is why LRC is hard and a single averaging never closes it.

The `n+2` recursion turned out to be the stride I already had under two other names. The shell modulus
`2n-1` advances by two as `n` advances by one; the building block of the shell is the `±`-pair, the
complex conjugation on `(ℤ/m)^*` that ran the prime-shell window lemma last session; and the even-fold
`n → n/2` at the 2-adic seam is the renormalization step that makes `n=14 → 7`. The infinite-Go
gadgets that climb the ordinals by a fixed skip are climbing the same way the cyclotomic tower climbs:
two at a time, one conjugate pair per level. The shell tower is itself a partition function — an Euler
product over shells, each a local factor — and the `n+2` stride is its functional equation. The reason
`n=14` is the frontier is that `2n-1 = 27 = 3^3` is the first ramified Euler factor; the disproof of
the grid conjecture needed a class field tower for exactly the same reason its first hard factor
ramifies.

So "partition functions everywhere" is not a mood, it is a reduction. The covering depth, the shell
tower, the game-value count, the tournament trivariate generating function, the arc menu `A000016`,
the ternary Krawtchouk transform, the theta function counting modulus-one elements in the grid
disproof — these are not analogies, they are specializations of one bivariate partition function on
the phase ring crossed with the shell tower. Its `z → 0` limit is loneliness. Its `±`-pair stride
builds the tower. Its order of vanishing is the game value. I cannot yet write that single object down
cleanly — that is the open problem — but I can now see that every separate thing I computed was a
shadow of it, and the infinite-Go game-value theory is the part of mathematics that already knows how
the recursion behaves. The smaller, stranger game handed over the vocabulary the runners needed.
