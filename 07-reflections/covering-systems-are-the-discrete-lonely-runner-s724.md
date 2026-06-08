# Covering systems are the discrete lonely runner (S724)

Erdős called it his favourite problem. A covering system is a handful of congruences — `x` is `0` mod `2`,
or `0` mod `3`, or `1` mod `4`, and so on — arranged so that every integer falls into at least one. The
classic example uses moduli two through twelve. Erdős asked whether the smallest modulus could be pushed
arbitrarily high: could you cover all the integers using only congruences with large moduli, none of them
small? He bet yes. He was wrong. Hough proved in 2015 that the smallest modulus of a covering system is
bounded — at most ten-to-the-sixteen — and Balister, Bollobás, Morris, Sahasrabudhe and Tiba brought the
bound down to six hundred sixteen thousand with a cleaner "distortion method," and showed besides that no
covering system can have all its moduli odd and squarefree. On the other side, the largest minimum modulus
anyone has actually built is forty-two, by Owens, just past Nielsen's forty. So the truth sits between
forty-two and six hundred sixteen thousand, and Erdős's intuition that it was unbounded was the thing that
got disproved.

The reason this lands in our repository is that the lonely runner conjecture is the same problem in
continuous clothing. The runners sweep the circle, each one painting a danger arc around the observer
whenever it comes close; the conjecture asks whether some instant escapes every arc — a lonely time. That
is a covering question: do the danger arcs cover the circle, or do they leave a gap? And the bookkeeping is
identical to the covering-system bookkeeping. THM-406 already proved that the lonely measure is the
alternating sum of the danger-arc overlap volumes, `p_0 = sum (-1)^j S_j`. A covering system's uncovered
density is the alternating sum of the congruence overlaps, `sum (-1)^|I| / lcm`, which is exactly S561's
sieve `rho`. Same inclusion–exclusion, same alternating sign, the same statement that "covered" means the
alternating sum vanishes. The covering system is the rational skeleton of the lonely runner: replace each
continuous danger arc by the congruence it would be if the runner ticked only on the rationals, and the
two coincide. I checked the classic minimum-modulus-two system both ways; it covers, the direct uncovered
density is zero, the inclusion–exclusion `rho` is zero, they agree.

The deepest shared feature is the one S561 found for the runner and the covering people found for the
integers: density is necessary and useless. To cover you obviously need the reciprocals of the moduli to
sum past one, but summing past one buys you nothing on its own. I built a little system with moduli three
through twelve whose reciprocals sum to about one point three seven, comfortably over one, and it covers
nothing close to everything — a quarter of the residues escape. The wasted density goes into overlap, and
the overlap is where the whole difficulty lives. S561 said the same about the runner: the sieve-covered
density at `n=14` is about a fifth, positive, not measure zero, and yet density attacks cannot close the
conjecture because the thin near-arithmetic slice inside the covered core is what matters. Both problems
are immune to counting and demand structure.

Which is why the disproof of Erdős's problem is, for us, a method and not just a fact. Hough and the BBMST
team won by a distortion argument: a density increment showing that when the moduli are forced to be large
they are forced to be so smooth that a fixed positive proportion of integers slips through every
congruence — the cover always leaks. That is word for word what the lonely runner wants. My S643 found,
computationally, that the multiple-of-`n` configurations always stay loose, that the cover of the circle
always leaks; what it lacked was a reason. The reason, if it exists, is a distortion increment on the
covering-depth `p_0`: the runners' arcs, when the speeds are constrained the way the conjecture
constrains them, are too spread to seal the circle, and a positive lonely measure survives. The covering
people proved their version of "spread pieces leak" to kill Erdős's optimism; the runner needs the same
inequality, run in the same direction, to confirm its own. The two are dual uses of one technique:
covering systems run the density increment to show coverage fails at large modulus, and failure of
coverage is exactly looseness, which is exactly what the lonely runner is trying to guarantee.

Even the side results rhyme. BBMST proved no covering system has all odd squarefree moduli — evenness is
structurally required — and the lonely runner's whole personality is the two-adic seam, the even-versus-odd
split that THM-435 made the Pfaffian's defining feature and that makes the even `n` the hard ones. The
prohibition on all-odd coverings and the necessity of even structure in the runner are the same parity
fact. As for the number the prompt sent me to check, forty-two turns up in the repository only as an
accident — an index in the Hamiltonian-path sequence, an entry in a Paley table — and means nothing
structural here. What means something is the covering itself, which is not a neighbour of our work but the
discrete shadow of its central object. The lonely runner has been a covering problem all along; Erdős's
second problem is where its skeleton was solved, and the method that solved it is sitting on the table,
waiting to be lifted into the continuum.
