# The global arithmetic wrapper has no nonterminal cycles

The [positive-slope generalized Collatz wrapper](FRACTRAN_ORBIT.md)
can be chosen so that every orbit either reaches 1 or diverges to
infinity, independently of the truth of Chan–Pak. Thus the equivalence
with global convergence puts the uncertainty entirely in divergence,
not in the possible existence of hidden finite cycles.

In a clean epoch, keep its budget b protected except while copying
it to the fuel register. Use a private scratch counter. A copy
iteration decrements the budget, increments scratch, then increments
fuel. Restore the budget by moving scratch back one unit at a time.
At every primitive step of these two macros, the sum of the budget
and scratch is at least the starting budget minus one. Other reset
phases and the simulation do not alter the budget; a restart
increases it by one. Scratch is cleared before each clean copy.
These requirements are compatible with the arbitrary-state recovery
argument: only the first partial epoch can enter one of these macros
with corrupted scratch.

Let p_min>1 be the smaller of the two primes encoding budget and
scratch. During a clean epoch of starting budget b, its even encoded
integer is at least p_min^(b−1) throughout. Every nonterminating orbit
reaches a clean epoch and starts infinitely many further epochs, each
of finite length and with growing budget. Its integers therefore
tend to infinity. Invalid and terminal states enter the odd reset
lane, which strictly decreases to 1. In particular no other cycles
are possible in this compiled map, including from arbitrary starts.
This is not a statement about the unmodified 3x+1 map: the protected
budget and prime encoding are extra coordinates.
