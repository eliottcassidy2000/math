# A literal Collatz-style orbit equivalent to Chan–Pak

This is an effective compilation to a **generalized** Collatz map,
not a reduction to the standard 3x+1 map or an independence result.
The search program is the bounded membership tester in
[Logic and bounds](LOGIC_AND_BOUNDS.md), run on each admissible odd
target in increasing order and halted on the first miss.

## Registers, prime markers, and a divergence clock

Any such finite algorithm can be implemented as a finite counter
program with increment and conditional decrement/zero-test
instructions. For example two stacks encode a Turing tape; stack
push and pop are fixed-base multiplication and division, implemented
with increment/decrement loops and finite remainder control.
Use a distinct state prime p_i for each control state and a
distinct prime r_j for each counter. Double the states into
two layers and alternate layers at every instruction, so no
transition uses the same source and destination prime.

Represent a configuration by
```
n = p_i · product_j r_j^a_j · c^t,
```
where c is one more prime. Increment of counter j and move i→k
becomes the fraction c p_k r_j / p_i. A zero-test/conditional
decrement becomes the ordered pair
```
c p_k / (p_i r_j)     (positive branch, subtract 1)
c p_l / p_i          (zero branch)
```
with the positive branch tested first. Different state primes
prevent all other instructions from firing; the doubling prevents
a canceled state prime from producing a spurious unconditional
fraction. All live steps raise the exponent of c by one.

On entry to the halting state switch to cleanup states. Drain each
counter and then c by ordinary decrement/zero-test macros without
the clock multiplier, and finally divide out the last state prime.
The distinguished starting integer is the start prime (zero
registers). Cleanup reaches 1. No fraction is applicable at 1;
make it fixed. From the distinguished start there are no other
prime factors. More generally we only assert the described
dynamics on this distinguished orbit.

Let D be the least common multiple of the reduced denominators
of this finite fraction list. At each residue modulo D the first
applicable fraction is known from the residue alone. Use that
rational slope on that residue, and slope one if none applies.
This is a total integer-valued piecewise-linear generalized
Collatz map of the kind in
[Conway's FRACTRAN compiler](https://iremi.univ-reunion.fr/wp-content/uploads/2016/10/Fractran_by_Conway3.pdf).

For this explicitly specified compiled orbit:
- Chan–Pak is false iff it reaches 1.
- Chan–Pak is true iff it diverges to +infinity, because every
  step then multiplies the clock valuation by adding one.
- There is no positive-length recurrent nonterminal cycle in
  either case: before cleanup the clock grows strictly; cleanup
  drains bounded counters in a fixed phase order.

This is stronger than just renaming a nonhalting event: it
gives a pure arithmetic map and removes the cycles-versus-
divergence ambiguity for its designated orbit. But the prime
clock is new state, not a consequence of the native Collatz
formula. Nothing in this compiler implies that ordinary Collatz
is universal or that the Chan–Pak instance is independent.

The included replay verifies the instruction/cleanup compiler on
finite counter programs and the elementary bounded checker on
small targets. It is not a listing of the full counter program
for the large bounded checker, nor an observed infinite orbit.

## A stronger analogue: global convergence, with the same truth polarity

There is also a fixed total **piecewise-affine** map G on the positive
integers with
```
Chan–Pak holds  iff  for every n, some iterate G^t(n)=1.
```
This wrapper matters because simply pointing to a nonhalting counterexample
machine gives the opposite polarity and only one start. The wrapper
uses a bounded simulation and can be entered in an arbitrary control
state, just as global Collatz convergence quantifies over arbitrary
inputs.

Use the total program D(x) that checks every admissible target up to x
by the elementary cutoff. It returns BAD if any fails, otherwise OK.
Maintain a protected input register x and a budget register b. An
epoch does the following.

1. Clear every work register. Copy x into the input register of D
   and b into a fuel register, restoring both original values.
2. Simulate D one primitive instruction at a time, charging one
   fuel unit before each instruction. A timeout restarts the epoch
   with b increased by one. An OK output goes to the terminal
   state. A BAD output also restarts with b increased by one.
3. The terminal state goes to the integer 1.

These are finite register macros: clearing and copying use countdown
loops on the register being drained, and no copy loop increments
the register it is draining. In an arbitrary internal state of a
reset macro, its finitely many phases still complete; a first
partial copy may corrupt the input or budget, but any later epoch
entered at the beginning has clean work registers and restores
those then-current values. Similarly an arbitrary simulation
state may begin just after one fuel check, but after at most that
one uncharged instruction it is bounded by its fuel. Thus every
arbitrary control state either terminates or reaches a clean epoch
entry in finitely many steps. From then on the same x is used
and b grows without bound until a definitive OK, if one exists.

To realize this as one integer map with **positive slopes only**,
choose all state and register primes odd and encode live configurations
as 2 times their prime-state word. A valid integer has 2-adic valuation
one, exponent exactly one on precisely one state prime, and exponent
zero on all other state primes. These tests are fixed residues modulo
4 and the product of the squared state primes. Extra odd prime factors
are harmless because the instructions neither inspect nor change them.
Register fractions on valid states preserve the valuation of 2.

Send invalid even integers and the terminal state to n+1. This is
odd. On all odd integers use the two reset branches
```
n -> (n+1)/2 if n=1 mod 4,
n -> (n−1)/2 if n=3 mod 4.
```
Both keep their output odd, strictly decrease n>1, and leave 1
fixed. No reset trajectory can re-enter the even simulation lane.
The combined residue modulus is a common multiple of 4, the state
test modulus, and all instruction denominators. This gives a
total piecewise-affine map with positive rational slopes and no
size tests outside the finite residue rules.

If Chan–Pak is true, D(x) returns OK for every x. Every arbitrary
valid starting configuration reaches a clean epoch or has already
terminated; eventually its budget is large enough for D(x), and it
reaches 1. If Chan–Pak is false, take a clean epoch with x at or
above a counterexample. D(x) returns BAD whenever it finishes,
so that orbit never reaches 1. The ordinary 3x+1 function is not
this wrapper: G uses a new finite modulus and has explicit
control-state and reset branches. The equivalence is a precise
target for comparing proof strategies, not an implication of
undecidability or a change to the usual Collatz conjecture.
