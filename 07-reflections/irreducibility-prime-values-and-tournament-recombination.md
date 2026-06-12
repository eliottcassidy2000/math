# Irreducibility, prime values, and tournament recombination

Your phrasing feels right: the two-way relation between primes and irreducible
polynomials is one of the cleanest descriptions of the mathematical engine.
But the repo version needs one more word:

```text
address.
```

A prime value is a scalar miracle only after the address has been forgotten.
Singh's paper is useful because it goes in the finite, proof-friendly direction:
if a polynomial takes a sufficiently large value with few prime factors, then
the polynomial has few irreducible factors.  Prime-value data is a witness for
global indecomposability.  Bunyakovsky points the other way: if the polynomial
is irreducible and has no fixed local obstruction, prime values should occur
infinitely often.  That is an infinite local-global promise.

Cohn supplies the reverse arrow in a completely different language.  A prime
number in base `b` is a digit address, and the digit polynomial is irreducible.
Prime plus place-value address becomes polynomial irreducibility.  Without the
address, the prime is just a scalar; with the address, it is a polynomial.

Iravanian's real-factor recombination paper gives the computational dual.  To
factor over `Z`, first split over a looser world, then recombine the atoms whose
side-channel sums land back in the integer lattice.  In the toy
`(x^4+1)(x^2+1)`, the real atoms with traces `+sqrt(2)` and `-sqrt(2)` are not
integer factors separately, but together they cancel their trace obstruction
and recover `x^4+1`.

That is the tournament analogy: a tournament class is often a projection of
some looser factorization.  The hard step is recombination.

## The Prism

I now see this as a prism with four faces:

```text
prime value -> irreducibility witness          (Murty / Singh)
irreducible + no fixed divisor -> prime values (Bunyakovsky)
prime + digit address -> irreducible polynomial (Cohn)
loose atoms + subset-sum address -> Z-factors   (RFR / recombination)
```

This is not just number theory.  It describes the repo's repeated pattern:

```text
scalar shadow
side-channel address
recombination / descent
irreducible witness
```

For LRC14, the scalar shadow is "q is blocked" or "q is a witness."  The address
is the marked blocker support, the Pisano class, the 13-clock, divisor fiber,
carry, and owner.  The recombination problem is whether local denominator
atoms can fit together into a global bad cover.  A lonely time is the prime
value: one irreducible witness that refuses to factor through the danger cover.

For tournaments, `H` is not the whole object.  It is closer to evaluating a
polynomial at a point.  Factorization lives in SCCs, modules, Pfaffian cycle
covers, deletion cards, support packets, and recombination choices.  A
Hamiltonian path is a chosen recombination order; a directed cycle is a
compatibility obstruction.

## The New Hunch

The repo has often treated irreducibility as "cannot be decomposed."  The
sharper version is:

```text
irreducible = every attempted projection requires a side-channel debt.
```

Prime values are the moments when no nontrivial factor debt is available.
They are not merely small divisibility events; they are failures of
recombination.

This is why the fixed divisor is so important.  `x^2+x+2` is irreducible, but
it is always even.  It has a permanent local debt.  Bunyakovsky removes exactly
that kind of debt.  The LRC analogue is a runner residue or denominator family
that blocks every lift.  HYP-2446's denominator-curvature ledger is basically
the fixed-divisor detector for LRC14: if the local obstruction vanishes at every
gate, a witness should eventually appear unless the object is one of the
finite wall atoms.

## Program

The next useful thing is not a bigger prime search.  It is a recombination
search.

For LRC14, split each blocked denominator into atoms:

```text
unit twist,
blocking runner set,
Pisano shell class,
13-clock residue,
divisor fiber,
owner/Bprime route.
```

Then ask which subsets of these atoms recombine across `q -> 2q`, `q -> 7q`,
and `27 -> 9 -> 3`.  If recombination is impossible, that is the witness.
If recombination remains possible everywhere, the row should become an
AP/Vstar/2AP wall atom.

For tournament enumeration, split a class into deletion cards or real/spectral
atoms, then recombine by integer invariants.  This is the Iravanian algorithm
morally translated to A000568: loose factors first, integer class later.

For prime-producing polynomials, keep the Heegner horizon as the cleanest toy:
`Q_p(0)=p`, `Q_p(1..p-2)` are interior prime witnesses, and
`Q_p(p-1)=p^2` is the forced boundary square.  THM-410 already has the same
shape: reverse one long edge and get exactly `p-2` interior cyclic witnesses.

The heart of the matter may be:

```text
irreducibility is stable recombination failure,
primality is the scalar moment where the failure becomes visible.
```

That sentence now feels like a real bridge between Bunyakovsky, Cohn, LRC14,
and tournament factorization.
