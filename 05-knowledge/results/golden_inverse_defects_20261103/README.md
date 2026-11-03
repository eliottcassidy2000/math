# Realizable inverse states, ternary offsets and quadratic transport

**PROVED scoped constructions and inverse / FINITE-EXACT certificates.
Same-agent proof review with independent computational paths. The global
claim that7 and21 are the only missing odd Hamilton counts remains OPEN.**

The inherited65-presentation bank now has a terminating inverse that
returns an actual source word or proves the proposed full state absent
at every source order. It retains all five Q4 presentations. This is a
membership algorithm for a specified language, not a proof that every
allowed Hamilton count belongs to that language.

The main positive construction is an infinite family with both ternary
offsets. In the inherited source alphabet set

    U  = a a r5,               v(U) =(31,14,25),
    U2 = a a q1 a a,           v(U2)=(95,41,75),
    U4 = a a a a a q0,         v(U4)=(97,40,75).

For every prefix w over triangle b and Q4 phase q2,

    H(w U2)=3H(w U)+2(-1)^(number of b symbols),
    H(w U4)=3H(w U)+4(-1)^(number of b symbols).

All targets are actual strong tournaments. Even triangle parity gives
both requested positive offsets. These operations replace a bounded
suffix; they do not preserve every old edge and are not claimed to be
single complement/re-marking edges. All other Q4 phases remain available
in the full inverse, with their different correction transports retained.

There are also decisive obstructions, rather than a proposed universal
attachment rule hidden inside the construction:

- For every bank state h>=3, all left word attachments reaching3h+2 or
  3h+4 are exactly one or two singletons with the guards in
  [THEORY I4](THEORY.md). No state admits both offsets by left attachments.
- All eleven source words with H25 fail to reach77 or79 by any two-sided
  wrapper. Singleton sources a^n, n>=6, fail under arbitrary insertions.
- The95 and97 output states fail the next four proposed affine updates
  at every order in the bank, although those four scalar counts occur
  elsewhere. A closed induction on realizable states is still missing.

## What the coordinates now mean

For v=(h,x,y), the positive coordinates

    (d,e,f)=(x+y-h,y-x,h+x-y)

count a single backward-step defect at the boundary port, paths ending
there, and paths ending elsewhere. Their sum strictly increases under
source extension. Inverting the integer matrices with strict descent
therefore supplies a realization guard and a finite generating polynomial
for every full-state fibre. The branch d=1 or2 has a complete all-size
classification:66 source words and41 states, maximum H155.

Incoming [universal block transfer](../golden_hamilton_transfer_20261103/RECURRENCE.md)
has already settled arbitrary-block three-coordinate compilation; it is
credited, not claimed as a result of this packet. We use its path-cover
interpretation to connect the defect to the incoming
[one-backward-step ear term](../session153_gadgets_20261103/cuts/README.md).

## Exact arithmetic transport

[TRANSPORT T3--T8](TRANSPORT.md) proves rational generating recurrences,
the quadratic clocks9+/-2sqrt(19) and4+/-sqrt(17), and their ramified
orders171 modulo19 and68 modulo17. Modulo7, q2 is literally multiplication
by a primitive element of F49. Frobenius has the phased rule
k->7k+36 modulo48. Each positive offset has exactly seven divisibility
positions in this48-step clock. Quotient primality remains an additional
condition.

Allowing all Q4 phases expands the correction module from a sign character
to the full lattice{(2u,x,y)}. Its exact minimum observer modulo m for
arbitrary integer corrections under left contexts is
(Z/(m/gcd(m,2)))^3. This minimum is not a statement about arbitrary joint
ports or full metagraph closure.

The positive prefix family supplies seven prime outputs and1253=7*179
below one million. It has at most O(X^(log_5 2)) outputs up to X and is
therefore too sparse to complete the odd spectrum by itself. No global
coverage ceiling or A000568 coefficient is extended here.

[RECOVERY.md](RECOVERY.md) identifies the requested ten-number list as
the already known H9 parent's complete nonconstant ear image, computes
the new source-fibre polynomials, and retains the existing correction of
the pasted Ogg/Fermat proposal. Equal numerals are not structural maps.

Reproduce from the repository root:

```powershell
python -B -X utf8 05-knowledge/results/golden_inverse_defects_20261103/verify.py
python -B -O -X utf8 05-knowledge/results/golden_inverse_defects_20261103/verify.py
python -B -X utf8 05-knowledge/results/golden_inverse_defects_20261103/seal.py
```

See [AUDIT.md](AUDIT.md) for exact universes and proof boundaries. The
next constructive target is an invariant family of actual source states
closed under all three ternary branches, or a different realization
theorem for the two prime lanes. The present inverse can certify or
reject proposed states in such a proof without silently resetting them.
