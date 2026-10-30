# Next reroute: retain precision where the current content can see it

**PROVED:** a one-sided nondeterministic simulation with O(E*5^E) states and
an exact refinement-credit identity. **VERIFIED:** all-word abstract potential
certificates at E=1..5, with all positive letters covered. **OPEN:** a uniform
potential for this simulation, useful higher-layer terminal fees, and all-layer
content control. These lower-precision checks do not improve the inherited
E11 theorem; they test a different scalable proof interface.

## Exact amount earned by delaying one further layer

Use the inherited epoch enclosure F_r. If h<=r it is already exact. If h>r,
write d_(r+1)=d_r+j*T_r, j in{0,1,2,3,4}, and T_(r+1)=5*T_r. The low-layer
counts through r agree in both enclosures. Subtracting their formulas gives:

* L<d_r: the improvement F_r-F_(r+1) is0.
* d_r<=L<d_(r+1): it is
  (h-r)+(lambda_h-lambda_r)(L-d_r).
* L>=d_(r+1): it is

      {(L-d_(r+1))/T_(r+1)} + j/5
        + (lambda_h-lambda_(r+1))*j*T_r.                 (A1)

Here braces denote fractional part. To check the last line, use
lambda_(r+1)-lambda_r=1/T_(r+1) and
n_(r+1)=1+floor((L-d_(r+1))/T_(r+1)). The constant h-r decreases by1
while the new exact count supplies that1; the remaining terms are exactly(A1).
Every term is nonnegative, proving F_(r+1)<=F_r for every epoch, with equality
and strictness explicitly determined. This identifies the actual gain from
using a digit: a later first return and a rounding credit, not the existence
of a larger state space by itself.
[delay_credit.py](delay_credit.py) checks7,863 exact boundary and periodic
cases against the separate closed-epoch formula; the subtraction above proves
the identity for all epochs.

## An adaptive abstract observer

Fix E>=1, q=5^E. Retain z modulo q exactly, set h=min(v5(z),E), and retain
the multiplier m only modulo5^ell with ell=max(1,h). A state stores a unit
residue m. Its current content is determined exactly:

    reward(z,m)=min(v5(m-1),h).

For a letter a, first compute the exact new z by

    z' = (3z+2^a-2)/2^a modulo q.

Let ell'=max(1,min(v5(z'),E)). If ell'<=ell, the new multiplier is determined
by3m/2^a modulo5^ell'. If ell'>ell, allow **all** lifts of3m/2^a from
modulo5^ell to modulo5^ell'. This is a branching upper enclosure, not a
deterministic quotient of actual history.

Every actual affine word has a path in this graph, with the same letter,
integer clock and exact content at every step: select its actual multiplier
lift whenever precision rises. This inductive simulation proves that a
potential bound for **all abstract paths** implies the same bound for every
actual word. Spurious paths only make the certificate harder, not unsound.

### Exact state count

For h=0 there are4q/5 choices of z and4 multiplier residues. For each
1<=h<E there are4q/5^(h+1) choices of z and4*5^(h-1) multiplier residues.
For z=0 there are4q/5 multiplier residues. Summing gives

    states(E) = ((16E+84)/25)*5^E.                       (A2)

This is O(E*5^E), compared with O(5^(2E)) for the full affine observer.
At E9 it has17,812,500 states. This is a bound on representation size, not
an end-to-end running-time theorem or a uniform feasibility proof.
The initial all-ones path retains its full precision at z=0, so no artificial
zero-cost reset of its phase or original clock is made.

### Why the branching cannot be silently removed

At q=25, the full states (z,m)=(1,3) and(1,8) have the same adaptive label.
After letter a=10, both have z'=0, but their multipliers are16 and1 modulo25,
giving contents1 and2. The old label does not determine that future answer.
The abstract observer must allow both. These are admissible affine states;
no assertion about a prescribed positive source realizing either state is
needed for the representation counterexample.

This is the same consumer boundary seen in the newly arrived rectangle-only
engine: erasing coordinates is exact only after proving equality of the
actual downstream table. Here exact equality is false, so we use a sound
upper enclosure and keep every possible lift.

## Cheap hostile probe and result

`adaptive_probe.py` constructs the abstract graph from the actual identity,
searches for a nonpositive potential at slope lambda_E, then independently
checks every edge by a separately organized update and a sum of divisibility
indicators. It checks letters1..64 and the all-letter tail
max V+8qE-65(5q-5)<=min V. All states have finite values and max V=0.

| E | states | audited finite edges |
|---|---:|---:|
|1|20|1,280|
|2|116|8,384|
|3|660|52,080|
|4|3,700|311,200|
|5|20,500|1,810,000|

The results in[adaptive_probe.json](adaptive_probe.json) are universal in
word length and positive alphabet **at each of these five precisions**.
They do not extrapolate to all E. No certificate from this exploratory
simulation is needed in the E11 proof, which uses the independently audited
F2 epoch certificate instead.

The next useful question is whether these abstract potentials have a
precision-compatible recurrence. A positive-cost abstract path would
identify an invented high-multiplier lift that must be linked to earlier
history; it would not refute actual Collatz. A successful all-E recurrence
would have to preserve the source identity, integer clock and every
precision-raising branch, or prove why certain branches are incompatible.
