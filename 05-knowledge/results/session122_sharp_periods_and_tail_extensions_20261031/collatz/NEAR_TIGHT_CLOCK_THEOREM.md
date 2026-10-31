# Restoring exact clocks on the first two slack layers

This extends `CRITICAL_CLOCK_THEOREM.md` from exact equality to integer
potential deficit at most two. It is a general finite-certificate method,
applied to the independently checked E10 certificate in `E14_CLOSURE.md`.

Let V be the integer potential and define edge slack

    d(s,a,t)=V(t)-V(s)-8*c(s,a,t)+5*a >=0.

For an actual prefix, its deficit is V(state)-(8*S_E-5*A). Deficits add
exactly along edges. Therefore all paths terminating with deficit b in
{0,1,2} are represented by the backward graph on augmented states
(state,b), with predecessor budget b-d(s,a,t). All compatible inverse
multiplier lifts and all allowed letters are tested, not sampled.

## Exact all-one phase loss

Write k=b*T_ell+s with0<=s<T_ell and T_ell=2*5^(ell-1). The canonical
entry has the same retained state as the long entry. Let s5(n) denote the
base-five digit sum. The exact difference between their integer scores is

    score(s)-score(k)
      =2*s5(b mod5^(E-ell+1))+10*floor(b/5^(E-ell+1)).

To prove it, put n=floor(k/2), epsilon=k mod2. Then

    5*k-8*sum_(j=0)^(E-1) floor(n/5^j)
      =2*s5(n mod5^E)+10*floor(n/5^E)+5*epsilon.

The low ell-1 digits are unchanged when replacing s by k; the remaining
ones are precisely the digits of b. Thus nonminimal phase costs at least
two. At total deficit at most two the entire phase bank is exactly:

- the canonical k=s;
- k=s+T_e for ell<=e<=E, each costing exactly two extra units.

The second family is available only when the canonical entry itself has
zero potential deficit. This is why an unbounded entry phase does not
create an unbounded search at these two slack layers.

## Finite graph and exact pair propagation

The reverse-edge equations are those in `CRITICAL_CLOCK_THEOREM.md`.
At each candidate predecessor, compute its exact nonnegative integer slack
and keep it iff it does not exceed the current budget. Canonical entries
are tested against exact floor sums; each permitted longer entry is checked
again by its literal floor sum and multiplier phase.

The graph is checked to be acyclic. This is necessary for the stated
finite path calculation; resource termination or a cycle is a failure to
certify, not an empty answer. On its verified topological order, each entry
starts with pair(A,L)=(k+a,k+1). Each edge adds(a,1). Sets of exact pairs
are propagated and deduplicated, preserving the existence of each possible
clock and length, without identifying word orders or translations.

Every omitted letter a>=65 has reduced slack greater than two whenever

    Vmax+8E-325 <= Vmin-3,
    8E-325 <= Vmin-3.

These strict bounds, not merely non-strict all-letter certificate tails,
are checked from the actual potential range before the graph is extracted.

For an actual word whose full endpoint is identity modulo5^K, K>=E, every
terminal pair must also satisfy3^L=2^A modulo5^K. If none of the propagated
pairs satisfies that equality, the entire deficit layer is excluded for
that actual endpoint class. A surviving pair would be only a necessary
condition: translation, order and native returns would remain to be tested.

## Pilot and validation

The Python and C++ implementations agree on the E7,R2 graph through budget
two:32079 nodes,32111 edges, all retained entries, all terminal pairs and
abstract path counts. They differ only in the number of scanned candidates
because C++ additionally skips required predecessor values outside the
verified potential range. The three score layers42,43,44 have respectively
5,0,70 abstract paths and2,0,12 distinct clock/length pairs. None passes
the pilot's identity multiplier test modulo5^10. This is a pilot for the
method, not the requested E14 conclusion.

## Completed E10 application

The full E10,R2 graph has3908743 augmented states and3908754 edges,
independently reconstructed by a separate Kahn-order implementation. Normal
and NDEBUG replays agree. Its score30,31,32 layers have respectively2,0,38
abstract paths and1,0,11 clock/length pairs. All twelve pairs already fail
the identity multiplier test modulo625. Together with the independently
checked all-letter potential, this proves the actual identity10 endpoint
bound H10>=33/8. The exact E14 deduction and open boundaries are recorded
in `E14_CLOSURE.md`.
