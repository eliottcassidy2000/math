# The exact doubling clock inside the odd-girth construction

PROVED corollary of the Session25 Walecki/paired-triples construction. This
is a conjugacy of its constraint skeleton, not of tournament Hamiltonian
counts or classical Collatz dynamics. No literature-priority claim.

Use n=4r+1, r>=1, k=2r+1, and M=r+1. Index each repaired pair of paths by
a in1,...,r. Its two path indices are k-a and r-a, while the third path
containing its repeated matching edge has index k-2a. The unpaired base
path has index r. If k-2a=r, the pair points to the base. If k-2a>r, it
points to pair2a; if k-2a<r, it points to pair2a-r-1. Therefore, identifying
the base with0, the functional constraint skeleton is exactly

    a -> 2a modulo M.

The base0 has no outgoing constraint edge; its formal fixed loop is omitted.
This explicit correspondence preserves parent, eventual-cycle membership,
period and entry depth. It forgets the actual path vertices, arc directions,
and edge-equation signs. Those must still be restored for a tournament cycle.

For period and depth statements below,0 is the formal fixed point of the
doubling clock. In the actual constraint skeleton its loop is omitted, so
a root-bound trajectory reaches the base rather than a retained cycle.

Write M=2^s u with u odd. CRT shows the eventual-cycle vertices are exactly
the u residues whose2-primary coordinate is0. A nonzero a reaches its
cycle after exactly max(0,s-v2(a)) steps, with v2(0)=infinity. On the odd
factor, an element of additive order d has doubling period ord_d(2).
There are phi(d) such elements. Hence the number of nonzero cycles is

    C(M)=sum_(d|u, d>1) phi(d)/ord_d(2).

The root component is a tree; each other component has one cycle. Thus
C(M) is exactly the cycle rank of the expanded path-reversal constraint
graph: it has2r+1 vertices and2r edges, and1+C(M) components. The orientation
argument must satisfy one independent binary compatibility equation for
each nonzero cycle. Some may already be satisfied; C(M) is not a count of
failed equations. The spare-path ordering repair resolves all of them.

For every compatible path order, the reversal equations then have exactly

    2^(1+C(M))

solutions, one free global flip per connected component. These are reversal
assignments for this fixed construction, not counts of distinct tournaments,
metagraph cycles, Hamiltonian paths or isomorphism classes.

In particular C(M)=0 iff M is a power of2. Exactly at the orders
n=4*2^s-3 (s>=1), this construction has no cycle-compatibility obstruction
for ANY path order and exactly two reversal assignments per order. If a
nonzero cycle exists, moving its spare path between its active endpoints
toggles its obstruction, so at least one of the two orders is incompatible.

Concrete arithmetic instances, with the quantities typed:
- M=21 gives n=81, C=5 and64 reversal assignments per compatible order.
- M=63 gives n=249, C=12 and8192 assignments.
- M=341=11*31 gives n=1361. The three nontrivial divisors11,31,341 have
  orders10,5,10 under2, contributing1,6,30 cycles. Hence C=37 and2^38
  reversal assignments.

The occurrence of21,63 or341 is now an exact modular-clock calculation in
a named construction. It does not identify this clock with Collatz, pay
any153N debt, or explain forbidden actual Hamiltonian counts7 and21.

Reproduce with python3 verify_doubling.py. The verifier compares the direct
path-index parent formula with doubling, independently counts functional
cycles and entry depths, and compares the divisor/order formula for every
2<=M<=1000. It does not enumerate the enormous labelled tournament graphs.
