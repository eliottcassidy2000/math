> **Curated publication copy.** Executable proofs and certificate data are in
> `code/`; independent audits are in `audit/`. Run `verify_all.py`.
> Historical exploratory scripts mentioned below are intentionally omitted;
> their finite searches are not dependencies of the certified theorems.

# Two fixed ears over a transitive spine: exact transfer and compiler obstruction

Status: PROVED HERE (elementary derivations); FINITE-EXACT verification below.
No assertion that the tournament spectrum conjecture is settled.

## Carrier and exact transfer

Vertices are x,y,v1,...,vm, with x→y and vi→vj for i<j.
Letter c∈{0,1,2,3} records bits (x→vi,y→vi), with x the low bit.
For each ear order (a,b)=(x,y) or (y,x), every Hamiltonian path is uniquely
P,a,Q,b,R, where P,Q,R are increasing subsequences partitioning the spine.

Maintain P and Q states 0=empty,1=good last vertex,2=bad last vertex;
P-good means last(P)→a, Q-good means last(Q)→b. A nonempty Q must start at
an a-outneighbor. R has states0=empty,1=started, and may start only at a
b-outneighbor. There are18 states per order,36 total. Start at(0,0,0)
in both order blocks. Append each vertex to exactly one of P,Q,R:

* P resets to1 if vi→a, otherwise2.
* Q may start only if a→vi; when allowed, resets to1 if vi→b, otherwise2.
* R may start only if b→vi; once started stays1.

These three alternatives are summed, not composed. The acceptance row sums
P∈{0,1}, R∈{0,1}, and Q∈{0,1} for order xy, Q=1 for order yx.
Thus middle-empty is allowed precisely when a→b. This proves the transfer
by a bijection, including the empty spine (count1).

`check.py` implements the matrices and compares them to an independent
adjacency/subset Hamiltonian-path DP for all 5,461 words of lengths0..6.
Normal and optimized receipts agree.

## Correct affine identities and a full-state synchronization obstruction

An exact scalar accumulator ray is v(N)=Nu+w, with output(u)=1, output(w)=0,
M_d u=2u, M_d w=w+d u, for d=±1. Common-fixed-offset identities of the form
M_-v=2v-z, M_+v=2v+z, M_±z=z do NOT supply that invariant ray.

For every actual nonempty-spine state, each order's Q=R=empty slice is a
single unit recording the last spine vertex's orientation to the first ear.
A letter acts on this slice by a rank-one reset, with spectrum{0,1}.
Consequently a common2-eigenvector u vanishes on both slices. The ray equation
(M_d-I)v=(N+d)u forces the current last-letter bit to equal the appended
letter's bit, in each ear order. Two distinct letters cannot satisfy this.
This is only a full-state obstruction; it does not by itself exclude
observable quotients or phase-dependent rays.

## Global scalar no-go theorem (all fixed equal-length complete digit banks)

**Theorem.** Fix L≥1 and B=2^L. There do not exist a seed on this carrier with
H>1 and B fixed length-L spine words, one assigned to each odd digit
D={-(B-1),-(B-3),...,B-3,B-1}, such that appending the word assigned d always
sends H to BH+d, on the seed and EVERY descendant formed from those words.
The obstruction is to scalar counts themselves: observable quotients and
phase bookkeeping cannot repair it while retaining these same fixed words
and their universal assigned digit actions.

**Proof.** Distinct digits require distinct words. After k templates, all
resulting counts are at most B^k(H(seed)+1), i.e. O(2^m) in the number m=Lk
of appended vertices.

Let A be the union of letters occurring in the templates. If both bit
coordinates vary on A, there are letters giving each required condition:
x-outneighbor, y-outneighbor, x-inneighbor, y-inneighbor. Periodically
concatenate all templates. In an early bounded prefix reserve distinct
vertices for first(Q), an x-outneighbor, and first(R), a y-outneighbor.
In a late bounded suffix reserve distinct vertices for last(P), an
x-inneighbor, and last(Q), a y-inneighbor. All old spine vertices go in P;
all intervening free vertices may be assigned independently to P,Q,R.
Early and late unused vertices can be fixed in P and R respectively.
This constructs at least3^(m-C) accepted xy paths, with C fixed independent
of the number of repetitions. This contradicts O(2^m).

Otherwise at least one bit is constant, so |A|≤2. Since there are B=2^L
distinct words, |A|=2 and the bank contains ALL length-L words over A.
There are four possible edges of the bit square:

1. A={1,3} or{2,3}: the bank contains3^L, which appends universal sinks
   and leaves H unchanged. BH+d=H would require |d|=(B-1)H>B-1.
2. A={0,1}: y beats no newly appended spine vertex. If y beats an old spine
   vertex, put that vertex first in R, start Q at a new1, and end P,Q at
   new0 vertices. Long intervening0 runs give3^(m-C) paths, contradiction.
   If y beats no old spine vertex, y is a universal sink throughout;1^L
   adds universal sinks immediately before y and leaves H unchanged.
3. A={0,2}: if x beats an old spine vertex, use ear order yx, that old
   vertex first in R, a new2 first in Q, and new0 vertices last in P,Q.
   Again long0 runs give3-growth. It remains that x beats ONLY y.

For the remaining case put S=H(spine+y), t=number of its paths ending at y,
and V=number of ordered partitions(P,R) of the spine into increasing
subsequences with R empty or y→first(R). A path in the full tournament
has either x last, counted by S-t, or x immediately before y, counted by V.
Thus H=S-t+V. After appending r≥1 letters2, the tournament spine+y gains
universal sinks, so S is unchanged, t=0, and V is multiplied by2^r.
Therefore repeating the template2^L k times gives

    H_k = S + B^k V  (k≥1).

If its assigned digit is d, the required recursion instead gives

    H_k = B^k(H_0+d/(B-1)) - d/(B-1).

Comparison for unbounded k forces S=-d/(B-1)≤1. But we may first append
2^L and then0^L, creating a directed triangle y→v2→v0→y in spine+y.
A tournament containing a directed triangle has at least three Hamiltonian
paths: the three triangle paths extend injectively under successive vertex
insertions. Hence S≥3 at this descendant, contradiction. ∎

### Scope and remaining escape routes

The theorem is NOT a nonattainment theorem for integers. It excludes a
specific universal compiler architecture. It does not cover variable-length
or count-dependent word choices, phase-dependent choices of DIFFERENT
physical templates, changing the two marked ears, adding nontransitive
blocks, or a compiler that does not realize every digit on every descendant.
In particular finite phase-dependent *state representations* are excluded
only when the physical template for each digit remains fixed.

## Historical finite searches (not shipped; not proof dependencies)

`actual_rays.py`: all21,845 words≤7; no full-state single-letter ray seed.
`observable.py`: row-space closure of the output under chosen matrices has
dimensions10 for{0,1},10 for{0,2},18 for{1,2}. All87,381 seed words≤8:
no affine accumulator in the exact future-output quotient.
`radix4.py`: for full-state four-template banks, the last-letter reset forces
a common last letter c, so the bank is{0c,1c,2c,3c}. All87,381 seeds≤8:
not even one scalar first-step digit set{4H-3,4H-1,4H+1,4H+3} occurs.
`rays.py`: unconstrained rational corrected ray systems can admit artificial
solutions; these are NOT reachable-state evidence.

Reproduce with Python3+SymPy using `python -B` and `python -O -B` on each
script. No external tournament code or data is used.

## Extension: finite-phase physical template choices are also excluded

**Corollary.** The same no-go holds if the physical length-L word for a digit
may depend on a phase s in a FINITE set, with next phase determined by(s,d),
and every reachable phase offers all B digits with the required scalar action
on every actual state arriving there.

**Proof.** Take a reachable bottom strongly connected component of the phase
transition graph. Every phase there has B distinct outgoing physical words,
since two distinct digits cannot give different H values for the same word.
Let A be their union alphabet. If both bits vary on A, a closed walk visiting
all edges exists by strong connectivity. Repeating it gives the same
3^(m-C) free-allocation lower bound, contradicting binary-rate growth.
Otherwise A is a binary edge, and cardinality forces every length-L binary
word to be available at EVERY phase in the component. The previous case
analysis therefore remains valid despite phase changes. In the sole-
outneighbor case, once S≥3, consecutive uses of the available word2^L satisfy
H_next-B H=-(B-1)S, outside the permitted digit range. ∎

This strengthens the earlier fixed-bank theorem. The earlier escape-route
paragraph applies to that theorem alone: the corollary additionally excludes
finite phase-dependent physical choices under its explicit hypotheses.
Variable-length words, infinitely many/count-dependent phases, changing
marked ears, nontransitive appended blocks, and partial digit domains remain
outside these results.

`check_no_go_identities.py` checks639 exact decomposition, repeated-letter,
triangle, and mixed-growth identities/inequalities by independent subset DP;
normal and optimized receipts agree. This tests the delicate sole-outneighbor
case rather than treating the formal transfer as its own independent oracle.
