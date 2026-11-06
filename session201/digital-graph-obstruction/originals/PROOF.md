# The source sequence is regular but its binary graph is not

Status: **RESERVED / UNREVIEWED PROOF CANDIDATE**. All-order arguments are
provided below, with exact finite controls. No independent audit has yet been
performed. The twelve-template occupancy assertion remains **OPEN**.

## 0. Target, inheritance, and literature boundary

We retain actual source indices u>=2 and

    g(1)=1,
    g(2h)=g(h),
    g(2h-1)=g(h)+h-(h mod 2).

The source constructor and complete arithmetic inverse are inherited from the
audited Session 201 central arithmetic packet. The factorization follow-up
was published at commit 15354a4b464d558b80eb5ea063c30ea622663cf4. The new
`../session201-digital-inverse/PROOF.md` is a candidate under review, not a
dependency of the proofs here.

A literature search did not identify this exact g sequence or a theorem
covering its J value range. This is a limited search outcome, not a novelty
claim or a claim that no such source exists; direct OEIS search endpoints
were inaccessible. The verified relevant standard distinction is between
2-regular sequences (finite-dimensional linear digit recurrences) and
2-synchronized sequences (finite-automaton recognition of the binary graph).
Definitions and source scope are recorded in `SOURCES.md`.

Concept board: source digits; weighted bit pairs; regular linear recurrences;
synchronized graph languages; actual-value support. The positive inheritance
is the finite-dimensional source recurrence. The hostile near miss is the
incorrect inference that this makes its unbounded-output graph a regular
language. The least-used sidecar is exact separation of graph and support.

## 1. A weighted binary-pair formula

Write m=u-1=sum_i b_i 2^i with b_i in {0,1}, padded by zeros above its highest
one. Then

    g(u)=1+sum_i (1-b_i)(floor(m/2^(i+1))+b_(i+1)).          (1)

Equivalently,

    g(u)=u-b_0
         -sum_(i<j) b_i b_j 2^(j-i-1)
         -sum_i b_i b_(i+1).                              (2)

These are exact digital formulas, not an identification with a named OEIS
sequence. The final correction in (2) matters: adjacent one-pairs receive
weight two in total, rather than the weight one supplied by the double sum.

Proof. Set f(m)=g(m+1). The original recurrence becomes

    f(2m)=f(m)+m+(m mod 2),
    f(2m+1)=f(m),        f(0)=1.

Peeling the least significant bits of m gives (1). Summing the contributions
without the factors b_i uses

    sum_i floor(m/2^(i+1))=m-s_2(m),
    sum_i b_(i+1)=s_2(m)-b_0,

where s_2 is the binary digit sum. Subtracting the b_i contributions gives
(2). All sums are finite.

### Sparse-bit specialization

For every a>=2 and b>=1,

    u=2^(a+b)+2^b+1  =>  g(u)=u-2^(a-1).                  (3)

Indeed u-1 has exactly two one-bits, separated by a>=2 positions; (2) has
exactly one pair contribution 2^(a-1), and no adjacent-pair contribution.

There is also a short recurrence-only proof. At b=1, the parent is even:

    g(2^(a+1)+3)=g(2^a+2)+2^a+2
                =g(2^(a-1)+1)+2^a+2
                =3*2^(a-1)+3.

Here g(2^k+1)=2^k+1 for k>=1, by direct induction. Increasing b replaces an
odd h by 2h-1 and preserves the deficit h-g(h), proving (3). The restriction
a>=2 is necessary: adjacent bits have deficit 2, not 1.

## 2. A finite-dimensional linear digit representation

For the sole purpose of describing a sequence indexed from zero, extend
g(0)=1. This adds no physical source. Put

    V(n)=(g(n),g(n+1),n,n mod 2,1)^t.

Then V(0)=(1,1,0,0,1)^t and

    V(2n)=T0 V(n),     V(2n+1)=T1 V(n),

with

    T0 = [1 0 0 0 0]       T1 = [0 1 1 1 0]
         [0 1 1 1 0]            [0 1 0 0 0]
         [0 0 2 0 0]            [0 0 2 0 1]
         [0 0 0 0 0]            [0 0 0 0 1]
         [0 0 0 0 1]            [0 0 0 0 1].

The first two coordinates use

    g(2n)=g(n),
    g(2n+1)=g(n+1)+n+(n mod 2).

These remain valid at n=0 under the stated extension. Induction on the digit
length proves the representation. Every sequence in the 2-kernel of g is an
integer linear combination of the five coordinate sequences, so g is
2-regular. This proves the assertion directly under the standard definition;
no classification or closure theorem is needed.

For each fixed modulus M>=1, reducing the entire vector modulo M gives a
finite-state automaton computing g(n) mod M. This finite residue-state fact
must not be mistaken for finite-state recognition of the exact graph.

## 3. The binary graph needs infinitely many distinct states

Represent a pair (u,y) by its binary digits read most significant first,
padding the shorter component with leading zeros. Introduce pair symbols

    A=(1,1), Z=(0,0), B=(1,0), C=(0,1).

For positive integers k,m define the prefix and suffix

    P_k=A Z^k,
    S_m=B C Z^(m-1) A.

Their concatenation encodes

    u_(k,m)=2^(k+m+2)+2^(m+1)+1,
    y_(k,m)=2^(k+m+2)+2^m+1.                              (4)

By (3), with a=k+1>=2 and b=m+1>=2,

    g(u_(k,m))=u_(k,m)-2^k,
    y_(k,m)=u_(k,m)-2^m.

Consequently,

    P_k S_m represents a pair (u,g(u)) iff k=m.             (5)

Suppose a deterministic finite automaton recognized the graph. Two distinct
prefixes P_k and P_l would end in the same state. Appending S_k would then
produce the same acceptance verdict, contradicting (5): P_k S_k belongs to
the graph while P_l S_k does not. Infinitely many distinct prefix states are
therefore necessary. The graph is not a regular language, and g is not
2-synchronized. Allowing arbitrary common leading zeros does not affect the
argument, since every displayed word starts with A and has no leading zero.

For completeness, the exact regular-language intersection is

    Graph(g) intersect A Z+ B C+ Z* A
       = { A Z^k B C^r Z^(k-1) A : k>=1, r>=1 }.

The state-separation argument above already proves the claim without relying
on a pumping-lemma convention or this alternative description.

A finite-height consequence is also exact. Any deterministic automaton
correct on every pair of binary length at most 2K+3 needs at least K states,
since P_1,...,P_K are distinguished by the suffixes S_1,...,S_K within that
length budget. Thus an exact bounded-height graph checker cannot stabilize
to one fixed finite automaton as the height grows.

## 4. Consequence for every fixed physical readout

Let alpha,beta,gamma,delta be fixed integers, alpha!=0, and suppose

    h(v)=alpha g(v)+beta v+gamma(v mod 2)+delta

is nonnegative on v>=2. Then h is not 2-synchronized on that domain.

Proof. For fixed integer coefficients, the relation

    t=alpha y+beta v+gamma(v mod 2)+delta                  (6)

is recognizable by a finite automaton on padded binary tuples. One direct
construction reads the digits from the least significant end, retaining
only the bounded carry, the finite initial parity choice, and the finite
constant's digit position. Reversing a regular language preserves
regularity, so the same relation is recognizable with most-significant-first
encoding. Negative coefficients do not cause unbounded carries: with fixed
coefficients the recurrence for the carry divides by two after adding a
uniformly bounded digit contribution.

If Graph(h) were regular, intersect its tuple lift with (6) and v>=2, then
project away t. Closure of regular languages under intersection and projection
would give the graph y=g(v), since alpha!=0 makes y unique. Extra leading
zeros caused by a longer t coordinate are harmless regular padding. This
contradicts Section 3. Removing or adding finitely many initial indices
cannot repair nonregularity.

In particular, every fixed-u row of each of the twelve physical templates
is nonsynchronized: its coefficient of g(v) is strictly positive, by the
reviewed physical readout theorem. For J, fixing u=2 already gives

    Psi_J(2,v)=g(v)+2v.

Hence the three-variable relation n=Psi_J(u,v), u,v>=2, cannot be recognized
by an ordinary finite automaton on their synchronized binary encodings:
restriction to u=2 followed by projection would give a forbidden row graph.
The same fixed-index argument applies to each individual template relation.

### Exact limits of the obstruction

This rules out an exact finite-state verifier for those graphs. It does not
rule out a digit algorithm retaining unbounded integer state, such as the
current digit-and-gcd inverse. It does not rule out finite-state recognition
of the one-dimensional VALUE SUPPORT: projecting a nonregular relation can
produce a regular set. Nor does it show that the unlabeled union of twelve
witness relations is nonregular. In particular it supplies no missing target
and does not refute eventual occupancy.

The familiar representation theorems for synchronized sequences therefore
cannot be imported merely from g's linear digit recurrence. Any different
constructor would need to avoid that unsupported graph-recognition step.

## 5. Reproduction and status

Run:

    python3 session201-digital-graph-obstruction/verify.py
    python3 -O session201-digital-graph-obstruction/verify.py

Both runs check the matrix transitions, digit formula, sparse-bit identity,
and positive/negative distinguishing words using exact integer arithmetic.
The all-order obstruction rests on (3)--(5), not a finite-state search or a
larger occupancy census. The verifier's finite limits are declared in its
JSON reports. The literature description and the new proofs are distinct.

No all-target constructor and no new all-height gap of the twelve-form menu
has been obtained. Uniform occupancy remains OPEN.
