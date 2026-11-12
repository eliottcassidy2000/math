# Two-place exact reconstruction and the precision carried by a recursive state

**PROVED IN NOTE; implementation controls are FINITE-EXACT.** No claim of
polynomial-time full enumeration, new Burnside formula, or new Chan–Pak coverage.

## Provenance and purpose

Session270 `enumeration/PROOF.md` already supplies odd-factorial clearing and
raw/projective response coalescence. THM-5098
`precision-cycle-cutoff-and-normalized-tournament-congruences` already gives
precision-first two-adic sector selection. Session270 `bridges/PROOF.md` gives
an even-order low-bit divisor query. Session271 supplies accurate huge-order
forecasts, but not exact integer recovery from relative error. Historical
`04-computation/a000568_hybrid_crt.py` reconstructs from many odd primes using an
estimated bit budget. Here the modulus is selected by a **proved absolute
interval**, not an estimated bit length plus a fixed safety margin.

The new deliverable is their compatible real/dyadic completion: a finite
explicit algorithm which gives exact integers at every finite order, with a
certificate of how much residue precision the approximation actually needs.
It is not necessarily faster than the existing full GMP counter. The structural
sidecar is precision debt: a compression valid in exact integers need not be a
valid quotient at fixed modular precision.

## 1. Termwise 2-adic integrality and one-inverse compilation

For an odd cycle partition lambda with multiplicities m_j let

    z = product j^m_j m_j!,
    e = sum m_j(j-1)/2 + sum j binom(m_j,2)
        + sum_(i<j) m_i m_j gcd(i,j),
    v = e - sum v2(m_j!),       d = z / 2^v2(z).

Burnside says a(n)=sum 2^e/z=sum 2^v/d. All d are odd. Each v>=0:
the within-equal-length term j binom(m_j,2) is at least v2(m_j!),
including j=1 and m_j=0,1,2; every remaining term is nonnegative.
Therefore each rational summand has an unambiguous residue modulo 2^b:

    0                              if v>=b,
    2^v (d^(-1) mod 2^(b-v))        otherwise.             (1)

Their sum is the ordinary integer a(n) reduced modulo 2^b. This does not say
each rational summand is an ordinary integer.

For D=odd(n!), every d divides D, because z divides n! (a permutation
centralizer order). Thus all surviving terms can instead be accumulated as

    S = sum_(v<b) (D/d) 2^v mod 2^b,
    a(n) = D^(-1) S mod 2^b.                              (2)

This performs just one inverse, and only exact integer divisions. It combines
inherited odd-factorial normalization with precision truncation. The direct
partition implementation supports both (1) and (2) as independent controls.

Its group-by-group branch pruning is safe: inserting k cycles of length j
adds valuation k(j-1)/2+j binom(k,2)-v2(k!)+k sum m_i gcd(i,j), which is
nonnegative. Once the accumulated valuation is >=b, every completion vanishes.

## 2. Real interval plus dyadic residue gives the exact count

Let L<=a(n)<=U be any certified rational bounds, and define

    lo=ceil L, hi=floor U, b=bit_length(hi-lo), M=2^b.

The integer interval is nonempty because it contains a(n). Its width is
strictly less than M, including b=0 if lo=hi. Given the exact residue r modulo
M, the unique answer is

    lo + ((r-lo) mod M).                                  (3)

Uniqueness follows because distinct integers congruent modulo M differ by at
least M. The candidate must lie <=hi; the program checks this consistency
condition, and rejects an undersized modulus or inconsistent interval.

`complete` uses the support-cut interval proved in sibling
`absolute/PROOF.md` and implemented by `certified_interval.py`. It obtains no
exact reference count while reconstructing. The Session270 cache is consulted
only by test programs after reconstruction. Thus this is a genuine exact
algorithm, not fitting an error correction to known answers. Its guarantee is
for each finite input, not bounded runtime independent of n.

## 3. Raw response recurrence is a safe modular state machine

The response-coordinate and degree-local normalization are precisely inherited
from Session270. We implement the raw recurrence over Z/(2^b). Its transitions
multiply by integers and shift left by

    k(j-1)/2+j binom(k,2)-v2(k!)+k q_j >= 0.

At the terminal 1-cycle block the shift is

    binom(r,2)+r q_1-v2(r!) >= 0.

Hence every operation descends to the quotient ring: reducing intermediate
weights modulo2^b and dropping zero weights preserves the final residue.
There is one final multiplication by odd(n!)^(-1), an existing unit.
States with identical retained response merge exactly. The implementation
uses full integer response coordinates; only weights are reduced.

This is a new finite-precision backend for the inherited recurrence, not a new
response theorem. It avoids the slower direct-partition inverses and makes
the real/dyadic method practical through the tested range250.

## 4. Why exact signed gauges require a precision sidecar

A tempting optimization is to take Session270's signed projective-gauge
recurrence and reduce every state modulo2^b. This is generally invalid: exact
right division does not descend to that quotient. For example 0 and 2^b are
congruent modulo2^b and both divisible by2, but their halves differ modulo2^b.
Exact divisibility of a full integer does not repair lost high bits.

**Precision-debt lemma.** In a finite acyclic arithmetic transition graph with
edges x -> u_e 2^s_e x (u_e odd, s_e any integer), assume every actual division
is exact. To compute terminal sums modulo2^b, assign terminal precision b and
back-propagate

    B(state)=max(0, max_edges [B(destination)-s_e]).         (4)

Keeping each state's residue modulo2^B(state) is sufficient. Along an edge,
congruence modulo2^B becomes congruence modulo2^(B+s_e); (4) makes that at
least the destination precision. Merging sums preserves congruence. Induction
on the DAG proves the claim. Additional even factors can be included in s_e.
This is sufficient, not necessarily minimal when cancellations are known.

A gauge w'=2^g(state)w illustrates the issue: precision b for w generally
requires precision b+g(state) for w'. A state compression can reduce the number
of states while enlarging the precision that each must carry. Here we choose
the safe raw backend; we do not advertise an unimplemented signed-gauge speedup.

The same abstract warning applies to Chan–Pak: an exact state quotient for one
readout need not preserve all later interfaces. The map here is concrete
(residue projection), the lost information is high bits before division, and
the sufficient sidecar is (4). It is an analogy to interface completion, not
a theorem that digit precision supplies missing Hamiltonian-count seeds.

## 5. What remains hard

The positive real truncation cannot silently discard all smaller sectors and
round: the sibling obstruction makes this quantitative. Conversely, low-bit
formulas alone do not determine the high bits. The two views together give
exactness, but the required residue precision is still large. The generated tables record the actual figures; no asymptotic savings claim
is made.

The exact printed answer has Theta(n^2) bits. A complete explicit-output
algorithm therefore needs Omega(n^2) bit-output work. This implementation
still enumerates a growing space of response states; no polynomial state bound
or universal finite structural alphabet is proved. Changing the real cutoff
trades more real-sector work against fewer residue bits; it does not abolish
that combinatorial work. The present advance is certified completion and
compatible recursion, not instantaneous arbitrary-index exact output.
