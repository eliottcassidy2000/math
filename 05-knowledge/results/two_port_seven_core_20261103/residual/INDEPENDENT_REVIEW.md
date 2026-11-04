# Independent review of the complete order-ten modular sieve

**Status: ACCEPTED / FINITE-EXACT.** Reviewed `check.cpp`, `run.py`,
`audit.cpp`, `audit.py`, the included prior census engine, and the full,
prefix and normal/optimized audit receipts. The producer's complete
order-ten result is:

    all classes tested                         9,733,056
    S=E modulo2 classes                           12,670
    also C=C^t modulo4 classes                         0

This is a necessary-condition sieve, not a reconstruction-collision census.
The complete run is supplied by `full10.json`; this review independently
checks its source and bank pins, all decoded addresses, the recurrence
proof, the integer audit, and the prefix replay. It does not claim a second
complete execution of the order-ten sieve.

## The fast parity filter cannot discard a valid candidate

For a nonempty vertex subset A, let E(A,v) be its integer Hamilton-path
count ending at v. The exact recurrence is

    E({v},v)=1,
    E(A,v)=sum_(u in A-v, u->v) E(A-v,u)  if |A|>1.

The producer stores, for every subset A, a ten-bit word whose v-th bit is
E(A,v) modulo2. Its predecessor mask is exactly
1023 xor (1<<v) xor outneighbors(v). Taking the parity of
`dp[A-v] & predecessors(v)` performs the displayed sum over F2. All proper
subsets are computed earlier because deleting a set bit lowers the mask.
The singleton case is implemented by `!rest`. Induction on |A| therefore
proves that every fast endpoint bit is exact. Repeating the same recurrence
on the converse gives the start bits on the original tournament.

Consequently `e==s` is exactly the endpoint-parity condition in B3, with
no possible false negatives from the fast predicate. Checking only
retained classes against integer DP would not by itself establish that
fact; the recurrence proof is essential. The additional independent
integer audit corroborates the actual implementation on2,123 controls:
every labelled core of orders1-5 embedded as an ordered prefix of an
order-ten tournament (1,099 controls), plus1,024 specified pseudorandom
order-ten tournaments. It compiles the actual producer function extracted
from `check.cpp`, rather than a separately transcribed fast routine.

For all12,670 retained bank classes the complete run additionally computes
integer S,E,C and checks every parity bit against the fast result before
testing C skew modulo4. Each C_uv is the sum over nonempty proper subsets
A of E(A,u)S(V-A,v), precisely the required ordered two-path cover count.
No endpoint marginal is substituted for the joint matrix.

The integer ranges are safe. The sum of all two-path-cover counts is at
most9*10!=32,659,200, since a permutation has nine cut positions. Individual
endpoint and cover entries and each nonnegative partial accumulation
therefore fit uint32. The unsigned subtraction in `(C_uv-C_vu)&3` gives
the correct residue modulo4 even if subtraction wraps, because4 divides
2^32. It is not a signed-order comparison.

## Universe, decoding, and source integrity

Completeness is inherited from the
[consecutive-flags bank proof](../../consecutive_flags_order10_20261028/GENERATION.md)
and its full alternate canonicalization receipt
`../../consecutive_flags_order10_20261028/verification10.json`.
That audit establishes uniqueness and labelled mass2^45; agreement with
the scalar class count alone would not suffice. This review does not
rerun the inherited canonicalization of the complete bank.

The bank contains191,536 sorted nine-vertex upper-triangle codes in
little-endian u64 format. An address splits as
`(parent_index<<9)|root_outmask`; a root-mask bit means new vertex->old
vertex. This agrees with the included `decode` implementation. Addresses
are stored as unsigned base128 deltas, first measured from zero.

An independent Python decoder in this review checked every bank address,
strict increase after the first, parent bounds, exact end of stream,
parent-code order/range, and the complete count9,733,056. The last address
is98,066,431. Both compressed and uncompressed bank SHA256 values matched
the manifest, and the compressed pins matched `full10.json`. The full
producer likewise requires191,536 parents, exactly9,733,056 records, valid
varints and addresses, and EOF immediately after the last full-run record.

The property being tested is invariant under simultaneous relabelling of
S,E and C. One representative of every ordinary class therefore suffices
to exclude every named order-ten core. No representative is selected by
an unproved modular or degree filter.

Current source hashes match every full/prefix receipt:

- `residual/check.cpp`: `7fbf7c6a67a8329af55ed82d7f6f2d0d3e820d27dd7108cb1e090d6fdd8e0fc3`
- `residual/run.py`: `4ba69bdfbdd4d8bae463270279dacee6dd5bd1d5ed3af15f5d9c24dfb9d46f7a`
- prior `census/census.cpp`: `9a3eabd8c5a85869ae22899c5f5de99682d670e777bc46567f5979d428de22b0`

Bank compressed hashes are
`ef5cf3d13f5b679890cbe69bcd1f78189161953e726e12d169fb576845677ef7`
for the flags and
`92cb19bbcdc0bcb58b0ded354aacdd1809c72251e2e29d28d23f92de2d8bbe05`
for the parents. The wrappers check their sources both before and after
execution and check compressed and uncompressed bank hashes.

## Replayed controls and the exact consequence

Independent normal and optimized executions of `audit.py` both matched
their frozen receipts, including all2,123 integer comparisons and all
source pins. An optimized `run.py --limit 10000` replay matched both
normal and optimized prefix receipts:10,000 classes,10 endpoint-parity
balanced, zero combined survivors. No check depends on Python assertions.

Combine the complete order-ten exclusion with
[B3-B4](../algebra/PROOF.md), accepted in the independent
[algebra review](../recovery/ALGEBRA_REVIEW.md). The inherited complete
[endpoint census](../../two_port_beyond_nine_20261103/audit/endpoints.json)
has zero endpoint-parity balanced classes at orders2,4,6, and exactly two
at order8, codes39887362 and72392840. Direct integer counts independently
checked here give, respectively,

    C_02-C_20 = 162-12 = 150 = 2 mod4,
    C_02-C_20 =  66-44 =  22 = 2 mod4.

Thus both order-eight possibilities fail B4. Every odd module of size at
least three is excluded at every size by the inherited
[odd-module current theorem A1-A2](../../two_port_beyond_nine_20261103/algebra/PROOF.md).

**Corollary (PROVED, with the stated FINITE-EXACT dependencies).** For any
named tournament M with2<=|M|<=11 and any ambient quotient Q, substituting
M at a vertex and reversing exactly M changes the full two-port response.
The exterior quotient has no size bound and need not be strong or
self-converse. More generally, a module-converse counterexample at any
size must have an even, nonselfconverse core of size at least12, satisfying
both intrinsic congruences. Even self-converse cores are excluded at all
sizes by B3 and inherited P4.

This statement reverses **one uniform module while keeping its exterior
fixed**. It does not establish arbitrary two-port reconstruction at
order11, does not exclude unrelated changes of arcs, and does not silently
extend the result to simultaneous reversals of arbitrary even modules.
The full two-port response contains additional three-path information;
this necessary-condition exclusion already works through its H,S,E,C
sidecar. Universal two-port reconstruction remains OPEN.
