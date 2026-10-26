# Independent audit: actual ordinary reverse valuation pass

2026-10-26. **REVIEW PASS / FINITE-EXACT.** No mathematical or implementation
gap was found in the frozen ordinary engine `reverse_ordinary.py`, SHA256
`04d8ac415d8bfec7dd818702e4cee8e718c528248f4f3da2cd100d89e7496312`.
The reverse gauge and precision-budget filter have the claimed consequences.
The structural recurrence is not made smaller by this artifact; the precision
filter reduces modular edge updates after constructing and processing the full
DAG. Optional Newton and skip-contraction artifacts are outside this report.

## Independent routes, not only the producer's assertions

`audit_reverse.py` imports the code under test but reconstructs the following
from literal cycle multisets, factorials and permutation actions:

1. At every stage through N=32, enumerate every unmerged prefix multiset `p`.
   Recompute `S=N-sum(p)` and every affordable defect directly as
   `D(d)=sum_(c in p)(d-gcd(c,d))`, rather than incrementing the implementation's
   stored defect vector. Recompute its raw rational prefix weight

       w(p)=2^[sum_(c in p)(c-1)/2 + sum_(a,b in p; a before b)gcd(a,b)]
             / product_d[d^k_d k_d!].

   Reconstruct its integer gauge contribution as
   `O(sum(p))*2^[S*(len(p)-b(D))]*w(p)`, aggregate only by the complete
   `(stage,S,D)` key, and compare the resulting Y and Z at every vertex.
   Every allowed multiplicity transition and destination key is also compared.

2. Enumerate each complete odd partition through N=32 independently. Its raw
   rational weight gives the ordinary count; its summand valuation is
   `E(p)-sum_d v2(k_d!)`. Trace the unique implementation path and verify that
   the sum of J equals this independently computed valuation and that the
   product of odd T equals `O(N)/product_d[d^k_d O(k_d)]` exactly. For each
   K in {1,2,3,8}, explicitly take the union of edges on complete partitions
   with valuation below B+K and compare it to the proposed edge filter.

3. Independently apply Burnside directly to **all permutations** of N vertices
   for N=0..7. Follow the permutation action on unordered vertex pairs and its
   orientation-flip parity. A pair orbit with odd flip parity admits no fixed
   tournament; otherwise it admits two orientations. Average the exact fixed
   tournament counts over S_N. No prime, connected, clone-free, core, or
   automorphism restriction is imposed. This is independent of the producer's
   labelled-tournament orbit marking and odd-partition oracle.

4. Through N=50 and additionally N=60,80,100, independently recompute all exact
   integer T/J values, Bellman equalities for V and P, structural edge identities,
   and odd-unit modular factors. Compare every modular coordinate against both
   exact Z and an unfiltered pass at K in {1,2,3,8,17}.

All audit gates are explicit exceptions; neither audit relies on `assert`.

## Why the operation and mergers are valid

Every nonterminal vertex records the decreasing odd stage, remaining degree,
and all affordable future defect coordinates. Future stages and remaining
degree only decrease, so a discarded coordinate can never be requested later.
The displayed direct prefix reconstruction proves the intended merger on the
finite tested range; the inherited scalar defect-gauge identity supplies its
general continuation proof. A zero residue is never a structural reachability
criterion. Every generated vertex has an all-singleton completion.

The ordinary multiplicity choice commits `ell*k` degree, not `ell*(k+e)`;
there is no excess parameter or signed beta coefficient in this engine.
Setting e=0 is **not** a clone-free restriction: the unrestricted tournament
Burnside formula still includes tournaments with clones. In particular it gives
one class at N=2 and two at N=3, not the inherited core counts zero and one.

The edge's odd-factor ratio is a positive odd integer. Indeed
`(u+ell*k)!/[u!*ell^k*k!]` counts the choice of the new cycle support and its
permutation with k cycles of length ell, so its odd denominator divides the
numerator's odd part. The nonnegative integer J is the inherited ordinary
specialization of the integral defect-gauge formula. Every edge in the finite
test range was checked separately with factorial-derived odd parts and
factorial-derived valuations.

Let V be the backward minimum J-cost to the terminal and B=V(initial).
Then `J'=J+V(target)-V(source)>=0`. A prefix path of cost c admits a suffix
of cost V at its endpoint, so `c+V>=B`. Its transformed contribution is an
odd integer times `2^(c+V-B)`, an integer **before** merging. Thus

    Z(v)=2^[V(v)-B]Y(v)

is valid on the reachable prefix lattice, including vertices where the
displayed coordinate exponent is negative. The potential is state-defined,
so it commutes with every exact `(stage,S,D)` addition. It is not claimed to
preserve integrality on arbitrary unrelated integer vertex assignments.

Every complete path corresponds to exactly one odd cycle partition. Its
factor products telescope to O(N) times that partition's rational Burnside
summand. Thus B is the minimum summand valuation, not necessarily the final
count's valuation, and the terminal equalities are

    Y(t)=O(N)*a(N),
    Z(t)=O(N)*a(N)/2^B.

Final multiplication by the inverse of the odd O(N) modulo 2^K returns the
normalized answer without dividing an already truncated absolute residue.
The inherited theorem identifies B=floor(N/2) for N>=3 and B=0 for N=0,1,2;
the implementation computes the graph minimum and then checks this identity.

## Why the precision filter is exact

P is the forward minimum reduced prefix cost. A minimum original suffix has
reduced cost zero by telescoping. Consequently the minimum reduced cost of a
complete path traversing a fixed edge `v -> w` is exactly

    P(v)+J'(v,w).

Both directions hold: the minimum prefix and zero-cost suffix can be joined
through the edge because the state merger is continuation-compatible; and no
prefix or suffix has cost below the corresponding minimum. Keeping precisely
the edges with this quantity less than K therefore gives precisely the union
of complete paths whose original valuation is below B+K. Every discarded
contribution vanishes modulo 2^K. This does **not** imply that every retained
contribution, or their sum, is nonzero modulo 2^K.

This also explains why the implementation's filtered intermediate values agree
with all exact Z residues: Z(v) is divisible by 2^P(v), and the dropped incoming
transfers are individually divisible by 2^K. No cancellation heuristic enters.

## Concrete hostile controls and wording boundary

* **Absolute-zero information loss:** at N=3, vertex `(stage=1,S=0,D=())`
  has Y=2, V=0, B=1, and Z=1. Y is zero modulo 2 while Z is not. At N=50,
  a(50) is zero modulo 4 but a(50)/2^25 is 2 modulo 4. Both are exercised.
* **Missing stage:** at N=6 the actual edge
  `(3,1,(0,)) -> (1,1,(0,))` has k=0. Erasing stage naively identifies its
  endpoints, producing a self-loop in the supposed decreasing DAG. The two
  vertices have different stage-local transitions. Their *completed* suffix
  sets coincide after the deterministic skip; an explicitly proved contraction
  of unit zero-cost skips is therefore a different, potentially valid operation.
  This audit does not overstate the hostile as different complete suffix sets.
* **Subtracting B only initially:** actual initial zero-cost skips would get
  negative shifts. The order-50 control explicitly verifies this failure.
* **State merger:** the receipt records a nonterminal exact merger of distinct
  cycle multisets. The checker audits all mergers, not just this witness.
* **Boundary inputs:** negative, Boolean, noninteger, null, and string orders or
  precisions are rejected; empty and singleton orders have their exact values.

The producer clarified during review that “K-bit arithmetic” means **K-bit
stored residue precision**. Python multiplication and left-shift temporaries
can exceed K bits before reduction. This clarification is present in the
frozen API docstring and implementation document; no algorithmic change was
needed. The modular engine stores neither a B+K-bit absolute-state array nor
the full exact tournament count.

## Finite receipts and honest accounting

The independent audit passes normally and under `python -O`, with identical
JSON receipts. Its coverage is:

| Check | Coverage |
|---|---:|
| Exact edge T/J and structural checks | 105,864 edges |
| Exact prefix coordinates | 32,875 |
| Modular full-coordinate comparisons | 270 |
| Independently reconstructed unmerged prefix histories, N=0..32 | 8,198 |
| Independently reconstructed merged prefix coordinates | 2,563 |
| Complete partition-path checks, N=0..32 | 2,765 |
| Exact budget edge-support equalities | 132 |
| Unrestricted pair-action Burnside | N=0..7 |

The independent unrestricted class counts are
`1,1,1,2,4,12,56,456`; at N=7 the labelled universe has 2,097,152 tournaments.
The independent audit computes fixed-point counts over all 5,040 permutations
there, not an explicit enumeration of all those tournament masks.

The producer's complete program was also replayed normally and with `-O`.
Those replay outputs match each other and the supplied producer receipts.

At N=100 there are 13,445 allocated vertices including the terminal and 50,044
generated edges. The filtered modular pass performs 568 edge updates at K=2
and 2,576 at K=8. These are genuine modular-update counts. They exclude full
graph generation, reverse/forward potential preprocessing, per-vertex traversal,
filter tests, odd-factorial preparation, and allocation of the full value array.
There is no elapsed-time, peak-memory, asymptotic state bound, end-to-end speedup,
or smaller unrestricted structural recurrence conclusion in this artifact.

## Reproduction

From this audit directory (the default sibling source path is recognized):

    python audit_reverse.py > audit-normal.json
    python -O audit_reverse.py > audit-optimized.json
    cmp audit-normal.json audit-optimized.json

For an installed copy, pass `--source /path/to/directory/containing/reverse_ordinary.py`;
the two producer receipt files must be in that same directory. Producer replay:

    python ../session64-counting/reverse_ordinary.py --repo ../session64-repo > producer-replay-normal.json
    python -O ../session64-counting/reverse_ordinary.py --repo ../session64-repo > producer-replay-optimized.json
    cmp producer-replay-normal.json producer-replay-optimized.json
    cmp producer-replay-normal.json ../session64-counting/reverse-normal.json

No repository file was edited by this audit.
