# Independent review of the depth and clock obstructions

2026-11-03. **REVIEWED: no mathematical blocker in the stated scopes.**
Reviewed [DEPTH.md](../DEPTH.md) and [CLOCKS.md](../CLOCKS.md) separately from
their derivation. Their byte hashes at the targeted follow-up are retained
in `collar.json`; this review does not promote a global theorem identifier.

For DEPTH(D1), splitting a Hamiltonian path immediately after the selected
actual vertex is injective into the ordered covers counted by E_w(2).
Therefore each zero-ready layer multiplies the previous H by at least2^s,
so the positive lengths have sum at most L=floor(log2 X). Positive j-tuples
of bounded sum are counted by binom(L,j), and each history has at most
M+L+d physical vertices available for each mark choice. The resulting
upper bound counts construction histories before any quotient, so arbitrary
marking and coincident counts do not invalidate it. Finite phase labels
only select among already allowed operations in this grammar.

The 2^k ternary descendants are distinct: subtract the all-2 word and use
uniqueness of base-three expansions with digits0 or2. Their sizes have
logarithm O(k). Fixed d therefore gives a polynomial versus exponential
contradiction. For varying d, replace it by min(d,L); logarithms of the
stated bound give k log2 <= O(1)+O(log(d+1))+O(d log(k+2)). This establishes
the claimed necessary Omega(k/log(k+2)) depth bound for a fixed finite
bottom bank. It does not exclude changing the bottom bank or operation.

For CLOCKS(C1/C2), subtracting at two distinct s values isolates the 2^s
coefficient. Substitution recovers the constant coefficient, giving both
necessity and sufficiency. The odd-F positive-shift exclusion follows
integrally from 2^k F'=3F. Negative k needs the allowed s+k>=1 domain,
which is already included in the statement about allowed s values.

Each explicit all-positive resynchronization split was checked:

- index0: v2(932)=2; equality of exponents is impossible by sign,
  and the unique surviving possibility is(s,t)=(2,4).
- index1: v2(32744)=3; equal exponents give(2,2), while the unequal
  possibilities at minimum exponent3 fail the stated sign/size bounds.
- index4: v2(39484)=2 and v2(97312)=5; equal shifted exponents fail,
  t=1 fails, and s=2 gives t=3.

For the general compiler, after removing the powers of two from a,b,
unequal shifted exponents have valuation equal to their minimum. This
leaves precisely the two minimum-exponent cases. Equal shifted exponents
give the third direct equation. The zero-right-side case is a separate
unbounded line exactly when the odd parts agree. No bounded experiment
is needed for completeness of this case analysis.

The separate inner-collar experiment does not establish these proofs;
it checks a relevant boundary they deliberately allow. The +4→+4 child
escapes every zero-ready marking of its fourteen-vertex suffix, while its
outer root39 remains certified zero-ready. Thus nested depth, the chosen
carrier, and the actual port incidence must not be conflated.
