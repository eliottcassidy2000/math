# Exact one-candidate-per-cylinder inverse for the S83 fixed bank

Status: PROVED elementary algorithm, with FINITE-EXACT implementation tests. Chan–Pak and F2 completeness remain OPEN. No new global coverage range is claimed.

## Inheritance and actual obstruction

The carrier is exactly THM-4875(D), strengthened by HYP-9478's S83 addendum: F_y(v)=a g(v)-b e(v)+c(v+1), v>=0, y=(a,b,c) in the 2,732 distinct rooted triples of order at most seven. The inherited coverage through 2^35 is a finite computation, not induction. Here g(0)=1, g(2u+1)=g(u), g(2u)=g(u)+u+(u mod2), and e(v)=1-(v mod2), INCLUDING v=0. Physical bank coordinates satisfy c>=a>=b>=0 and a odd positive.

The closest old numerical inverse is Session217 `inverse/SUFFIX_INVERSE.md`: square-root inversion of the different fixed-three-core family. The low-bit-cylinder technique and its limitations are inherited, not claimed new. This packet specializes it to the actual strongest S83 family and supplies a particularly simple one-candidate certificate without a length loop.

The decisive inherited hostile is Session220 `selectors/PROOF.md`, Sections 1–7. It already applies to THIS EXACT F2 family: every rational finite-state witness relation in ordinary binary (including asynchronous transducers, unlabelled finite-bank choices, and finite-state recodings) has bounded-run sources and only polylogarithmically many outputs. Thus residue automata or rational selectors cannot settle this bank's all-target assertion. Its positive mixed-period drift, not just a technical choice of synchronous reading, is the obstruction. An unbounded arithmetic register or a different non-rational encoding is needed.

Concept board: exact physical bank; finite-state drift; suffix coefficients; magnitude-constrained inversion; modular-vs-integer occupancy. The new certificate retains the actual source integer and exact final equality. It never replaces an integer witness by a congruence witness.

## 1. Exact suffix identity

For k>=1, M=2^k, 0<=r<M, let rev_k(r) reverse exactly k binary digits, and t=[r>=M/2]. Then for every u>=0,

    g(Mu+r)=g(u)+(M-1-rev_k(r))u+(1-t)(u mod2)+g(r)-1.       (1)

This is the inherited append law, now in closed reversal form. To prove it directly, append the k digits of r from left to right. Every appended zero adds the previous prefix plus its parity. The coefficient of u is the sum of 2^(j-1) over zero positions j, namely M-1-rev_k(r). Incoming parity matters only if the first appended bit is zero, giving 1-t. All later parity contributions are fixed by r. Evaluating at u=0 identifies the constant as g(r)-1. Leading zeros and u=0 are valid, since appending a leading zero to (0,1) preserves that state.

Put C=M-1-rev_k(r), E=1-t, and

    L=cM+aC,
    Q=a(g(r)-1)-b e(r)+c(1+r).

Then exact target membership in this cylinder is

    x:=N-Q = L u+a g(u)+a E(u mod2).                       (2)

Here L>=cM>0 and Q>=0 for the physical bank. All terms other than the already-subtracted constant are nonnegative. No unphysical source state is introduced: every nonnegative integer u defines the valid source v=Mu+r.

## 2. At most one integer remains

Induction gives 1<=g(u)<=u+1. Equation (2) consequently implies

    max(0,ceil((x-2a)/(L+a))) <= u <= floor((x-a)/L).       (3)

If x<a, there is no candidate. For x>=a the unclipped real interval has length

    (x-a)/L-(x-2a)/(L+a)=a(x+L-a)/(L(L+a))
       <= aN/L^2+a/L.                                    (4)

Choose M to be the least power of two >=2 satisfying

    c^2 M^2 >= 4aN,     cM>=4a.                           (5)

Then (4)<=1/2<1. Thus each low-bit suffix r has AT MOST ONE integer u in (3). Test it in the exact equation (2). The union of all M suffixes is a complete membership decision, because every v>=0 has exactly one decomposition Mu+r. No maximum source length needs to be guessed.

Per bank label, M=O(1+a/c+sqrt(aN)/c). A fixed finite bank therefore admits O_Y(sqrt N log N) arithmetic operations using direct recurrence evaluation of g; all integers have O_Y(log N) bits. Suffixes can be streamed with O_Y(log N) working bits apart from a fixed bank. This is numerical, not bit-polynomial complexity. The square-root strategy is inherited from Session217; this exact bank-specific implementation is the present adaptation.

The C++ executable restricts inputs to 1<=N<=10^15. The supplied bank has a<=189, b<=45, c<=432. For N<=10^15, (5) implies M<2 max(2,4a/c,2 sqrt(aN)/c), hence M<1.8e9 and L,Q<1.2e12. Candidate sources obey v<=N/c+M<1.000002e15, so even 189(v+1)+432(v+1)<6.22e17 stays below 2^63. These bounds place signed 64-bit readouts safely in range for the supplied bank; threshold multiplications and final verification use signed 128-bit arithmetic. Abstract theorem (1)–(5) uses unrestricted integers. Bank provenance is the repository's `04-computation/experiments/rooted_triple_bank_order7_20261104.txt`, read from the parent's frozen HEAD; physical realizability of its entries is inherited from THM-4875 and its audit, not reconstructed here.

## 3. The precise remaining all-target statement

For every odd N other than 7 and 21, at least one y and one r in the explicitly bounded list (5) must have a nonempty interval (3), whose unique candidate satisfies (2). This is an exact equivalent arithmetic formulation of F2 completeness, with all source and magnitude constraints retained. Proving merely that a candidate interval is nonempty, or that (2) holds modulo any modulus, is insufficient.

Three attempted routes and stopping reasons:

1. Finite-state residual descent is blocked in ordinary binary by Session220's positive mixed-word drift theorem. No new no-go is claimed.
2. Fixed-bank radix closure fails because appending a suffix leaves the g(u) coefficient a unchanged; an identity to 2^j F_z(u), j>=1, would require the physical odd coefficient a to equal 2^j a_z. This is also inherited Session220 Section5, and does not exclude target-dependent arithmetic.
3. Target-adaptive low-bit descent genuinely evades those restrictions because L grows without bound. It yields the complete finite certificate above, but supplies no reason that one branch has zero exact residual. This is the surviving occupancy obligation, not an algorithmic blocker.

## 4. Reproduction

    g++ -std=c++17 -O3 inverse.cpp -o inverse
    python3 -I -B verify.py > CHECK.json
    python3 -I -B -O verify.py > CHECK.optimized.json
    cmp CHECK.json CHECK.optimized.json

`CASES.txt` records selected sparse target probes, including values above the inherited 2^35 cutoff. They are isolated witness checks, NOT a larger interval-coverage claim. The checker independently brute-enumerates the bank through 2048 and compares every odd target with the inverse; it independently evaluates all displayed witnesses and 51,510 suffix identities. It checks bank size, uniqueness and required coordinate inequalities. No asserts are used.
