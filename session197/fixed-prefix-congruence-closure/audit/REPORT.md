# Independent audit: complete congruence closure of a fixed rooted prefix

**PASS, no correction required.** The exact modular image, compatible odd lifts, explicit Euclidean construction, cardinality and common-parity-branch requirement are correct. This is complete local congruence admissibility for a fixed physical prefix, not exact-value sufficiency or a new cofinite spectrum theorem.

## Frozen source binding

Source: `session197-prime-atoms/fixed-prefix-congruence-closure/`.

- PROOF.md: `2754e0de2ab154fa47918969cce43040866683914143047d996c5db2f95fe360`.
- verify.py: `755df735d670ea66666547e17f2ed92fefe939b72812fe6f86bfbc154a1c5516`.
- Both supplied and independently copied normal/optimized source receipts: `8ca02a9c30554ccee256e1f95bcc5c1be2057a515bbbff8e90682de2fd363d8d`.

The original files were not changed. Source replays ran in an audit-only copy and matched exactly.

The named joint modular gadget theorem was read at its exact published commit48e10c6a4773de23f9ca19ef120b17012a0d7812, path session197/modular-gadget-universality/PROOF.md, SHA74df1989d75ff743177937c108ac36dce133412e48701e4d412297d8c47b7d63. The inherited physical arithmetic law was checked in ARITHMETIC_INVERSE.md, SHA4ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7.

## All-order argument

For a fixed physical prefix, h is odd, R=e modulo2, R>=h>=1 and0<=e<=h. Put d=gcd(h,R), D=gcd(d,m). Both d and D are odd. The joint gadget theorem modulo2m genuinely allows independent v,w modulo m in L=2v+t and g(L)=2w+1, with t fixed. Thus the lifted image modulo2m consists exactly of the two cosets

    h+2D Z,  h+R-e+2D Z.

Their constants are odd. Their reductions modulo D are0 and−e. If m is odd, every residue modulo m has exactly one odd lift modulo2m. If m is even, precisely the odd residues have odd lifts; adding m preserves parity and, since2D divides m, preserves the same lifted coset. Therefore the projected image is exactly the displayed compatible-parity union of cosets. It is not the larger set obtained by choosing unrelated branches at different prime divisors.

For sufficiency, choose an odd representative N and t with D|(N+t e). With K=h+t(R−e), N−K is even and divisible by D. Since D is odd, b=(N−K)/2 is still divisible by D. The congruence h w+R v=b modulo m is soluble exactly under this condition. Bezout coefficients for h,R reduce it to d s=b modulo m. The inverse of d/D modulo m/D exists; when m/D=1, s=0 is valid. L and g are then reconstructed as the prescribed parity-compatible residue pair modulo2m. This proves the constructor's sign and normalization, not merely a numerical residue count.

Arbitrarily large witnesses follow from the already proved modular gadget loop. Their output counts are unbounded because physical prefixes satisfy h>=e, g(L)>=1 and R>=1, hence F(L)>=LR. The argument would require these physical bounds; it does not silently treat arbitrary signed triples as physical prefixes.

Testing the moduli2 and d proves necessity of the full congruence-closure criterion. Its sufficiency for every modulus follows because every D divides d. One global choice of t works whenever N belongs to a specified coset modulo d, though the constructed integer L may depend on the modulus. This does not construct a single L with F(L)=N.

The two cosets merge exactly when D|e. Their sizes after projection give the stated formula, including modulus one. Every power of two has D=1 and hence the complete odd image. The physical prefix(33,1,33), target143, correctly passes the tests separately modulo3 and11 but fails modulo33: its residues are2,0,11. This is a concrete counterexample to independent primewise parity-branch selection.

## Independent physical and modular controls

The standalone checker imports no source implementation. It computes full subset path counts and anchored two-cover counts on all labelled tournaments through order five, giving5405 literal rooted profiles. These reduce to42 distinct triples; a separate literal seven-vertex one-back-edge prefix supplies(33,1,33), for43 distinct physical profiles in the modular controls. This bank is a physical stress test, not a claim that every one of its roots belongs to the narrower inherited grammar.

It reconstructs the actual one-apex gadget for every L=1,...,128. It checks H=g(L), exactly one path ending at the root, and ending-two-cover count L. That ending response is kept distinct from the prefix's starting-two-cover response R. It then constructs472 actual seam graphs and independently checks their Hamilton counts against the physical formula.

For48 moduli, including all1 through40 and several larger prime powers and composites through101, a breadth-first search enumerates the complete binary-word state automaton(n modulo2m,f(n) modulo2m). Transitions are actual appended binary digits from the singleton initial state, rather than arbitrary independent residue variables. Across130838 reachable states, every expected odd-f state occurs. Projecting these reachable states through each physical prefix gives2064 exact modular-image/cardinality checks, all matching the theorem.

The checker also recovers986 actual integer gadget indices from shortest automaton paths and evaluates their counts by a separate left-to-right binary recurrence. These are different witnesses from the source's Euclidean/long-loop constructor. Their residues agree exactly. The named examples and the composite-branch hostile are included explicitly.

## Reproduction and hashes

Run the standalone checker without inputs:

    python3 check.py
    python3 -O check.py

Both outputs are byte-identical.

- check.py: `dceaf5bf878621cee6a8b3e0a75c8aa8281a2e3eb4c056758417600b92f43454`.
- independent-normal.json and independent-optimized.json: `b70b2ea1c93bca4470c4101bae06984f3fa1f510cba4712fcfd54d8f198846c3`.

The supplied source's7712 residue/cardinality checks and986 separate Euclidean-constructed witnesses also replay exactly. The finite tests support the implementation; the physical joint universality theorem and coset proof establish the all-order statement. No source correction or additional scope restriction is required.
