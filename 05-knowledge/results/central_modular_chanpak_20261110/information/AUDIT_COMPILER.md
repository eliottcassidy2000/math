# Independent audit of the address compiler

**VERIFIED proof audit; FINITE-EXACT computational controls.** This audits
[the compiler proof](../compiler/PROOF.md), not unrestricted Hamilton-count
coverage or a complete Chan–Pak induction.

The suffix obstruction has the correct quantifiers. Fixed suffix appends
and flank interchange pull a bounded low-bit predicate back to another
bounded low-bit predicate. Any fixed finite sequence of policy choices can
therefore be fixed on a finer initial binary cylinder. The four-point
stencil then extracts the deletion-core coefficient, and forces its
multiplication by 3 at every uniform ternary return. No assertion about
unbounded guards or arbitrary full-address automata follows.

The odd-multiplier carry formula is exact: for input digit a and channel r,
the decomposition ma+r=2j+b has 0<=j<m. The parity term is
(mn+j) mod2=(n+j) mod2 because m is odd. Its one-hot coordinate is
epsilon_(1-(j mod2)); constants are supplied by epsilon_0+epsilon_1=1.
Thus no extra affine coordinate is missing. The m+3 coordinates retain
unbounded integer values; they are not m+3 finite automaton states.

`audit_compiler.py` extracts only the producer's `digit_matrix` function
through its syntax tree. It does not import the compiler module, recursive
g, state constructor, reader, or inherited graph routines. A separate
iterative implementation evaluates g, an independently written binary
reader applies the tested matrices, and the carry formula is checked
channel by channel against direct iterative g values.

The explicit universe is ten odd multipliers
1,3,5,7,9,15,17,21,25,35 and eight addresses
0,1,2,7,18,123,1024,170807. This gives 80 full-state comparisons,
160 digit-transition comparisons and

    2 * 8 * (1+3+5+7+9+15+17+21+25+35) = 2208

channel checks. Leading-zero invariance is checked for all ten multipliers.
Prime and composite controls both pass. These finite checks corroborate
the proof; they do not replace its all-order quantifiers.

Dropping oddness has a small explicit hostile:
m=2, n=1, digit a=0. The odd-only matrix formula gives
(6,3,2,1,0), while direct evaluation gives (5,3,2,1,0).
The first failed implication is (mn+j) mod2=(n+j) mod2, at j=0.
The strongest survivor is the stated formula for every fixed odd m;
an even-m formula would instead need its actual parity correction.

Reproduce from the repository root:

```bash
python3 -B 05-knowledge/results/central_modular_chanpak_20261110/information/audit_compiler.py
python3 -B -O 05-knowledge/results/central_modular_chanpak_20261110/information/audit_compiler.py
```

Both runs give identical JSON, with source and auditor SHA-256 hashes.
The receipt is `audit_compiler.json`; the optimized receipt is
`audit_compiler.optimized.json`. The original compiler checker separately
recounts its physical graph examples by an independent predecessor DP;
this audit does not duplicate that graph enumeration.
