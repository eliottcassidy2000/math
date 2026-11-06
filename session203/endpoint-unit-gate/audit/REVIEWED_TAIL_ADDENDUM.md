# The rooted unit-gate obstruction survives every source-count cutoff

**PROOF CANDIDATE, pending independent audit.** This strengthens the source quantifier in PROOF.md without changing its context class.

Fix any integer t>=1. Let U_ab^(t) be the one-apex tournament with transitive private chain and root-incidence word

    1, 0^t, a, b, 0,      (a,b) in {0,1}^2.

The first and last bits still supply a spanning directed cycle and a unique root-starting Hamilton path. The two varied bits are adjacent and interior, so the identical twin-contraction proof applies to any two fixed rooted-module contexts. In particular the mixed difference of child minus an even multiple of parent remains2 modulo4.

All four input H counts tend to infinity with t. One may see this without using a scalar fit. Their common induced prefix has word1,0^t and Hamilton-path count2^t+1. It is a cyclic triangle whose middle block is the transitive t-chain: a Hamilton word places that chain either in one run (three choices of block order) or in two nonempty runs separated by the other two vertices (2^t-2 assignments), giving3+(2^t-2)=2^t+1. Induced tournament path insertion therefore gives H(U_ab^(t))>=2^t+1.

Consequently, for any cutoff N_0, choose t with2^t+1>=N_0. The four inputs all have H at least N_0, but the proposed constant-correction identity fails on at least one of them. Thus the theorem rules out an even-multiplier constant gate on all sufficiently large strong unit-start sources, not only on the four small inputs in the main proof.

The cores can even be chosen separately for each input length t, provided the same parent and child contexts are used for the four assignments at that length. Their sizes do not enter the mod-four contradiction.
