# Scope check beyond the complete degree-twelve census

**Status:** Structural equality range proved from the inherited pronormal-base theorem. Four selected numerical probes are FINITE-EXACT and are not a complete subgroup census at degrees thirteen or fourteen. The probes use the extension's compressed engine and its global-sign/isomorphism primitives; they have normal/optimized agreement, but no separate uncompressed replay is claimed here.

Use the groups D and P from README §5, on p+3+f vertices. The all-f relative identities in §5.5 do not themselves imply equality. A separate inherited theorem does extend equality past the clean f≤2 multiplicity range:

**Proposition.** If 0≤f≤p−4, then tau_T(D)=tau_E(D) and tau_T(P)=tau_E(P).

**Proof.** The degree n=p+3+f is less than 2p, so the displayed Cp is an ambient Sylow-p subgroup and is pronormal. Both D and P lie between Cp and its ambient normalizer, and their quotients by Cp have orders 3 and 9 respectively. The pronormal-base/one-quotient-prime theorem applies with q=3. This does not require P itself to be pronormal. QED.

For p=7 the range includes f=3, or n=13. The selected degree-fourteen probe is the first possible loss of this particular ambient-Sylow sufficient argument.

| n | embedding | R | tau_T | tau_E | F_T | F_E | Z |
|---:|---|---:|---:|---:|---:|---:|---:|
| 13 | P=F21×C3 | 176 | 172 | 172 | 4 | 2 | 2 |
| 13 | diagonal F21 | 352 | 172 | 172 | 180 | 2 | 178 |
| 14 | P=F21×C3 | 1408 | 1376 | 1375 | 32 | 24 | 9 |
| 14 | diagonal F21 | 2816 | 1376 | 1375 | 1440 | 24 | 1417 |

Thus p=7 shows that the larger equality range is sharp for this particular prime: both embeddings first fail at f=4. This is not a claim that every prime p first fails at f=p−3, nor that these are all failing groups at n=14.

At n=13 the product is no longer fusion-free: F_T=4 and F_E=2, while Z=2 compensates the difference. This is why the clean statement m_T(D)=2 and m_E(D)=1 was restricted to f≤2. At n=14 the diagonal group's large escape count 1417 is almost entirely compensated by source fusion 1440 and target fusion 24; the actual defect is only one.

Reproduce with:

    python3 family_scope_probes.py > family_scope_probes.json
    python3 -O family_scope_probes.py > family_scope_probes.optimized.json
    cmp family_scope_probes.json family_scope_probes.optimized.json

All four selected rows explicitly verify the all-f relative identities. These probes test the boundary of a proved family law; they are not used to infer any complete classification beyond degree twelve.
