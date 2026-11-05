# Independent audit: exact C2 spectrum with feedback number at most two

**PASS. No correction requested.** The exact class is strong tournaments with exactly two Hamilton cycles and feedback vertex number at most two. The result is neither a classification of all C2 tournaments nor of all fixed-port tournaments with arbitrary cycle count. In particular its H533 exclusion does not exclude a fixed-port H533 graph having three cycles.

## Frozen source binding

The source is `session197/fixed-two-c2-spectrum/`. All eleven entries of `RELEASE_FILES.sha256` were checked before copying. Key hashes:

- `PROOF.md`: `5703612b53a46453b1947f0e3c209c01608d3eb88dc585e36213b7d5f46ad5e5`.
- `verify.py`: `d68091d9d0232fbb574f299bd11db0a97dc78233125c2a6b0a3a490734fa8e9c`.
- `certificate.py`: `565642b3aaa854db80a50d6a86806b2a3c132fdd3ceedab316f32029f32ba47c`.
- `CONTROL.json`: `789d4a346b70aee2ead2c5db495271e146c0ae00b978984c66cc6cf31079ed9d`.
- `VALIDATION.json` and `VALIDATION_OPTIMIZED.json`: `ed87faff25868359e9d6335f0d97a7c7f5fb78f52c2d19b8db41d7de6dd23562`.

Copied normal and Python-optimized replays reproduce the complete frozen validation receipt byte for byte. Original source files were not changed, and their final key hashes remain unchanged.

The mathematical inputs are the independently reviewed physical N/M/P families, the fixed two-port cover transfer and strong-induced cycle monotonicity, and the pair-entry/marked-clone identities. I read the exact named proofs and checked their use here. Their frozen hashes are recorded in the source proof. This audit does not turn the separately proved strong-induced cycle monotonicity into a claim about arbitrary weak intermediate graphs.

## All-length language completeness

The fixed-port representation is exactly the claimed graph restriction: deleting two retained vertices leaves a transitive graph, and any feedback deletion set of size at most two can be extended to two vertices and used as the retained ordered ports. The remaining transitive ordering supplies the four mask letters. Strong C2 graphs have enough vertices for this extension. Conversely, every fixed-port word has that deletion certificate.

The converse-word involution is correct. After graph converse, swapping the two ports and reversing the nonport order sends masks0 and3 to each other and fixes1 and2. It preserves H and C; it need not preserve any particular distinguished response root.

Coordinate clipping at3 commutes with the nonnegative integral path-cover transfer and the cycle readout. Thus the controller distinguishes cycle counts0,1,2 and at least3 exactly at arbitrary depth. A rejected prefix has a Hamilton cycle and is strong. Strong-induced monotonicity compares it directly to any proposed final strong C2 graph, so pruning remains valid even if an intermediate graph is weak. Acceptance must be remembered as the last edge's output, as the source correctly does.

I independently checked the regular-language equality without using the source's Thompson NFA or its minimized product:

- `check_language.py` rechecks every one of the 1,644 raw controller transitions, including the absorbing reject state, using explicit readout indices and sparse multiplication from the unchanged audited transfer.
- It verifies reachability of all411 raw states, the physical seed, and the source's15-class edge-output quotient congruence.
- The seven displayed languages and their optional branches are expanded into24 direct glob patterns including converse patterns. Stars are loops at the current state; the construction has110 states and no Thompson epsilon machinery.
- This independent NFA is crossed with the full raw controller, not its minimized quotient. The complete product closes at425 states and1,700 transitions. At every reachable state, NFA acceptance agrees with the remembered last-edge output. All four outgoing transitions are included.

This closed finite product proves the displayed language equality for every finite word. It is not a comparison through some maximum length. Independent normal and Python-optimized outputs agree.

## The physical Q family and its marked root

The proposed first-clone profile is the inherited adjacent first-clone profile. Applying J to that root gives exactly

    H=3D+6u-2p-2,
    e=D+2u-p-1,
    R=3D+7u-2p-3,
    d=D+3u-p-2.

The marked vertex is the newly added J vertex b. It is not either old clone. The subsequent retained-root binary suffix has the same scalar value as the ordinary one-seam gadget on the pre-J source: substituting J's four rows into the retained suffix formula cancels the constant terms and leaves

    g(v)(D+3u-p-2)-(v mod2)(u-1)+v(2D+4u-p-1).

This also handles v=1, where the actual graph is the J seed and no suffix is present. Its start count is two, and cutting a Hamilton cycle at the last seam gives C=2. Thus Q is an all-order physical family rather than a formal profile assignment.

The independent C++ checker builds V(u) directly from its apex incidence bits, adds the two old clones, then adds the new b beating the first clone. It also constructs the fixed-port word independently and checks the explicit relabelling, preserving the new b as port1. On u=2,...,128 and v=1,...,8, all1,016 marked maps and70,559 unordered arc checks pass. This construction does not import the source's graph-building or map functions.

Separate endpoint inventories verify63 Q seed profiles through u=64. On132 Q graphs (u=2,...,12 and v=1,...,12), endpoint Hamilton DP and literal full Hamilton-word enumeration agree with Q, with exactly two cycles and an explicit transitive deletion of the two ports. These finite controls corroborate the all-order graph and suffix argument above.

## Terminal reductions and scalar scope

The A0^k3 reduction is correct: each sink leaves h fixed, doubles the two reserved-cover responses, makes the selected ending counts zero, and preserves the two directed K counts by appending the sink to the other cover slot. The final pair birth gives precisely N(u,2^k).

For B1, the retained-apex O then Z calculation gives R_a=4u-1. Deleting b and reinserting it at z gives R_b=2D+2u-p. In K_ab the b-start path can only be b or b,z because the other path ends at a, so continuing through a is forbidden. These give1 and p. For K_ba, z must end the a-start slot; after deleting b,z the only forbidden reserved-a cover is the singleton a, giving u-1. Hence K_ab+K_ba=u+p and the second terminal reduction gives M(u,2^k).

These two reductions are correctly stated as value reductions, not marked-graph identities. The empty M suffix reduces to N(u,2), while P's empty prefix/suffix boundaries correspond to parameter1. The isolated word013 contributes the regular-five value15. All converse words preserve the scalar value.

An additional230 independently built terminal graphs (u=2,...,24 and k=1,...,5, both terminal forms) pass both endpoint and literal Hamilton counting, with C=2. This does not replace the all-order cover argument, but checks every gate and both terminal formulas independently of the source verifier.

## Exact inverse and separation

The N/M/P symmetry and parameter bounds are valid inherited physical statements. Their weaker rectangular regions in the code still include every solution: N gives uv<=(h-1)/3; M uses u+v>=4; P uses u+v>=2. They correctly allow parameter1 only for P.

For Q, the physical pre-J source has h-e>=2u and R>=4u+1. The odd-u case uses g(u)>=3; the even case permits g(u)=1. Therefore Q>=2u+v(4u+1), giving exactly the finite u and v limits in the code. Q is not incorrectly symmetrized. The target1 boundary yields no strong C2 witness, and regular5 is added at15 independently of the four families.

My C++ checker enumerates an independent full rectangle u,v<=1861, without applying the source inverse regions. It uses literal binary-apex path assignments for g rather than the source recursive implementation. Across8,655,511 admissible family-parameter tests it records3,223 exact solutions with h<=1861, including regular5. `check_inverse.py` compares every resulting parameter set with the source inverse for every integer target1,...,1861, including even targets. All sets agree, in normal and Python-optimized modes. In particular533,821,1861 have no parameter solution.

The explicit11-vertex H533 graph is independently counted by endpoint DP and literal Hamilton words: H=533 and C=2. Every one of its55 pair deletions contains a directed triangle, while deleting {2,5,8} is transitive. Its feedback number is therefore exactly three. Together with the all-order inverse exclusion at533, this proves strict separation of the two-feedback C2 value set from the full C2 value set. No C1 exclusion, least-gap assertion, or unrestricted-spectrum conclusion is required or inferred.

## Independent artifacts

All paths below are relative to this audit directory:

- `check_language.py`: `b61d86159f01b357d7ecec40b1c8c2caafbd3f8d828aceb3c45dc0603e496e96`.
- `language-normal.json` and `language-optimized.json`: `3d68a72844383d89a5c09905016853041bbc1591e1375b956666968256460498`.
- `check_graphs.cpp`: `7a630791129974706b33d4cbf140cf705beb5f28a0c9aedf2d0b86c9d00170e8`.
- `graphs-normal.json` and `graphs-debug.json`: `95f9e43476c77cd1bab1e4e50d8293023391dcc747e931e2a4bc282509904dc8`.
- `inverse-normal.csv` and `inverse-debug.csv`: `0bae5d0e076e715d5546fcc7737689409a5ee4fb53b37265e4f9787480f8bc8b`.
- `check_inverse.py`: `fca88fb0603b2077ef954704f86f924d83d49483ec45b94c8a08dd1d0db5b7fd`.
- `inverse-check-normal.json` and `inverse-check-optimized.json`: `354a9e6440e94aa5094de796789b483c8c85fe6c9527d71949c49b5d6ae26773`.

The C++ checker takes the output CSV filename as its argument. The Python checks use the preserved `source/` snapshot. Its two-port transition module is an explicit audited dependency of the language check, not a newly independently rederived transfer. The graph and arithmetic checker imports no source code. C++ optimized/debug and Python normal/optimized receipts agree. No unresolved issue or corrective edit remains.
