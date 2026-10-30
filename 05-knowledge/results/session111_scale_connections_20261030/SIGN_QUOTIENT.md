# All four numbers occur in one exact action

**PROVED for the signed and unsigned normalizer actions; FINITE-EXACT check at modulus25.**

This is the sharpest connection found in the session. The numbers25,50,100 and200 are not merely parallel examples: they are kernels and images of two representations of the SAME5000-element group.

## The two representations

Let R=(Z/25)[δ]/(δ²−5), μ=(1+δ)/2, H=<μ⁴>, |H|=25. Let N be the additive normalizer of H, of order5000.

Its action on orientation choices retains both the permutation of the24 antipodal clock-orbit pairs and the signs that flip their orientations. The kernel of this signed action is H, and its image E has order200.

Now forget the signs and retain only the permutation of the24 pairs. The kernel grows to100, while the image shrinks to50:

| representation | kernel in N | image |
|---|---:|---:|
| signed orientation action |25|200|
| unsigned coordinate action |100|50|

Thus5000=25·200=100·50. The additional factor4 is exactly two independent binary signs, not an unexplained numerical match.

## Exact kernel proof

From PISANO_GROUPS, N/H=E≅C20×C10. Write c for ring conjugation and −1 for scalar negation.

Every nonzero vector can be written δ^d s g^j, where s is a scalar unit and g=μ⁴. Conjugation sends it to
(−1)^d δ^d s g^(−j).
Therefore c preserves every unsigned ±H-orbit. Scalar−1 also does so.

These induce independent signs:
- scalar−1 flips every orientation coordinate;
- conjugation flips precisely the odd-valuation coordinates.

Both even and odd valuation layers occur. Hence they generate a four-element kernel K≅C2² in E.

There are no additional unsigned-kernel elements. Fixation of the pair containing1 forces the multiplier u into±H. Remove that multiplier. Fixation of the pair containingδ then gives
(b−epsilon g^j)δ=0
for epsilon=±1, where the remaining ring automorphism is δ↦bδ. The scalar b−epsilon has even δ-valuation; a nonzero epsilon(g^j−1) has odd valuation. They cannot cancel before both vanish moduloδ³. Thus b=epsilon modulo25, leaving only identity or conjugation.

Consequently the unsigned kernel in N is exactly four cosets of H, of order100, and
E/K≅C10×C5
has order50. The unsigned-kernel group itself is C2×D50, with D50 of order50.

## Exactly two bits restore the action

Every unsigned permutation in the image has four signed lifts. If the task is to recover WHICH signed normalizer action occurred, two bits are necessary and sufficient.

Two named probes suffice: retain the orientation signs on the orbit of1 and the orbit ofδ. Global negation changes both; conjugation changes only the latter. They distinguish the four lifts.

This is a scoped information theorem. It does not say every downstream decision requires two bits, nor does it identify orientation assignments with tournament isomorphism classes.

## All-depth version and the boundary

For the golden tower R_k, k≥2, signed kernel order is5^k and unsigned kernel order is4·5^k. Signed image order is8·5^k and unsigned image order is2·5^k. The same valuation proof and the same two probes work at every depth.

The corresponding order200 clock/conjugation group C4×D50 from PISANO_GROUPS remains a different group from the signed image C20×C10. The present diagram uses the full additive normalizer and explicitly named representations, so those two meanings of200 are not conflated.

[sign_forgetting.js](sign_forgetting.js) enumerates all200 effective signed actions, verifies exactly50 distinct unsigned permutations and lists all four kernel sign patterns. Their binary rank is2. No full graph automorphism census is claimed.
