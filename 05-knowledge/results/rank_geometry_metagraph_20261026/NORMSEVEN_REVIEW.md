# Parent integration review of the separate norm-seven continuation

2026-10-26. **ACCEPTED within its stated closed-fibre scope.** The separately
requested chat completed and published commits `830271623e` and
`a3b6805487`. Parent /root read
[RESOLVED_FIBRES R1-R4](../normseven_descent_moments_20261026/RESOLVED_FIBRES.md)
after publication and independently replayed `verify_packet.py --replay`:
20 manifest files,19 source pins,4 paired receipts and8 normal/optimized
replays passed. This is a proof review and exact computational replay,
not a formal proof-assistant certificate.

The incoming normal-surface/tangent-cone mechanism is inherited from the
separately audited arithmetic packet. Fixing the rational coarse point,
the new twists have four closed nodes with residue K=Q(sqrt(13)). The
exceptional conic is A²+B²=delta C² over that actual field. The review
checked the following new steps independently:

1. For delta2, [1:1:1] gives a K-point; a point field must contain K, so
   degree2 is necessary and sufficient.
2. For delta3,13 has a simple square root modulo3 and hence in Q3.
   Thus K embeds in Q3. A nonzero solution of A²+B²=3C² over Q3 can be
   scaled to integral coordinates with a unit. Reduction modulo3 forces
   A,B divisible by3, then the equation forces C divisible by3, a
   contradiction. Hence there is no K-point.
3. Every point field contains K and therefore has even degree over Q.
   Degree2 fails in the second case. Q(sqrt(13),sqrt(2)) has degree4 and
   supplies [1:sqrt(2):1], proving the exact minimum4.
4. At p outside {2,3,13}, each nonsingular exceptional conic has q+1
   points over F_q. The residue torsor supplies eight rational nodes
   exactly when13 is a square. Therefore both twists have identical
   counts over every good finite-field extension. Summing those counts
   gives the two displayed local zeta factors; these are factors of the
   specified resolved fibre, not whole-surface zeta functions.
5. A decomposition into nonempty split K-tori would supply a K-point
   from the identity of any piece. Its absence obstructs that specific
   relative torification. Over Q even the first fibre has no point;
   the statement cannot be moved from K to Q by forgetting the residue
   field. No assertion about all F1 formalisms follows.

This is consistent with arithmetic local-global principles: the two
models differ at a bad place excluded by the good-prime count statement.
The useful connection to the other lanes is a precise limitation on
summary data. Full good-prime counts of these fibres do not retain a
tangent norm obstruction; scalar tournament statistics do not retain
cycle incidence; finitely many centre distances do not retain an omitted
affine pullback. These are separate demonstrated losses. No map between
the corresponding open problems, or completeness-style equivalence,
has been proved.
