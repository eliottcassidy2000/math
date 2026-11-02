# Independent structural audit: light-mass dual closure at n>=526

Status: ACCEPTED structurally, conditional only on the explicitly inherited
Session138 inequalities and separately audited exact constants in
`../stability/PROOF_N526.md`. This audit does not claim to reprove all prior
matching and binary equality results.

1. The assigned greedy groups partition their edge sets. Counting selected-
edge wedges and within-endpoint-star wedges therefore creates no cross-group
double count. The inequality q_e<=w_i gives the simultaneous invoice
`W_i>=w_i S_i/2+S_i^2/4`. It is compatible with the separate `W_i>=w_i S_i`.
2. Keeping B in the nonlinear margin, rather than discarding it, correctly
leads to `V-D>c/kappa-n kappa/2` by a completed square. Removal of light mass
U then leaves a positive heavy matching. The integer inference r_h>=18 from
r_h>17 is valid and guarantees m=R-L>=18.
3. The stronger retained-matching coefficient applies because n>2m. Its
monotonicity follows from its representation as an average of the decreasing
geometric sequence. The replacement of r_h bym therefore has the correct
inequality direction.
4. Retaining U exactly yields `cL-4U<4 epsilon V+4D`; the strict inequality
uses kappa V>c and is not lost when U or individual discarded masses vanish.
5. With m=R-L, the admissible interval is1<=L<=R-18. The logarithm of
`(3/4)^(R-L)/L` is convex, so its maximum is at an endpoint. This maximizes
the quantity itself because logarithm is increasing. `Q_R=R+33/2+8/R`
is exactly the odd-order envelope for the upper bound on4V, and it also
bounds even order. Both endpoint expressions decrease forR>=263 by the
stated ratio and rational-excess arguments. Thus a finite endpoint check
atR=263 legitimately supplies the all-order tail bound.
6. In the local dual identity, lambda*d=t^2 and4lambda=t/2. Expanding the
right side gives exactly the left side. For0<=w<=t the remainder
`w(2t-w)/4` is nonnegative, with no unstated assumption onS>=0 beyond its
meaning as mass. Heavy groups use the separate selected-edge bound.
7. `dL-4U>0` forcesD>0 in the light-edge case. The strict budget bound then
producesW>hD>kappa D>=B>=W. With no light edges, D>0 produces the analogous
strict contradiction, so D=B=0. The accepted matching equality theorem
then yields the advertised fair-pair equality classification.

Accordingly, once the exact endpoint signs are independently confirmed, this
proof establishes sharp weighted stability for all n>=526. It makes no claim
for orders10..525 and uses no finite graph census or spectral analogy.
