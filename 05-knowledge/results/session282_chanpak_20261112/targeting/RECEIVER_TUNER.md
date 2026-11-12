# A physical scalar-neutral receiver tuner, and its unavoidable content cost

**Session282. PROVED, independently reviewed.** This addendum uses the exact
source-root twins of [PROOF§6](PROOF.md), not the erroneous root-deletion/full
cover substitution corrected there. Chan–Pak remains OPEN.

Write their common H=h as a, and their Q/U values as q_0,q_1, with
δ=q_1−q_0≠0. For every m>=1, a=9H(T((000111)^m)), so9 divides a.
The full source-root carriers are

    X_i:(H,h,Q,U,V,d,e)=(a,a,q_i,q_i,a,0,a).

Their converses are actual sink-root objects

    Y_i:(H,h,Q,U,V,d,e)=(a,a,q_i,a,q_i,a,0).           (1)

Taking the converse exchanges U,V and d,e and leaves H,h,Q unchanged. This
is path reversal, not an identification of arbitrary scalar twins.

## 1. Positive receiver tuning at a fixed physical count

For any actual pointed tournament B, left fusion Y_i*B has

    H'=a H_B,  h'=a h_B,  V'=q_i V_B,  e'=0.          (2)

Thus the two choices have exactly the same scalar H and deletion h but
receiver difference δV_B. This is a physical operation, not subtracting
path counts. They also have the same zero e. Formula(2) follows immediately
from the inherited rooted-fusion matrix product; equivalently Y_i is a sink
root, so its deletion precedes the entire B in an ordinal sum. Scalar count
and deletion therefore multiply, while the root-ending two-cover mass uses
q_i instead of a.

At m=1, a=9, q_0=21504, q_1=18432, so the receiver change is−3072 V_B.
This can distinguish equal-H sources or tune a composite receiver. It is
not a +2 gate, nor a repair of the initialized e=h corridor without further
operations. Applying the known C1 preparation changes H by2V, so it exposes
rather than erases this difference.

## 2. Exact neutrality classification in every common fusion sandwich

Let L,R be arbitrary actual pointed tournaments, singleton included. Their
h,U,V are strictly positive. The inherited endpoint-cover injection gives

    U_L=h_L iff L's root is a sink;
    V_R=h_R iff R's root is a source.                 (3)

For nonextremal roots, the corresponding difference is at least h. These
are physical facts; they are not assumptions on freely chosen integer matrices.

For the source-root choice the matrix difference is
δ (e_1+e_2)e_2^T, and for the sink-root converse choice it is
δ e_2(e_2+e_4)^T. Therefore exact common-context differences are

 ΔH(L*X_i*R)=δ U_L(V_R−h_R),
 ΔV(L*X_i*R)=δ Q_L(V_R−h_R);                          (4)

 ΔH(L*Y_i*R)=δ (U_L−h_L)V_R,
 ΔV(L*Y_i*R)=δ Q_L V_R.                              (5)

In both cases Δh=Δe=0. Matrix indices here are one-based. The singleton
matrix is the identity, so the formulas do not introduce an empty-graph
endpoint convention.

Since δ≠0, equations(3)–(5) prove:

* Source-twin common sandwiches preserve H iff R has a source root. In
  that case the receiver is also unchanged: they do not supply H-neutral
  receiver repair.
* Converse-twin common sandwiches preserve H iff L has a sink root. In
  that case ΔV≠0, so a genuine scalar-neutral receiver tuner survives.

Repeated common left/right fusions collapse to such a sandwich by associativity.
Root changes, choice-dependent contexts and internal surgery are not classified.

## 3. The unavoidable arithmetic price of this tuner

If R has a source root, the whole source-twin sandwich is the ordinal sum
L -> (X_i−root) -> (R−root). Consequently

    H(L*X_i*R)=H_L a h_R.                             (6)

If L has a sink root, the converse-twin sandwich is
(L−root) -> (Y_i−root) -> R and

    H(L*Y_i*R)=h_L a H_R.                             (7)

Thus **every H-neutral use of these pairs through common rooted-fusion
contexts produces a count divisible by9**. This is an exact family-specific
content cost, not a theorem excluding arbitrary scalar-neutral gadgets.
It explains why the positive receiver tuner cannot simply repair the prime
or unit-mesh Chan–Pak corridor. Combining it with the h169 central gate in
such a common fusion still retains the factor9; a later root change or
non-fusion operation would need a new analysis before claiming escape.

The mechanism parallels the inherited reset warning of Session246: an exact
physical reset/tuner can carry an unavoidable multiplicative bill. Here the
bill is proved for the precise two-sided sandwich family, with the full
neutrality iff, not for every rank-one ideal or every tournament operation.

`verify_tuner.py` reuses the standalone literal helpers only to construct
actual graphs and endpoint profiles. It checks both converse/source cases,
the singleton boundary, source/sink/nonextremal contexts, formulas(4)–(7),
and a direct instance of(2). No new numerical coverage claim is made.
