# The sharp-defect boundary has zero bounded-start density

Status: RESERVED / proof candidate awaiting independent audit.

This addendum depends on the complete physical normal form in PROOF.md. It concerns **all** positive-start physical sources with D=s at a fixed start value, not just a selected seed orbit. It does not exclude general sources with D>s and does not disprove F(8,137).

## 1. Exact readout for the optional flipped arc

Use the normal form: U has new vertices v,a, edges v->a and T->v, and a dominates T except possibly x->a, where start_x(T)=end_x(T)=0. Write

    t=H(T), c=C_2(T), l=L_x(T), r=R_x(T), z=E_3,x(T).

Here C_k counts k labeled spanning paths, allowing empty bins. E_3,x has a designated nonempty path ending at x and two other labeled bins. For the flipped case, direct insertion gives

    h=t+c+l-r,  s=t,  e=d=t+l,
    b=2c-2r+2z,
    B+t=3c-3r+2z.                             (1)

To count h, compare with the unflipped cyclic construction. A path using a->x has v,a consecutive; deleting these vertices yields an ordered two-cover of T whose second path starts at x and whose first path is nonempty. There are r such paths because start_x(T)=0. A newly admitted Hamilton path using x->a must end at v: after a has been used, v has no possible successor. Its T portions form an L_x(T) cover, and the second portion is nonempty because end_x(T)=0. Thus there are l new paths. The baseline count t+c proves h.

For d=H(U-v), start at a to obtain t paths, since start_x(T)=0. Inserting a immediately after x gives l further paths. The end count e=d follows from end_a(U-v)=end_x(T)=0. The start count is t.

For b=C_2(U-v), begin with C_2(a->T)=2c. Flipping a->x removes covers in which a immediately precedes x, counted by 2r: choose the labeled bin containing a. It adds covers in which x immediately precedes a, counted by 2z: splitting that bin at a gives the designated path ending x, its optional suffix, and the other labeled bin. These operations are exact inverse deletions and preserve the bin labels. Therefore b=2c-2r+2z. At D=s the already proved identities L=b-d and B=h+L-s finish (1).

In the unflipped case the inherited exact formula is B+t=3C_2(T).

## 2. Separating the marked component

For a flipped source, decompose T=P->M->Q, where M is the strong component containing x. Both P and Q are nonempty, because x has zero start and end counts. Path-cover counts multiply across an ordinal sum: restrict each labeled path to each component, or concatenate its successive component pieces. Hence

    c=C_2(P) C_2(M) C_2(Q),
    r=H(P) R_x(M) C_2(Q),
    z=C_3(P) E_3,x(M) C_2(Q).                  (2)

For r, the designated x-start path cannot meet P, so P fills the ordinary bin. For z, the designated x-end path cannot meet Q, so Q fills its two ordinary bins. The remaining concatenations are unrestricted and labeled. Substituting (2) into (1) gives

    B+t=C_2(Q)[3 C_2(P) C_2(M)
               -3 H(P) R_x(M)+2 C_3(P) E_3,x(M)].  (3)

## 3. A finite family of two-parameter expressions

Fix a positive integer t. In any tournament with H(T)=t, the product of the Hamilton counts of its strong components is t. Every nontrivial strong component has order at most its Hamilton count, since each vertex starts a Hamilton path, and has Hamilton count at least three. Thus there are only finitely many possible ordered lists of nontrivial strong component types. Singleton components may occur arbitrarily often.

Fix such a list and, for the flipped case, fix the marked component M and mark x, or a position for a marked singleton. Delete the other singleton components, retaining the nontrivial components in P and Q. Let k and l be the numbers of deleted singleton components in P and Q. Since C_2(singleton)=2 and C_3(singleton)=3, (3) becomes

    B+t=2^l (alpha 2^k + beta 3^k - gamma),    (4)

for fixed positive integers alpha,beta,gamma depending only on the retained component list and mark. Here

    alpha=3 C_2(P_core) C_2(M) C_2(Q_core),
    beta=2 C_3(P_core) E_3,x(M) C_2(Q_core),
    gamma=3 H(P_core) R_x(M) C_2(Q_core).

Empty core cover counts are one. Since C_2(P_core)>=H(P_core) and C_2(M)>=2R_x(M), alpha>=2gamma. Also beta>0. Thus the right side of (4) is at least beta 2^l 3^k. Every allowable B<=X therefore has k,l=O(log(X+t)). Each retained list supplies at most O(log^2(X+t)) possible B values. Their number is finite for fixed t.

For unflipped sources, the same component decomposition gives B+t=alpha 2^k with finitely many positive alpha values. Combining both cases proves:

**Theorem.** For every fixed positive integer t,

    #{B<=X : some actual positive-start source has s=D=t}
        = O_t(log^2(X+2)).

For every fixed start cap K, the union over 1<=t<=K obeys the same logarithmic-square bound, with a constant depending on K. Its density among odd B is zero.

The conclusion is about the exact equality face D=s. Most targets in any eventual bounded-start profile construction must therefore use strict defect D>s. Finite shifts B=X-2j do not change this zero-density conclusion. This is an obstruction to using only minimal-defect sources; it is not a profile-spectrum gap or a failure of F(8,137).
