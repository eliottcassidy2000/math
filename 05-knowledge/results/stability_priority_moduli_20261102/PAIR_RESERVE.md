# A local joint-profile reserve for the remaining small orders

Status: **PROVED**, elementary triangle identity. This is a possible input
for orders 10..275, not a proof of their sharp stability coefficients.
It is not needed by `STABILITY_N276.md`.

Let a weighted triangle have p=A_uv, a=A_vx, b=A_ux, q(t)=t(1-t), and
alpha=min(p,1-p). Its expected cyclic indicator is

    tau_uvx = p a(1-b)+(1-p)(1-a)b.

Therefore

    tau_uvx >= alpha [a(1-b)+(1-a)b]
             = alpha [q(a)+q(b)+(a-b)^2].                 (J1)

For a fair selected edge p=1/2 this is an identity. For arbitrary p the
coefficient alpha is best possible in an inequality
tau>=c(p)[q(a)+q(b)] valid for all a,b: if p<=1/2, set b=0 and let a
decrease to zero; if p>=1/2, interchange a,b. The endpoint cases alpha=0
are included. In terms of the selected-edge variance w, alpha is
(1-sqrt(1-4w))/2, so alpha>=w and alpha=2w for a fair edge.

For any matching M of selected edges, each triangle contains at most one
selected edge. Hence, without overlap,

    tau(A) >= sum_({u,v} in M) min(A_uv,A_vu)
                    sum_(x not u,v)
                    [q(A_vx)+q(A_ux)+(A_vx-A_ux)^2].      (J2)

Together with the inherited B=E Delta(T)>=tau(A), this supplies a lower
bound on the binary baseline retaining *joint external profiles*. The
unlabelled multiset of edge variances cannot recover the last summand:
with p=1/2, the choices (a,b)=(0,0) and (1,0) have exactly the same three
edge variances, but cyclic probabilities 0 and 1/2. The missing coordinate
is alignment of the two endpoint rows at the *same* external vertex.

This identifies a concrete reason that variance-only certificates lose
information. Identically oriented deterministic ports around a fair pair
make the reserve zero; discordant deterministic ports cost 1/2 per triangle.
That is compatible with the fair ordinal-pair equality locus. It is not a
classification of general strong-atom attachments, nor a way to drop the
actual complemented relation or the twelve-vertex 81/27 hostile control.

The next test is to combine J2 with the existing five-vertex sharp feasible
region and disjoint cycle terms on shared ports. Simply adding two overlapping
triangle budgets would be invalid. The accepted n=7 graph-only hostile
C4 disjoint-union K3 and the Pfaffian coefficient-3 cancellation remain
mandatory controls for any resulting small-order certificate.
