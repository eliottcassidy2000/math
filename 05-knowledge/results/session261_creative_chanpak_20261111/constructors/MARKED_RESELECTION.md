# A literal same-multiplier return invisible to every marked module context

**PROVED**, with explicit finite certificates. This supplies a real same-q
reselection outside the L/R/G carrier. It does **not** supply q→2q+2 or prove
universal rearming. The activated mixed ladder itself is inherited from
Session260; the new contribution is the context-preserving reselection and
its retained/lost-state description.

## 1. Complete marked module interface for H,W,U

Let (M,O) be a nonempty selected tournament. For each r>=1 let

- C_r count ordered nonempty spanning r-path covers;
- S_r count those whose first path starts in O;
- E_r count those whose first path ends in O;
- K_r, r>=2, count those whose first path ends in O and second path starts in O.

Put K_1=0 and all entries beyond |M| equal to zero. Call the four-sequence
array F(M,O) its marked cover tensor. Slot permutations make the choice of
first and second distinguished slots immaterial, as long as they are distinct.

**Marked-context theorem.** Replacing a selected module (M,O) by (M',O') with
the same tensor F preserves the host's (H,W,U), for every quotient tournament
and every fixed set of exterior selected vertices. This allows arbitrary
marks inside both the old and new module.

Proof: compress maximal module runs in each path of the host object. For H,
each skeleton uses C_r. For U, the required second-path start is either an
exterior selected vertex, giving C_r, or lies in a module run, giving S_r.
For W, its two required endpoints are both exterior (C_r), only the first
path's end is interior (E_r), only the second path's start is interior (S_r),
or both are interior (K_r). In the last case they belong to different global
paths and hence to distinct module runs. The ordering of the local run slots
is fixed by the skeleton and can be permuted to the distinguished slots in
the definitions. An empty first path is allowed only in U and is handled by
the same decomposition. The finite signed-free skeleton coefficients depend
only on the quotient, exterior orientations, marks, and path breaks. Equal
four-sequence tensors therefore give equal values of all three statistics.

The tensor intentionally does not preserve a distinguished vertex's identity
or its future substitution response. It is sufficient for all current marked
contexts, not for operations that look inside the module.

For a whole module viewed as a standalone seed,

    H=C_1,    W=K_2,    U=S_1+S_2.

The last identity accounts for the empty first path separately. These
conventions are checked directly by the supplied permutation certificate.

## 2. Actual same-q return on a strong five-vertex tournament

Take outgoing-bit rows

    T=(2,4,9,19,7),
    O={0,3},    O'={1,4}.

The two different selected sets have exactly the same marked tensor. Listing
(C_r,S_r,E_r,K_r) in increasing r gives

    r=1: (9,4,4,0)
    r=2: (66,28,28,8)
    r=3: (186,76,76,20)
    r=4: (240,96,96,24)
    r=5: (120,48,48,12).

Thus both are positive unit seeds with

    H=9, W=8, D=1, U=32.

The operation is merely to reselect the two vertices, leaving the tournament
unchanged. It is a literal exact-q return, and the marked-context theorem
shows that no surrounding module context measured by H,W,U can distinguish
the two presentations. No postulated state label is assigned.

Yet their future operation menus differ. Under O, vertex0 is selected. Clone
it to length m and retain only its first clone as selected. At m=1,2,3 the
exact triples (H,D,U) are

    (9,1,32), (15,1,86), (27,3,230).

The selected-first continuation therefore fails to stay a unit seed at m=3.
After the same-q reselection to O', vertex0 is unselected. Its entire
transitive-clone ladder has

    H_m=3*2^m+3,
    D_m=1,
    U_m=4*3^m+9*2^m+2,          m>=1.

The first multipliers are32,74,182,470,1262,... . This is exactly the mixed
family in Session260 `constructors/RESEARCH.md` / `module_ladders.py`, now in
an explicit vertex gauge that exhibits a context-invisible seed return.
The formula is credited there, not counted as a new ladder.

The certificate records its finite skeleton coefficients, including the
separate deletion boundary m=0:

    H: (5,9,3,0,...),
    D: (-1,1,0,0,...),
    U: (15,32,21,4,0,...).

These prove the all-length formula in the surjection basis. A finite table
of m-values alone would not suffice.

## 3. Infinite physical transport and the missing digit

Appending b unselected exterior sink vertices transports the exact same-q
return to q=2^(b+5), for every b>=0. The usual sink-padding law gives the
post-return family

    U_(m,b)=2^b(4*3^m+9*2^m+2),  D=1.

Every such multiplier has actual tail counts1+2^t U_(m,b), t>=1. This is a
constructive infinite family of returns and descendants, not a bounded
same-q coincidence. The scalar recurrence for the m-direction needs both
2^m and3^m; it does not become a binary digit recurrence by forgetting one.

In particular the first activated birth is32→74, not32→66. Thus the binary
rearming criterion of Session260 remains unproved. The useful change is
structural: exact return operations exist without preserving the future port
menu, even when they preserve the complete H/W/U module-response tensor.
A search for a universal rearm should therefore retain root-mark incidence
and response jets, rather than quotienting by the tensor alone.

## 4. Connection to the first module wave

The unmarked cover content of T is gcd(9,66,186,240,120)=3. The theorem in
MODULE_REBASE.md excludes T from being a wholly **unselected** module in any
unit seed. Nevertheless the internal selected response has S_1=E_1=4, and
the marked tensor's content is1. Carrying the marks through the interface
literally removes the earlier divisibility obstruction. The distinction is
not semantic: T is an actual unit seed with the two markings above.

Transfer ledger: source = boundary-count tensor / target = seed state return /
map = same-host reselection O→O' / preserved = all marked-context H,W,U /
lost = selected status and rooted continuation menu / sidecar = a distinguished
port and its clone response / decisive test = exact tensor equality plus the
first continuation that separates the menus. All are supplied here.

## Reproduction

    python3 marked_reselection_certificate.py

The tensor equality is enumerated literally over vertex permutations and cut
sets. The activated ladder is separately checked in the finite skeleton basis
with the inherited Session260 helper. Output: `marked_reselection_receipt.json`.
The result does not depend on the wider exploratory census.

## 5. Mechanism recovered: balanced reversal, not an accidental collision

The involution rho(i)=4−i is an anti-automorphism of T and carries O to O'.
There is no automorphism carrying O to O' (the certificate checks the120
permutations). Reversing every local path sends

    (C_r,S_r,E_r,K_r) to (C_r,E_r,S_r,K_r).

For K_r, reverse the two distinguished paths and then exchange their slots.
Thus **any** marked module with S_r=E_r for every r admits a full-tensor
preserving reversal. When an anti-automorphism is available, it can be
implemented as a reselection on the same underlying tournament. The displayed
five-vertex state is a concrete balanced example. This explains its tensor
equality conceptually rather than leaving it as a numerical coincidence.

This also identifies the exact second response needed by a simpler reversal
return. Define V=E_1+E_2. Full reversal always preserves D and interchanges U,V;
it is a same-q return precisely when U=V. Equality of every S_r,E_r is stronger:
it ensures the return works in all surrounding module contexts as well.

On the old marking O, cloning its unselected vertex4 gives the two responses

    U_m=16*2^m,
    V_m=4*3^m+9*2^m+2,             m>=1.

The second formula follows by the anti-automorphism from the activated port0;
the first follows from the literal clone skeleton. They coincide at m=1.
For m>=2, U_m is a power of two while V_m is2 modulo4 and greater than2.
Therefore no amount of exterior source/sink dyadic padding can make these two
responses equal again. This is an exact failure of repeatable balancing on
this ladder, not a failure of the first physical return. A repeated reversal
route needs a new way to balance both endpoint responses.
