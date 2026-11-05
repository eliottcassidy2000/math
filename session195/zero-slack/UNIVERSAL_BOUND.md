# Universal bounded-deletion classification by finite cores and dyadic padding

Session 195. **PROVED refinement of the zero-slack packet.** The finite-bank hypothesis is unnecessary for the bounded-deletion conclusion, and the counting bound improves from `O(log² X)` to `O(log X)`. The original grammar-specific proof remains valid but is weaker than this result.

Use the conventions of `PROOF.md`: `(h,e,R,d)` is a rooted tournament state and `f=R+e-h-d`. Its Theorem 2 proves universally that f=0 exactly when the root has outdegree at most one.

## 1. Exact SCC transport

For a nonempty tournament S, write its ordered strong components as `C1 -> ... -> Ck`, and let a lie in Cj. Put

    Q(C)=sum_(U subseteq C) H(C[U])H(C-U).

Then

    R_a(S) = (product_(i<j) H(Ci)) R_a(Cj) (product_(i>j) Q(Ci)),
    e_a(S) = [j=k] (product_(i<j) H(Ci)) e_a(Cj).       (1)

To prove the first formula, a path starting at a cannot visit an earlier component. Its intersection with Cj is a nonempty path starting at a; its intersection with each later component is a possibly empty path. All these pieces concatenate in SCC order, and their unused vertices contribute independent Hamiltonian paths in the same SCC order. The decomposition and concatenation are inverse. The second formula is ordinary SCC multiplication, with the endpoint necessarily in the last component. In particular, `H(K1)=1` and `Q(K1)=2`.

## 2. A finite-core normal form at every bounded deletion count

**Theorem.** Fix D>=1. For every rooted tournament (T,r) with d<=D and f=0, either r is a universal sink and h=d, or there is a rooted tournament (C,a), of order at most D+1, and an integer t>=0 such that

    H(C)=d,
    h=d+2^t R_a(C)-[t=0]e_a(C).                       (2)

C has no singleton strong component except possibly the one containing a. Conversely, each such rooted core C and t constructs a physical zero-slack tournament of deletion count d with the count in (2).

**Proof.** If r has no out-neighbor, it is a universal sink. Otherwise Theorem 2 of `PROOF.md` says it has exactly one out-neighbor a. Set S=T-r, so H(S)=d. Delete from S every singleton strong component except the one containing a, if that component is a singleton. Let C be the resulting induced rooted tournament, and let t be the number of deleted singleton components after the component containing a.

Earlier singleton components contribute factors one to (1); later singleton components contribute factors two to R. Deleting singleton components does not change the nontrivial components, their order, or H. Formula (1) gives exactly

    H(C)=H(S)=d,
    R_a(S)=2^t R_a(C),
    e_a(S)=[t=0]e_a(C).                              (3)

If a still has later nontrivial components, both endpoint terms in (3) are already zero; otherwise t>0 is precisely what makes the old endpoint impossible.

Each nontrivial strong component B has at least one Hamiltonian endpoint path at every vertex. Thus `|B|<=H(B)` and `H(B)>=3`. If their Hamilton counts are h1,...,hm, their product is d. For m>=2, the elementary inequality `sum hi<=product hi` holds when every hi>=3 (start with x+y<=xy and induct). For m=1 it is equality, and for m=0 there are no nontrivial vertices. Therefore the total number of vertices in the nontrivial components is at most d, and retaining one marked singleton gives `|C|<=d+1<=D+1`.

The exact singleton-birth formula applied to r is `h=d+R_a(S)-e_a(S)`. Substitution of (3) proves (2).

Conversely, append t universal sink vertices to C as a transitive chain, retaining the mark a. Then adjoin r with a as its sole out-neighbor. The same formulas prove (2), d=H(C), and f=0. QED.

This is an exact parameterization, not only an upper envelope. Singleton components before the marked component may be deleted. Those after the marked component may be moved to one terminal chain. These specific operations preserve H(S), R_a(S), and e_a(S).

## 3. Consequences and scope

There are finitely many rooted tournaments of order at most D+1. Each contributes only one exceptional t=0 value and a dyadic family `d+C0*2^t`, t>=1, with C0=R_a(C)>0. Hence the number of distinct H values at most X representable by **any tournament** with d<=D and f=0 is

    O_D(1+log X).                                    (4)

The unit-deletion classification in the original proof is the case C=K1: H is 1 or `2^t+1` for t>=1.

Meanwhile, the physical A,D macros in the original proof attain every `3^a 5^b 17^c`, giving `Omega(log³ X)` distinct counts up to X. Consequently, for every fixed D, infinitely many already-attained counts admit **no rooted tournament representation whatsoever** with both d<=D and f=0. This strengthens the earlier grammar-only conclusion and allows arbitrary gadget banks or arbitrary reconstruction within that proposed normalization.

This is still not a counterexample to the full Hamilton-path spectrum or to the enlarged grammar's all-order scalar coverage. It rules out a bounded-deletion, zero-slack normalization architecture. Counts may use unbounded deletion data or a different marked interface.

## 4. Reproduction

Run `python3 verify_universal.py` and `python3 -O verify_universal.py`, comparing stdout. The verifier checks (1)-(3) and the finite-core bound on every rooted labelled tournament through order five, plus larger padded controls. It imports only the standalone `verify.py` included in this same packet. The all-order theorem follows from the SCC argument above, not the finite controls.
