# When finitely many precision sensors can survive a change of tree shape

2026-10-26. **PROVED:** the semantic closure criterion for a finite bank of
rational 2-adic proximity sensors; its common-center exception; and an
integer-cylinder family sharing the actual 86-macro history and all 67
specified entry sensors while regenerating unbounded periodic fuel.
**FINITE-EXACT:** the stated actual state, bank constants and bounded replays.
**OPEN:** a source-7/source-27 invariant controlling arbitrarily changing
tree shapes. The families below vary the original integer. They do not
prove failure, success or undecidability at the original integer 7.

## Inheritance and the new question

The closest proved mechanism is the signed center-switch identity in
[Session66, sections 2–6](../../session66_root27_aliquot_graph_interfaces_20261026/session66-root27/REPORT.md).
Its canonical actual hostile is a refill at native macro index 86: the
periodic fuel for `100` decreases from 5 to 2, then changing the reference
to `101000` gives fuel 9. The corrected near miss is that a finite bank
controlling repetitions of each shape need not control switching shapes.
The least-used sidecar is the **entire vector** of current center distances,
rather than one chosen center's distance. The board is: exact source / live
credit / affine carry / finite sensors / pulled centers / actual history.

The [independent center-bank appendix, Lemma A](../../session66_root27_aliquot_graph_interfaces_20261026/session66-audit-root27/CENTER_BANK_APPENDIX.md)
already proved a syntactic obstruction: a finite nonempty center set closed
under inverse block maps must be a common-fixed-center singleton. The new
step proves that this closure is **necessary even for nonlinear prediction
from the whole distance vector and any fixed congruence information**.
An omitted pulled center is not merely an inconvenience in one formula:
it produces a single observation class with unbounded successor distance.

This is a result about finite numbers of sensors, not finite-valued states:
each valuation itself can be arbitrarily large. No claimed map to tournament
clones, Frobenius defects or graph symmetry is needed for this obstruction.

## 1. Domains, legal blocks and what is observed

Use the inherited arithmetic macro

    T(z)=(3z+c)/2, c=z mod2,       q'=q+3c-1.

The machine sector has integer `z>=3,q>=1`. At `q=1,c=0` the macro is
blocked. A nonempty balanced word P has increment +2 per 1, -1 per 0,
nonnegative prefix increments and total zero. Write

    L=|P|, A=3^L, B=2^L, C=C_P,
    C_w=sum_(i<|w|) c_i 2^i 3^(|w|-1-i),
    F_P(z)=(Az+C)/B,     G_P(t)=(Bt-C)/A,
    r_P=-C/(A-B).

Its **exact cylinder** consists of integers satisfying `Az+C=0 mod2^L`.
These, and only these, have parity prefix P. Every such positive integer
source can execute P at every `q>=1`, and P leaves q unchanged. Restricting
to `z>=3` preserves every construction below. Positivity of all intermediate
states follows from literal iteration of T on positive integers; the native
guards follow from the balanced prefix inequalities.

Let S be a finite set of rational numbers of odd denominator, including
zero or negative numbers if desired. Put `h_s(z)=v_2(z-s)`, with
`v_2(0)=infinity`. For an arbitrary fixed positive integer M, observe

    O_(S,M)(z) = (z mod M, (h_s(z))_(s in S)).

The target is the vector of successor distances `(h_s(F_P(z)))_(s in S)`.
The extra residue is allowed as input information. We do **not** claim the
successor residue modulo M is determined by the same input: division by
`2^L` can require further binary digits. Initial q may also be fixed and
recorded; it does not change the theorem.

## 2. Semantic closure is equivalent to inverse-center closure

**Theorem 1.** Fix a nonempty balanced P, finite S and M>=1. On all
positive integer sources in P's exact cylinder (or all such sources >=3),
the successor distance vector is a function of `O_(S,M)` if and only if

    G_P(S) subseteq S.                                      (1)

If (1) fails, one observation class contains sources whose one fixed
successor sensor tends to arbitrarily large *finite* values. Consequently
no finite-valued upper bound on that successor sensor can be a function of
these observations. An identically infinite bound is of course vacuous.

**Proof, sufficiency.** Exact affine substitution gives

    h_s(F_P(z)) = h_(G_P(s))(z)-L,                          (2)

including infinity with the usual convention. Thus each required input
sensor is present under (1). The proof does not need the residue sidecar.

**Proof, necessity.** Suppose `s0 in S` but `rho=G_P(s0) notin S`. Then
`H=max_(s in S) v_2(rho-s)` is a finite nonnegative integer. Write
`rho=a/b` in lowest terms with b positive and odd, and `M=2^e m`, m odd.
For arbitrarily large integers `N>H`, `N>=max(L,e)`, impose

    bz-a = 2^N mod 2^(N+1),      z = t mod m.               (3)

Here t is any fixed residue modulo m. Since b is odd, the first condition
is one residue modulo `2^(N+1)`. The Chinese remainder theorem combines
it with the odd condition. Choose an arbitrarily large positive lift.
All resulting z have the **same compatible** residue modulo M: its
2-part is rho's residue modulo `2^e`, and its odd part is t. No claim is
made that an arbitrary independently prescribed incompatible binary
residue can be used.

Equation (3) gives `v_2(z-rho)=N`. For every s in S, the unequal-depth
ultrametric identity gives `h_s(z)=v_2(rho-s)`, independent of N. Also
`A rho+C=B s0`, so rho belongs to P's cylinder over the odd-denominator
rationals. Agreement through at least L binary digits puts z in the same
integer cylinder. P is therefore legal and (2) gives

    h_s0(F_P(z))=N-L.

This is unbounded within one observation class, proving necessity and
the stronger bound obstruction. Large positive lifts avoid equality with
any finitely many positive sensor centers. This also shows positivity
does not remove the obstruction. If S is empty, the target vector is
empty and (1) holds vacuously; the necessity argument is then unnecessary.

**Corollary 2 (complete finite-bank boundary).** Let Bcal be any nonempty
collection of nonempty balanced blocks. A finite **nonempty** sensor bank
can predict its entire successor distance vector under each legal block
exactly when it is a singleton `{r}`, with `r=r_P` for every P in Bcal.
Such a bank exists precisely when all the blocks are powers of one common
primitive bit word.

Indeed, Theorem 1 gives `G_P(S) subseteq S`. Injectivity and finite size
make this equality. For any P,

    diam(S) = diam(G_P(S)) = (2/3)^L diam(S)

as real sets, hence S is a singleton. Its member must be fixed by every
G_P. The common-center/commuting-word equivalence was proved in Session66,
section 3. Conversely a common-center singleton satisfies (2). The empty
bank is a separate vacuous case; an empty block collection would impose
no condition and is excluded here.

This supplies a positive, scoped rank: on a legal concatenation of
common-center blocks, the finite initial quantity `v_2(z-r)` decreases by
exactly each executed block length. Since r is negative and z is positive,
it is finite. Arbitrarily changing among powers is therefore paid by one
fuel account. Different lengths still have different maps and durations:
the center quotient does not permit their substitution as clocked macros.

## 3. A full observation vector at an actual root-27 state

Take

    P=100,      r_P=-9/19,
    Q=101000,   r_Q=-351/665,
    rho=G_P(r_Q)=-977/1995.

Here `v_3(rho)=-1`. Every balanced center has denominator coprime to 3,
because `3^L-2^L` is coprime to 3. Thus rho is absent from **every** bank
consisting solely of balanced periodic centers, even with zero included.
Inverse transport forces a genuinely different kind of center.

For a concrete finite bank S take zero and all centers of words obtained
by deleting the last leaf from full ordered ternary trees with 1–4
internal vertices. There are 71 words and 67 distinct sensors including
zero. Exact computation gives

    max_(s in S) v_2(rho-s)=11.

Literal native execution from `(z,q)=(7,1)` through **86 executed macros**
gives

    z*=10952205219944621, q*=47,
    F_P(z*)=36963692617313097,
    h_P(z*)=5, h_P(F_P(z*))=2, h_Q(F_P(z*))=9,
    v_2(z*-rho)=12.

Therefore the full 67-sensor vector at the actual state equals the vector
at rho. Its histogram is

    depth:  0  1  2  4  5  6  7  8  11
    count:  1 30 22  6  3  1  1  1   2.

The live credit 47 is not any of these precision values. Native macro
index 86 is the index of the single-bit arithmetic updates T, not the
index of odd returns or the number of P/Q blocks. These exact numerical
facts are bounded controls, not evidence of infinite continuation.

## 4. Exact-prefix hostiles with unbounded regenerated fuel

**Theorem 3.** For every m>=2 and every positive odd integer M0, there
exists a positive integer source `s_m>=3`, distinct from 7, with all of
the following properties.

1. Its first 86 executed native macros, from initial credit 1, have exactly
   the same parity word w and credit history as the source 7. Let its
   resulting entry be z_m, with credit 47.
2. `s_m=7 mod M0` and `v_2(s_m-7)=98` exactly. At the entry,
   `z_m=z* mod M0` and `v_2(z_m-z*)=12` exactly.
3. All 67 entry sensors equal the corresponding sensors at z*, while
   after the legal block P the Q-periodic fuel is `6m`.
4. Its next parity word is `P Q^m 0^47`. All but its final symbol execute;
   that final 0 is blocked. Its total executed-macro count is
   `86+3+6m+46`.

**Proof.** The rational rho has infinite parity word `P Q Q Q ...`:
P carries it to Q's fixed center, whose exact parity period is Q. Prescribe
the finite word `v=P Q^m 0^47`. At initial credit 47, its balanced blocks
keep the credit at 47, then its first 46 zeros reduce it to 1, and its
last zero is blocked. The exact cylinder for v is one residue modulo
`2^|v|`. Combine it by CRT with

    z_m = z* mod (3^86 M0).                               (4)

Choose a positive lift greater than z*. Since v and rho's itinerary
agree through `3+6m` bits and differ at the next bit (0 versus Q's initial
1), their difference has valuation exactly `3+6m>=15`. As
`v_2(z*-rho)=12`, the unequal-depth law gives
`v_2(z_m-z*)=12`, and every entry sensor is fixed because its depth from
rho is at most 11. Equation (2) yields Q-fuel `6m` after P.

Pull back through the **actual** prefix:

    s_m=(2^86 z_m-C_w)/3^86,
    s_m-7=(2^86/3^86)(z_m-z*).                            (5)

The congruence (4) makes this an ordinary integer. A sufficiently large
positive lift gives `s_m>=3`. Equation (5) gives the odd congruence and
the exact binary depth 98. The integer identity `F_w(s_m)=z_m` implies
the exact prefix cylinder for w, so every intermediate parity is right.
The original prefix's legal credit history is unchanged, and the
positive-integer map T keeps every intermediate integer positive.
This proves actual legality, not just rational endpoint agreement.
Appending v then proves the stopping count, with its last zero excluded.

The binary precision limit is essential. These entry sources agree with
z* modulo `2^12` but disagree modulo `2^13`; their initial sources agree
with 7 modulo `2^98` but disagree modulo `2^99`. Arbitrary odd sidecars
are preserved, not arbitrary binary precision at this fixed actual cut.
The sources vary and cannot be substituted for the single source 7.
They are not claimed to have integer ancestry root 27 (or root 7).

## 5. What remains, and the next useful certificate question

The source/target map here is the exact affine block F_P; the preserved
predicates are integer-cylinder membership, native guards, full finite
sensor values and the specified history. The observation quotient loses
distance to an omitted **pulled** center. Its cheapest decisive hostile is
the CRT family in Theorem 1; Theorem 3 puts a larger version next to a
known actual state, while marking its finite binary separation.

This rules out a universal exact refresh rule based only on a finite
rational sensor bank plus bounded congruence data, except for common-center
blocks. It does not rule out a root-specific invariant, a context-restricted
block language, an infinite bank with an independent bound, or a rank that
retains further arithmetic data. In particular, retaining the entire
integer source distinguishes all our hostiles. The prior full-bank lemma
says that maximizing over **all** balanced centers at initial credit 1
merely recodes the legal lifetime; that maximization cannot be assumed
as a free terminating oracle.

A sharper remaining target is therefore a **source-specific exclusion
certificate for newly pulled centers**, carrying exact binary depth and
the signed cancellation residue. It must say why the actual source 7
cannot follow the escaping center hierarchy at a dangerous guard; no
finite-bank total-transition theorem can provide that step. The present
result identifies the necessary missing coordinate, but does not bound it.

## Reproduction

From the worktree root:

```powershell
python 05-knowledge/results/rank_geometry_metagraph_20261026/rank/check.py
python -O 05-knowledge/results/rank_geometry_metagraph_20261026/rank/check.py
```

The standalone checker uses exact fractions, modular inverses and literal
native replay, with no prior producer imports. It tests arbitrary zero,
positive and negative sensors; exact infinity transport; distinct-center
failures; equal-center powers; the empty bank; arbitrary moduli including
one with 18 binary digits; the 67-sensor actual-state vector; and positive
integer pullbacks with m in `{2,3,4,8,16,32}`. It explicitly checks failure
of higher binary agreement and excludes the final blocked symbol from
the execution clock. The universal statements depend on the proofs above,
not on the tested m range. Normal and optimized receipts are saved beside
this file and should be byte-identical.
