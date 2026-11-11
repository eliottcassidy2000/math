# Re-marking turns the cofactor-21 trap into a reusable valuation controller

**PROVED application of the inherited physical fan/graft law, with new
literal packet and all-level valuation classification. Not exact target coverage.**

## Why this extension was tried

The previous note exhibits an actual H609 graph whose gadget-content quotient
is the forbidden 21. Rather than divide by 29, keep all nine vertices and
move the marked interface into the gadget. This asks whether *physical*
continuation, not rational normalization, can escape the content.

A reproducible unrelated source draw used Python `random.Random(27521)` on
sorted tracked `05-knowledge/results` Markdown paths containing dates
202603 through202608. One selected source was
`spectral-H-and-the-LRC-floor-duality-opus-20260629.md`.
Its opening comparison between a fixed second moment and a varying
concentration statistic suggested looking beyond the two aggregate
coordinates H,Q. Only that methodological contrast is used; this note does
not rely on its broader maximizer equivalences or historical open claims.
The concrete replacement statistic here is the rooted response U−E.

## 1. The actual escape, with no division

Keep the graph from PROOF.md,

    T=(2,508,433,165,9,17,189,49,249), H(T)=609.

At its vertex7, direct exact counting gives

    (D,E,U,d)=(609,108,1066,241), f=U−D−d+E=324.

A singleton fan birth b→7, T−7→b has count

    D+U−E=1567 ≡ 1 (mod29).

The full root menu is in `escape_receipt.json`. At the two outer triangle
vertices, the same operation gives1189 and1711, both divisible by29.
At each of the seven internal gadget vertices it escapes29. Thus the loss
was not an unavoidable property of the scalar609; the choice of physical
root changes the available continuation. No assertion of minimum size or
novelty of the scalar1567 is made.

## 2. Full retained-root packet at every fan length

For any nonsink-rooted source with the notation of PROOF.md, replace the
singleton birth by a transitive k-vertex fan, k>=1. Its gadget packet is
(H,Q)=(1,2^k), so the inherited graft identity gives

    D_k = (D−E+f)2^k+d−f.                         (1)

The root is still nonsink. Deleting it leaves the ordinal join of A−a and
the fan, so d_k=d. All Hamilton paths ending at the root consist of a path
through A−a followed by the fan and then the root, so E_k=d. Adding another
ell>=1 fan is the same physical graph as the (k+ell)-fan. Comparing (1)
with the proper-graft expression from the new source, for two distinct
values of ell, gives f_k=f. Therefore the **full four-coordinate packet** is

    (D_k,E_k,U_k,d_k)=(D_k,d,D_k+f,d).             (2)

This derivation makes the recurrence's boundary explicit: k=0 need not
obey either (1) or (2). It is the same physical fan mechanism as
Session274, not a new generic fan theorem.

For the concrete marked T,

    (D_k,E_k,U_k,d_k)=
      (825·2^k−83, 241, 825·2^k+241, 241), k>=1. (3)

Thus every subsequent single birth has the reusable exact law

    D_(k+1)=2D_k+83,

while d=241, E=241, f=324 remain constant. These are real tournaments,
not formal matrix states. The original T is strong; adding a vertex that
beats its nonsink root and loses to the remaining graph preserves
strong connectivity here, and every fan has the same direct reachability.

## 3. Exact all-level valuation control

Modulo29, (3) reads

    D_k = 13·2^k+4.

The order of2 modulo29 is28. Exactly one class of exponents is bad:

    29 | D_k  iff  k ≡ 26 (mod28).                (4)

In particular every k≡1 (mod28) gives D_k≡1 (mod29), an infinite
renewable family that has escaped the original content29.

There is a sharper all-level statement. For each r>=1 there is exactly
one residue class k_r modulo28·29^(r−1) such that

    29^r | D_k  iff  k ≡ k_r mod 28·29^(r−1).     (5)

Proof: (4) is the base. Direct arithmetic gives
v_29(2^28−1)=1. The binomial theorem then gives
v_29(2^(28·29^(r−1))−1)=r. For P=28·29^(r−1),

    D_(k+jP) = D_k + 825·2^k·j(2^P−1) mod29^(r+1).

After dividing by29^r, its coefficient ofj is nonzero modulo29. Exactly
one of the29 lifts of k_r therefore survives at the next level. Conversely
the order of2 modulo29^r is28·29^(r−1), by the same binomial calculation,
so these are all solutions. This proves (5).

For every r>=0 there are consequently infinitely many actual fan outputs
with v_29(D_k)=r: take a class surviving levelr but not its unique lift to
levelr+1, with r=0 interpreted using (4). The reusable physical interface
controls arbitrary valuation depths while the finite-height source is
fixed. The first four bad classes are recorded in the receipt.

## 4. Productive conclusion and exact remaining gap

The normalized21 fibre does **not** imprison every continuation in its
content. Re-marking supplies a concrete coprime escape, and the retained
interface is renewable at all lengths, with an explicit complete packet
and exact valuation schedule. This is stronger than a single accidental
gcd-one output or a finite modular experiment.

It still does not select arbitrary odd totals. For the fixed interface (3)
the output set is the single exponential progression 825·2^k−83. Its
renewal is rigid. To turn it into a Chan–Pak constructor one must add a
physically available control that changes the intercept/response while
preserving enough of (2) to invert prescribed absolute targets. Dividing
counts, taking differences, or asserting that arbitrary29-adic depths
fill integers would repeat the normalization error this example exposes.

## Reproduce

    python3 05-knowledge/results/session275_exact_selection_20261111/descent/escape.py

Controls compute all nine initial root packets, independently construct
and recount their singleton fans, and verify the full packet (3) for
k=1,...,5 by subset dynamic programming on the actual arcs. Four arithmetic
lifting levels are checked as controls of the all-level proof, not as its
basis. Normal and optimized Python outputs agree. Required checks use
explicit exceptions, not removable assertions.
