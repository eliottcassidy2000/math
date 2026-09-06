# A retained cancellation factor converts a collective gcd into a pair gcd

**Status: PROVED ELEMENTARILY + INDEPENDENTLY AUDITED; FINITE-EXACT controls.** This is a precise conditional bridge. It does not prove
that every remaining actual eleven-component has a qualifying path and
does not close LRC(14).

The inherited suppliers are the audited twelfth signed-box/gcd theorem,
`05-knowledge/results/overnight12_20260906_lrc_gcd_semigroup.md`, and the
audited connected decoder descent,
`05-knowledge/results/overnight12_20260906_lrc_decoder_descent.md`.
They retain the actual THM-3818 graph, the physical box, and W=V_dec. The
current joint-shadow caps are30,90 at sizes8,7. Its namespace reservation
is not treated as a proved dependency. The completed proof note is
`05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.md`.

The attempted direct implication from a five-coordinate gcd to a pair
gcd loses support information. A walk through actual decoder edges restores
an exact sidecar: the cancellation between its accumulated numerator and
denominator. The live board is collective gcds, actual edge orientations,
prefix denominators, returning primes, endpoint pairs, and bounded crossing
certificates. No multicoordinate Bezout identity is turned into a relation
on three labels.

Let w_0,...,w_m be any positive integer walk, repetitions allowed. Write

    w_i/w_(i+1)=a_i/b_i in lowest terms,
    A=product a_i, B=product b_i, C=gcd(A,B),
    L=lcm of the reduced denominators of w_i/w_0, 0<=i<=m,
    d=gcd(w_0,...,w_m), D=gcd(w_0,w_m).

**Theorem.** The exact packet is

    L=w_0/d,
    J=L/(A/C),
    D=dJ,                       J divides C.          (1)

**Proof.** L is the smallest positive integer for which every Lw_i/w_0
is integral. This integer vector is primitive: a common divisor greater
than one would permit a smaller multiplier, while its first coordinate
is L. Hence it is the original vector divided by d, proving L=w_0/d.
The endpoint ratio is w_m/w_0=B/A=(B/C)/(A/C). Its reduced denominator
A/C divides L, and coprimality of its numerator and denominator gives
D/d=L/(A/C). Every prefix denominator divides its prefix product of a_i,
which divides A. Their lcm therefore divides A. Since A/C divides L,
the integer J divides C. This proves (1) without factoring any speed.

The weaker certificate C is often easier to bound, but J is the exact
lost factor. There is no claim that the entire cancellation C is realized
as loss of the collective gcd.

For a hypothetical unsafe actual 11+2 entry, let its physical core be tV,
with gcd(V)=1, K=max V. The current inherited hierarchy gives t<=2. Every
subset of the physical eleven-core of size5,6,7 has gcd at most

    M5=11,342,250, M6=31,950, M7=90.

Take an actual decoder walk starting at the physical maximum tK and
ending at another core vertex, visiting exactly r distinct vertices for
one of r=5,6,7. The gcd of the walk is the gcd of that subset, even if
vertices repeat. Equation(1) gives

    gcd(K,w_m/t) = D/t <= M_r J/t.                    (2)

Consequently the actual branch is safe if

    J <= floor(tH/M_r),             H=76,388,115.     (3)

Indeed (2) produces exactly the endpoint pair required by the proved
twelfth cutoff, with its actual coefficient budget. In a proof by
contradiction this violates the supposed strict failure. The thresholds
for r=5,6,7 are respectively

    t=1: 6, 2390, 848756;
    t=2: 13, 4781, 1697513.

The same sufficient test may replace J by C because J divides C. A
coprime accumulated numerator and denominator has C=J=1. Thus a suitable
walk with no returning prime immediately closes this branch. More generally
the thresholds specify the allowed exact cancellation, not a claim that
such a walk always exists. Arbitrary endpoints can instead be passed to
the full native two-coordinate radius criterion; retaining the global
maximum is what makes (3) a direct endpoint-cutoff consumer.

The cheap hostile is the actual atlas path 6→2→3. Its edge ratios are
3:1 and2:3, whose sums4 and5 are admissible. Its collective gcd is one,
but the endpoint gcd is three: A=6, B=3, C=3, L=6, J=3. Dropping J is
therefore false even on an actual connected decoder path. This is a
minimal hostile for this endpoint-versus-collective claim with distinct
endpoints, requiring at least three distinct vertices,
not a counterexample to the existence of some other small-gcd pair.
A positive control is729→243→81→27→9→3→1, seven distinct vertices with
actual1:3 atlas edges and J=1.

The three-distinct-vertex actual path6→4→1 has J=1<C=2, with admissible
edge sums5,5. A stronger five-vertex example is18→12→3→9→6: J=2 while
C=12. It passes the t=1 five-subset threshold J<=6, whereas the sufficient
replacement C<=6 fails. These vertices lie in the connected primitive core
V=(1,2,3,4,5,6,7,9,10,12,18), with18 its full maximum. An actual spanning
tree uses edges1–3,1–4,1–9,1–10,2–3,4–6,3–7,5–12,4–12,6–18, whose
primitive sums belong to{4,5,10,11,17}. Adjoining g,3g with g=36Q+1 gives
an actual finite-box equality entry by the same dominance argument as the
twelfth controls. This is a strict gain for the exact packet over its
cancellation-content bound, not an example of an unsafe row.

The bounded probe also investigated primorial-style five-coordinate sets.
They can have a small collective gcd and much larger pair gcds, but making
their enclosing eleven-graph connected can introduce extra vertices whose
seven/eight-subset gcds violate the current90/30 cuts. No such rejected
construction is asserted to be a remaining actual core. The general
implication sought by the parent remains unresolved: connectedness and
the inherited scalar/profile constraints have not yet been shown to force
a path satisfying(3), nor has a full hostile passing all of them been found.
The concrete progress is the exact factor that must now be bounded or
retained, and the resulting proved conditional endpoint criterion.

The companion standard-library source tests all7776 positive five-entry
integer walks over{1,...,6}, comparing the packet with direct gcds. It also
checks both named controls and the exact integer thresholds. Universal
claims follow from the proof; the finite bank tests the normalization.
The [standalone source](../../04-computation/overnight14_20260906_lrc_endpoint_walk.py)
and [frozen output](overnight14_20260906_lrc_endpoint_walk.out)
are reproduced by:

    python -B 04-computation/overnight14_20260906_lrc_endpoint_walk.py
    python -B -O 04-computation/overnight14_20260906_lrc_endpoint_walk.py

Both runs pass31,122 always-active gates with byte-identical LF output.
Source SHA256:
`e38ce143de236a35d089f0e7d4110fae5200ea819868e7b925a705f50428664a`.
Output and optimized-output SHA256:
`4178d510f86a4cdb4cca86afcb6405439e5179cddd924e48fffd1654a1729561`.

The [independent audit](overnight14_20260906_lrc_endpoint_walk_audit.md)
accepts the exact cancellation identity and its conditional LRC consumer.

**Filing:** root read these proofs and audits and integrated the fourteenth
checkpoint. Reproduction commands are relative to the repository root.
