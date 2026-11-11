# Independent analytic audit of the all-k triangle matching theorem

**PASS.** The proposed argument in
[TRIANGLE_ALL_K_PROOF.md](../arithmetic/TRIANGLE_ALL_K_PROOF.md) establishes
the claimed all-k extension-content formula for its full declared triangle-
chain family. It does not establish specified target realization or a
general Chan--Pak construction. This audit is independent of the producer's
finite two-birth coefficient certificates.

## Factorization and local quantities

The coefficient grouping can be refined to retain the selected/unselected
flag of each used edge factor, its old-path role, its block, and its selected
endpoint label. The expansion sign and new-new factors are then constant.
Unselected endpoint labels remain summed. For each old segment, its start
and end blocks determine its availability interval. Every compatible local
cover concatenates uniquely along the old ordinal order, so the grouped
count factors. This refinement changes no valuation argument.

If a nonzero group has an availability-zero region, that region contains
no old vertex of either type. A free triangle there, just like a free
singleton, would make the local count zero. Treating zero regions as
possibly empty, rather than assuming they contain a positive gap, is
therefore necessary and handled correctly.

The free rotation of an unmarked triangle preserves every grouped datum.
A triangle with at least two contacts needs only the elementary bound
loss<=1, paid by its contact budget. This covers contacts at the same
vertex, starts and ends of one old path in the same triangle, and multiple
old segments; no multiplication of separate one-event factors is used.

For a triangle with exactly one selected contact, all other available bins
are optional within that block. The three possible sizes of the fixed-
endpoint bin give `(q-1)^2`, `q-1`, and1, so the local count is
q²-q+1. Its valuation is exactly one for q=2 mod3 and zero otherwise.
An unprotected single contact therefore has exactly one adjacent region
whose availability is divisible by3. This uniqueness is the key invariant.

## The three external payment rules

A contact-free singleton at positive availability q contributes q.
When q is divisible by3 its one valuation unit pays the at most two
half-requests at that region's ends. Distinct regions use distinct
singletons, preventing duplicate expenditure.

A singleton with one contact has availability changing by one, so only
one adjacent region can be divisible by3. It can receive at most one
half-charge. With two contacts its one old vertex is a complete old path
segment: one start and one end, with no net availability change. There are
at most two adjacent requests, paid by its full unit of contact budget.
Allowing a singleton to carry contacts of two distinct old paths would be
illegal, but the stated upper bound would not become smaller; the proof
requires only the legitimate two-contact case.

At an exterior boundary, q is0 or1 according as the full Hamilton path
starts/ends new or old. A terminal request can use only q0 and corresponds
to a new full-path endpoint. Thus the two possible terminal half-charges
are paid by precisely the endpoint slack in
`m+2e_NN=2k-e_end`. No unjustified availability slack is introduced.

## The odd-run saturation step

After those payments, a request's entire region has triangle endpoints
and no singleton. Consequently it lies within one maximal zero-gap run.
Let B be the number of unpaid requesting triangles in a run of length L.
Every triangle requests at most once, so B<=L. For even L this suffices.

For odd L, the only problematic case is B=L. Then every triangle in the
run has exactly one contact and an unpaid internal request. No intermediate
contact-free triangle remains. Each requested region is between two
consecutive triangles. Its availability is a multiple of3, and the
unique-side property forces BOTH endpoints to request into that same
region. Thus these request edges are disjoint pairs covering all L
triangles, impossible when L is odd. Therefore B<=L-1. This argument also
covers a length-one run and makes no assumption about event order inside
a multiple-contact triangle: such a triangle simply prevents B=L.

The remaining half-charges are consequently at most
`sum floor(L/2)=nu`. Combining this with contact and endpoint budgets
gives v3(group)>=r-k-nu. The earlier untouched-triangle bound supplies
r-2k, yielding exactly `max(r-k-min(k,nu),0)`.

The reverse divisibility uses disjoint H41 pair extensions, H5 single
extensions and sink padding. An all-sink extension has H=3^r, so no other
prime can divide the full gcd. Since nu<=floor(r/2), the first k at which
the displayed exponent vanishes is r-nu. Small r and k0 are consistent:
unused births are sinks, and without births the gcd is exactly3^r.

## Independent hostile controls

The legal walk SSESESEE has availability2,3,3,3,3,3,3,2 at its eight
single events. Marked triangle factors alone provide only two factors3,
which would be insufficient at k4 with all positive gaps. The three
constant-q3 separator regions each supply a singleton factor3. This
explicit hostile verifies why gap valuations are essential rather than
optional bookkeeping. The audit retains this control in its receipt.

`audit_all_k.py` separately checks7,916 exhaustive relaxed patterns through
five contacts, including arbitrary grouping of events into triangle or
singleton blocks, both endpoint states, and free-singleton gaps. It also
checks1,562 seeded larger patterns, including contact-free triangle blocks,
high availability and odd zero-gap runs. The relaxed model admits some
unrealizable path assignments, making these hostile budget controls; it
does not replace the analytic factorization proof.

Run `python -B audit_all_k.py` and `python -O -B audit_all_k.py`.
`all_k_audit.json` contains the universe, payment counts, high-q hostile
and source SHA256. Later edits to the producer proof require a refreshed
pin or explicit retention of the earlier audited version.
