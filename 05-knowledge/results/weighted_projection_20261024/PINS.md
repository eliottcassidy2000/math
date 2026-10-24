# The weighted bound retains the moving shared owners

Status: **PROVED conditional weighted extensions** of `PROOF.md`, instantiated
by the literal native equations and receipts in
`../three_target570_20261024/{SHARED-PIN,TWELVE-PIN}.md`.
These are two DIFFERENT pin families; their 176/368 equations must not be
conjoined. Statistical weights are never substituted for physical currents.

## One moving owner

Retain exactly `path_176[2]=rho(path_368[1])`, with the incoming literal rho.
Its bit equation is z=b_176=1-b_368. Define Z^z[r,g] by summing actual R
weights over subsets with the specified bit176. For BOTH z, all coefficients
are positive at every r=7,...,13. For outside state J the forced z is known.
Use six auxiliary reference bits, and

    F(J)=Z^z[r,q-ell(J)] / binom(6,r-7).

Let kappa be the maximum divided by the minimum over z, g and 7<=r<=13.
This is an overbound if some combinations do not occur, which is harmless.
The same pair heat bath, Metropolis tilt and conditional refresh give

    gap >= 1/(2 kappa binom(m+6,2))                        (1)

on the complete interior band. Every cross-fibre event costs <=22 sources /
54 clocks. The forced change of176 is included among the twenty R records.
Shared-owner equality is checked at both complete endpoints. Inverse and
converse use the same reflected literal paths and integer clocks as the
inherited certificate; no new capacity is assigned to those lifts.

For the full fibre, take A_r=max_{z,g} Z^z[r,g], and B_7,B_13 the minima over
both phases and all charges at those grades. Replace the two sums in H by

    sum_{d=1}^7 (A_{7-d}/B_7) e_d(w),
    sum_{d=1}^7 (A_{13+d}/B_13) e_d(w^{-1}).

Nearest-band admission changes at most14 sources and costs at most63 clocks.
The bit z is recomputed from J at every arrival. Parent mass comparison uses
maxima/minima over both phases, so it remains valid even when admission
changes368 and176 together. Thus the full-fibre gap is (1)/(8H), for
7<=k<=m+13. At the remaining extreme total grades, complete refresh changes
at most twelve sources and costs at most84 clocks. Overall budget: **22/84**.
Under uniformly bounded outside odds this is an inverse polynomial of degree
at most nine. It does not claim the original 20/63 routing budget is optimal.

## Twelve literal owner equations

Use the exact incoming twelve-pin family frozen at
`2bc5d106db075752dd28090cb5637555ede01d66`, receipt
`guarded_component_completion_20261024/native/receipt.json`.
Ten controlled sources follow or oppose root176, two follow or oppose1532,
and2948 stays U. The controls always contribute six V records. Their total
flip costs are61 for the first ten and11 for the last two. No control is
erased from the weight, charge, owner or clock ledger.

The visible state is (J,z), where J is a subset of the m freely varying
outside sources (48 or47 in the two native binary branches), and z=b_176.
Let t=b_1532, determined by J. C(z,t) is the FULL prescribed twelve-control
assignment. Put

    Y[r,g,z,t] = weight(C(z,t)) * Z^z[r,g],
    g=q-ell(J)-ell(C(z,t)).

The constant weight of fixed2948 can be canceled. The ordinary V count k
includes J and R; the physical total V count is k+6. Its integer value is
not replaced by a residue. On r=7,...,13 all phases have complete charge
coverage. Use the product reference on a size k-7 subset of m+6 coordinates,
AND a fair independent z bit. The reference chain spends half its updates
on the pair heat bath and half on refreshing z. Its gap is at least
`1/(2 binom(m+6,2))`. Tilt by `Y/binom(6,r-7)` and lazify. With kappa the
max/min of that finite potential, the interior gap is at least

    gamma = 1/(4 kappa binom(m+6,2)).                      (2)

Every reference transition has a native lift:

* an outside swap changes at most20 R records, two outside records and the
  two1532 controls: <=24 sources /65 clocks;
* a z flip changes at most20 R records and the ten176 controls:
  <=30 sources /101 clocks;
* a one-bit outside change costs <=22 sources /56 clocks;
* a visible hold refresh keeps both roots fixed and costs <=20/40.

Hence **30/101** pays the entire interior weighted chain. This is a larger
budget than the earlier route through grade7; it allows direct weighted
phase refresh at every grade, including the central grade10. The clocks
have not been divided by an occupancy or identified modulo anything.

For boundary attachment keep z fixed, admit to7 or13, and recompute the
two1532 controls from J. Work is at most74 clocks and support16. Use
`A_r=max_{g,z,t}Y[r,g,z,t]` and `B_a=min_{g,z,t}Y[a,g,z,t]` in the same
seven-term H formula. The full-fibre bound is gamma/(8H). At ordinary
grades k<=6 or m+20-k<=6, complete full-state refresh changes at most twelve
ordinary records (cost84) and all twelve controls (cost72): <=24/156.
Thus **30/156** is an all-grade weighted budget. Under bounded free-source
odds, fixed reservoir and control weights, the bound has degree at most nine.

## Why the shared-frame information matters

This construction uses all actual control weights and their charge shifts
before comparison. A charge-only trace cannot replace Y. Incoming separator
sections in `guarded_component_completion_20261024/theory/SEPARATOR_SECTIONS.md`
provide the general owned-once bookkeeping rule. Here the complete conditional
band gives the additional positivity and bounded-event hypotheses that a
separator count alone does not supply. Arbitrary larger pin components need
their own conditioned coefficients and native work bounds.
