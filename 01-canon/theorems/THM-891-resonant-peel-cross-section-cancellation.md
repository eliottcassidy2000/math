---
id: THM-891
title: Exact cross-section cancellation for a resonant far peel
status: CLAIMED / WRITE-UP IN PROGRESS. The limiting identity and its seven-residue reduction are proved; the machine-exact referee and canonical constant table are being added in this session.
source: codex-2026-07-16-S17
depends_on: [THM-727, THM-883, THM-884, THM-887-uniform-maxS-and-affine-witness-coordinate]
related: [THM-888, THM-889, HYP-7021, HYP-7024]
---

# THM-891 — exact cross-section cancellation for a resonant far peel

Namespace reserved for the following proved mechanism.  Let a fixed six-offset core
`E` contain `0`, append a far speed `t`, and evaluate the two-scale error at the owner
resonance `w=at`.  If `B_s` is the measure on which `E` misses exactly `{s}` and
`A_{s,c}` the measure on which it misses exactly `{s,c}`, then the limit as `t -> infinity`
is an explicit signed sum of the seven source-section integrals

`J_a(c,s) = integral_[c/7,(c+1)/7] (1_[s/7,(s+1)/7]({ay}) - 1/7) dy`.

The exact microcell count gives

`a J_a(c,s) = N_a(c,s)/7 - a/49`,

so `a` times the limiting error depends only on `a mod 7`.  This retains the
cross-section signs that are discarded in the sectionwise absolute-value bound of the
other `THM-887`, and it explains the factor-of-several gap between that bound and the
exact audits in `THM-884`.

Still being added: the full proof, exact residue table for `E={0,1,2,3,4,5}`, finite-`t`
referee, Tournament Analysis diagnostic, and the precise boundary between this theorem
and the open universal extremal inequality recorded as `HYP-7024`.
