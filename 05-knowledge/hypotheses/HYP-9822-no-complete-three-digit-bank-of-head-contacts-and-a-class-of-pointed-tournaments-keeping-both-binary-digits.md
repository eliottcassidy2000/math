---
id: HYP-9822
title: "After THM-5272: (a) no tournament base has a complete 3-bank (three head contacts whose 8 on/off patterns give 8 consecutive odd counts); true for every base of order at most 9, although strictly submodular quadratic bijections of the 3-cube exist; (b) there is a class of pointed tournaments, closed under tail extension, in which ports of defect -1 and +1 both persist from some count on, which would give every large odd number as a count by binary expansion; (c) the longest run of consecutive counts on a single base of order n is unbounded (24 at order 8, 29 at order 9)"
status: >
  OPEN (CONJECTURE / QUESTION), mac-mini-2026-11-15-carry. (a) FINITE-EXACT through base order 9 (191,536
  classes); by THM-5272 (B) a complete 3-bank needs H(B - e - f) <= 14 for its three pairs of contacts.
  (b) is a QUESTION, not supported beyond the port spectra: both digits exist at 26 of the 50 counts up to 104 and
  fail at 43, 85, 87 simultaneously; no closed class is known. (c) NUMERICAL trend only.
  Chan-Pak remains OPEN.
session: mac-mini-2026-11-15-carry
date: 2026-11-15
source: 05-knowledge/results/chanpak_carry_banks_ports_20261115/README.md
scripts:
  - 05-knowledge/results/chanpak_carry_banks_ports_20261115/scripts/bank.c, ext_runs.c, port_spectra.py
related:
  - THM-5272, THM-5262, HYP-9812, physical_radix_carry_20261116/BUDGET.md
---

# HYP-9822 - banks of head contacts and persistent binary ports

Source: [the note](../results/chanpak_carry_banks_ports_20261115/README.md), sections 2-4.

**(a)** For every tournament `B`, context `X0` and three vertices `e, f, g` outside `X0`, the eight counts
`H(B + x_{X0 ∪ Y})`, `Y ⊆ {e, f, g}`, are not eight consecutive odd numbers. Falsifier: one base with a complete
3-bank (search space by THM-5272 (B): bases with `H(B - e - f) <= 14` on the three pairs).

**(b)** Is there a family of pointed tournaments `(A, u)`, closed under adding a tail vertex at a port, such that
every member beyond some count has a port of defect `-1` and a port of defect `+1`?

**(c)** `max_B (longest run of consecutive odd numbers among the counts H(B + x))` tends to infinity with the
order of `B`.
