#!/usr/bin/env python3
"""HYP-2673 scout: split the LRC14 open constant into proof currencies.

The incoming HYP-2653c work first corrected the old "one open constant" story:
the raw far-element quantity

    C(k) = sup_{E',w} w*|Delta_w|

is not a tiny non-resonant scalar.  The later HYP-2653d correction is sharper:
`w*Delta_w` is not the proof currency at all, because Delta_w has a small
nonzero floor along resonant dyadic families and so w*Delta_w grows with scale.
The proof target is instead the uniform far-tail cap

    sup_{max(E') > B} Delta_w(E', w) <= cap_k - Q(k-1).

The post-shell HYP-2671 work, meanwhile, localizes a different constant: the
shell-full boundary-tax ratio Delta_w^+/p1(E').

This script records the exact constants side by side.  Its goal is not another
large search; it is an exact proof-route ledger showing which normalization
each constant belongs to and which remaining lemma it asks for.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Constant:
    name: str
    value: Fraction
    currency: str
    role: str


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def ceil_frac(q: Fraction) -> int:
    return -(-q.numerator // q.denominator)


def print_constant_table(constants: list[Constant]) -> None:
    print("constant stack")
    print("name | value | currency | proof role")
    for c in constants:
        print(f"{c.name} | {fmt(c.value)} | {c.currency} | {c.role}")
    print()


def print_gap_table() -> None:
    thr2 = Fraction(426, 35035)
    drop6 = Fraction(7, 858)
    tower = {
        "delete 1": Fraction(1333, 30940),
        "delete 2": Fraction(27, 1547),
        "delete 4": Fraction(335, 23023),
        "delete 8": Fraction(6163, 336336),
    }
    finite_leader = Fraction(997, 2562)
    new_leader = Fraction(1371, 4319)
    tail_leader = Fraction(932669, 4085893)

    print("exact slack ledger")
    print(f"drop-6 mouth floor                 = {fmt(drop6)}")
    print(f"AP one-hole second threshold        = {fmt(thr2)}")
    print(f"second threshold - drop-6           = {fmt(thr2 - drop6)}")
    print()
    print("tower-deletion margins above 426/35035")
    for name, value in tower.items():
        print(f"  {name:8s}: value={fmt(value)}, margin={fmt(value - thr2)}")
    print()
    print("shell-full p1-tax margins")
    print(f"  finite B13 leader below 2/5       = {fmt(Fraction(2, 5) - finite_leader)}")
    print(f"  new-speed leader below 1/3        = {fmt(Fraction(1, 3) - new_leader)}")
    print(f"  far-tail leader below 1/4         = {fmt(Fraction(1, 4) - tail_leader)}")
    print()


def print_uniform_tail_table() -> None:
    # The tight k=9 far-peel row in the sector route compares cap_9 to Q(8).
    cap9 = Fraction(1979, 4004)
    q8 = Fraction(621, 1715)
    margin = cap9 - q8

    print("corrected far-tail target")
    print(f"cap_9 - Q(8) margin                 = {fmt(margin)}")
    print("target")
    print("  prove sup_{max(E')>B} Delta_w(E',w) <= cap_k-Q(k-1)")
    print("  with a finite check for max(E')<=B, empirically B around 20.")
    print("HYP-2653d scout")
    print("  B=14 is already below margin but tight at k=9.")
    print("  B=20 gives about 2.3x safety in the KPS adversarial scan.")
    print("  the B=14 worst row is the HYP-2671 dyadic-block extremizer.")
    print()


def print_superseded_span_diagnostic() -> None:
    # These are no longer a proof target.  They explain why the old constant
    # framing looked plausible before HYP-2653d corrected the normalization.
    cap9 = Fraction(1979, 4004)
    q8 = Fraction(621, 1715)
    margin = cap9 - q8

    candidates = [
        ("old non-resonant target", Fraction(39, 20)),
        ("KPS scale-count example", Fraction(293, 100)),
        ("HYP-2655 multiscale sample", Fraction(2804017, 717360)),
        ("linear k=9 toy bound", Fraction(9, 1)),
    ]

    print("superseded w*Delta span diagnostic")
    print(f"cap_9 - Q(8) margin                 = {fmt(margin)}")
    print("candidate C | needed w >= ceil(C/margin)")
    for name, c in candidates:
        print(f"  {name:28s}: C={fmt(c)}, span >= {ceil_frac(c / margin)}")
    print()
    print("Interpretation")
    print("  These span numbers are historical diagnostics, not proof obligations.")
    print("  HYP-2653d shows Delta_w has a nonzero resonant floor, so w*Delta_w")
    print("  grows with scale.  The far-tail proof must bound Delta_w itself")
    print("  uniformly after a finite max(E') cutoff.")
    print()


def main() -> None:
    constants = [
        Constant(
            "shell-damage threshold",
            Fraction(426, 35035),
            "fixed-observer safe mass",
            "HYP-2661/HYP-2666 gate: damaging the dyadic-1 tower pays this floor",
        ),
        Constant(
            "finite packet tax",
            Fraction(2, 5),
            "Delta_w^+/p1(E')",
            "shell-full finite pocket; B13 leader stays below this",
        ),
        Constant(
            "new-speed packet tax",
            Fraction(1, 3),
            "Delta_w^+/p1(E')",
            "HYP-2671 local constant; dyadic m=4 leader is the sharp scout row",
        ),
        Constant(
            "far-tail packet tax",
            Fraction(1, 4),
            "Delta_w^+/p1(E')",
            "suggested tail lemma for max(E')>24 after the shell-full quotient",
        ),
        Constant(
            "k=9 uniform Delta margin",
            Fraction(129643, 980980),
            "Delta_w",
            "HYP-2653d far-tail target: sup Delta_w below cap_9-Q(8)",
        ),
    ]

    print("HYP-2673 LRC14 constant stack")
    print("exact Fraction arithmetic")
    print()
    print_constant_table(constants)
    print_gap_table()
    print_uniform_tail_table()
    print_superseded_span_diagnostic()
    print("Tournament Analysis")
    print("  vertices: shell_damage_gate > finite_packet_tax > new_speed_packet_tax > far_tail_packet_tax > uniform_delta_tail > raw_runner_vertices")
    print("  observable: exact proof currency attached to each gate, not raw runner labels")
    print("  switch/gauge: normalize endpoint discrepancy either by p1(E') or by the absolute Delta_w tail cap")
    print("  Hamiltonian path: shell_damage_gate > finite_packet_tax > new_speed_packet_tax > far_tail_packet_tax > uniform_delta_tail > raw_runner_vertices")
    print("  challenged assumption: there is one scalar constant to prove.  The data")
    print("    supports a stack of currencies: local boundary mass p1, plus a")
    print("    uniform far-tail Delta cap after finite max(E') cutoff.")
    print("PASS: constant-stack scout complete.")


if __name__ == "__main__":
    main()
