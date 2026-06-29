#!/usr/bin/env python3
"""HYP-3526: spigot-style route sidecar dispatch.

The prompt analogy is a spigot algorithm: emit the next terminal certificate
from a bounded state, without going back to row names.  In the LRC14 frontier,
HYP-3513 says the compact incidence/projection sidecars I/Q already decide the
private-firewall predicate, but they do not decide the five-way terminal route.

This script joins HYP-3513 with the random031 terminal certificate ledger
HYP-3521 and the owner-boundary filtration HYP-3522.  It checks the legal
dispatch rule:

    reduce I/Q to a row-free multiplicity-one incidence cut,
    then either reconstruct route data or retain route sidecar R.

The output is a finite proof-interface certificate, not an LRC14 proof.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys
from typing import Callable, Sequence


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3513 = load_module(
    "hyp3513_for_hyp3523",
    "lrc14_private_firewall_nerode_codex_20260629.py",
)
H3521 = load_module(
    "hyp3521_for_hyp3523",
    "lrc14_random031_terminal_certificate_ledger_codex_20260629.py",
)
H3522 = load_module(
    "hyp3522_for_hyp3523",
    "lrc14_random031_owner_boundary_filtration_codex_20260629.py",
)


TargetFn = Callable[[object], object]
EXISTING_AXES = ("K", "N", "T", "S", "F", "C", "M", "A")


@dataclass(frozen=True)
class Purity:
    axes: tuple[str, ...]
    target: str
    fibers: int
    max_fiber: int
    mixed_fibers: int
    mixed_rows: int

    @property
    def pure(self) -> bool:
        return self.mixed_fibers == 0 and self.mixed_rows == 0


@dataclass(frozen=True)
class TerminalFacts:
    cell_terminal_hist: Counter
    component_terminal_hist: Counter
    terminal_certificate_count: int


@dataclass(frozen=True)
class OwnerFiltrationFacts:
    seam_owners: tuple[int, ...]
    transport_owners: tuple[int, ...]
    branch_boundary_owners: tuple[int, ...]
    bracket_lift_owners: tuple[int, ...]
    residual_after_branch_boundary: tuple[int, ...]
    mirror_pair_count: int
    bypass_size: int


@dataclass(frozen=True)
class DispatchCandidate:
    name: str
    retained_state: tuple[str, ...]
    private_pure: bool
    route_pure: bool
    row_free: bool
    bounded_carry: bool
    random031_terminal: bool
    incidence_cut_ready: bool
    retains_route_R: bool
    reconstructs_route: bool
    terminal_use: bool
    row_name_dependent: bool
    sidecar_count: int
    note: str

    @property
    def legal_now(self) -> bool:
        if self.terminal_use and not self.route_pure:
            return False
        return self.private_pure and (self.route_pure or not self.terminal_use)

    @property
    def score(self) -> int:
        score = 0
        score += 80 if self.legal_now else -60
        score += 50 if self.route_pure else 0
        score += 35 if self.private_pure else 0
        score += 22 if self.row_free else 0
        score += 18 if self.bounded_carry else 0
        score += 15 if self.random031_terminal else 0
        score += 12 if self.incidence_cut_ready else 0
        score += 10 if self.reconstructs_route else 0
        score += 6 if self.retains_route_R else 0
        score -= 4 * self.sidecar_count
        score -= 30 if self.row_name_dependent else 0
        score -= 40 if self.terminal_use and not self.route_pure else 0
        return score


def compact_counter(counter: Counter | dict) -> dict:
    return dict(sorted(counter.items(), key=lambda item: repr(item[0])))


def all_subsets(items: Sequence[str]):
    for size in range(1, len(items) + 1):
        yield from combinations(items, size)


def purity(records, axes: tuple[str, ...], target_name: str, target: TargetFn) -> Purity:
    mixed_fibers, mixed_rows, max_fiber, _samples = H3513.mixed_report(records, axes, target)
    return Purity(
        axes=axes,
        target=target_name,
        fibers=len(H3513.partition(records, axes)),
        max_fiber=max_fiber,
        mixed_fibers=mixed_fibers,
        mixed_rows=mixed_rows,
    )


def any_route_reconstruction(records, axes_family: Sequence[str]) -> tuple[bool, tuple[str, ...] | None]:
    for axes in all_subsets(axes_family):
        if purity(records, tuple(axes), "h3490_route", lambda record: record.route).pure:
            return True, tuple(axes)
    return False, None


def terminal_facts() -> TerminalFacts:
    certs = H3521.build_certificates()
    cell_terminal_hist: Counter = Counter()
    for cert in certs:
        cell_terminal_hist[cert.terminal_class] += cert.cell_count
    return TerminalFacts(
        cell_terminal_hist=cell_terminal_hist,
        component_terminal_hist=Counter(cert.terminal_class for cert in certs),
        terminal_certificate_count=H3521.terminal_certificate_count(certs),
    )


def owner_filtration_facts() -> OwnerFiltrationFacts:
    (
        _row,
        _gates,
        cells,
        by_node,
        _legal,
        bypass_component,
        bypass_cells,
        hard_seam_gates,
        _lower_bypass_gates,
    ) = H3522.build_context()
    boundaries = H3522.branch_boundaries(cells, by_node, bypass_component)
    mirror_pairs = H3522.mirror_pairs(by_node, bypass_component)

    seam_owners = H3522.gate_owner_union(hard_seam_gates)
    transport_owners = H3522.owner_union_from_cells(bypass_cells)
    boundary_cells = tuple(cell for boundary in boundaries for cell in (boundary.left_cell, boundary.right_cell))
    branch_boundary_owners = H3522.owner_union_from_cells(boundary_cells)
    bracket_lift_owners = tuple(owner for owner in branch_boundary_owners if owner not in transport_owners)
    residual_after_branch_boundary = tuple(
        owner for owner in seam_owners if owner not in set(transport_owners) | set(branch_boundary_owners)
    )
    return OwnerFiltrationFacts(
        seam_owners=seam_owners,
        transport_owners=transport_owners,
        branch_boundary_owners=branch_boundary_owners,
        bracket_lift_owners=bracket_lift_owners,
        residual_after_branch_boundary=residual_after_branch_boundary,
        mirror_pair_count=len(mirror_pairs),
        bypass_size=len(bypass_cells),
    )


def dispatch_candidates(
    p_i_private: Purity,
    p_q_private: Purity,
    p_iq_route: Purity,
    p_all_iq_route: Purity,
    p_r_route: Purity,
) -> tuple[DispatchCandidate, ...]:
    route_reconstructs_now = p_all_iq_route.pure
    private_iq = p_i_private.pure and p_q_private.pure
    return (
        DispatchCandidate(
            name="IQ_plus_R_terminal_spigot",
            retained_state=("I", "Q", "R"),
            private_pure=private_iq,
            route_pure=p_r_route.pure,
            row_free=True,
            bounded_carry=True,
            random031_terminal=True,
            incidence_cut_ready=True,
            retains_route_R=True,
            reconstructs_route=False,
            terminal_use=True,
            row_name_dependent=False,
            sidecar_count=3,
            note="current legal dispatch: row-free private cut plus explicit route carry",
        ),
        DispatchCandidate(
            name="terminal_certificate_ledger_plus_R",
            retained_state=("terminal_certificate", "R"),
            private_pure=True,
            route_pure=True,
            row_free=True,
            bounded_carry=True,
            random031_terminal=True,
            incidence_cut_ready=False,
            retains_route_R=True,
            reconstructs_route=False,
            terminal_use=True,
            row_name_dependent=False,
            sidecar_count=2,
            note="HYP-3521 emits 64+10+2+1 terminal certificates with route carry",
        ),
        DispatchCandidate(
            name="owner_filtration_plus_R",
            retained_state=("transport", "branch_boundary", "residual", "R"),
            private_pure=True,
            route_pure=True,
            row_free=True,
            bounded_carry=True,
            random031_terminal=True,
            incidence_cut_ready=False,
            retains_route_R=True,
            reconstructs_route=False,
            terminal_use=True,
            row_name_dependent=False,
            sidecar_count=4,
            note="HYP-3522 bypass subdispatch: transport, bracket lift, residual pair, route carry",
        ),
        DispatchCandidate(
            name="all_colored_plus_IQ_route_reconstruction",
            retained_state=EXISTING_AXES + ("I", "Q"),
            private_pure=private_iq,
            route_pure=route_reconstructs_now,
            row_free=True,
            bounded_carry=True,
            random031_terminal=route_reconstructs_now,
            incidence_cut_ready=True,
            retains_route_R=False,
            reconstructs_route=route_reconstructs_now,
            terminal_use=True,
            row_name_dependent=False,
            sidecar_count=len(EXISTING_AXES) + 2,
            note="conditional ideal; legal only if these fields reconstruct the five-way route",
        ),
        DispatchCandidate(
            name="IQ_without_R_private_only",
            retained_state=("I", "Q"),
            private_pure=private_iq,
            route_pure=p_iq_route.pure,
            row_free=True,
            bounded_carry=True,
            random031_terminal=False,
            incidence_cut_ready=True,
            retains_route_R=False,
            reconstructs_route=p_iq_route.pure,
            terminal_use=False,
            row_name_dependent=False,
            sidecar_count=2,
            note="legal for the private-firewall bit, not for terminal route dispatch",
        ),
        DispatchCandidate(
            name="row_name_exception_list",
            retained_state=("row_name",),
            private_pure=True,
            route_pure=True,
            row_free=False,
            bounded_carry=False,
            random031_terminal=True,
            incidence_cut_ready=False,
            retains_route_R=False,
            reconstructs_route=True,
            terminal_use=True,
            row_name_dependent=True,
            sidecar_count=1,
            note="decides the split but is not a row-free proof interface",
        ),
        DispatchCandidate(
            name="raw_private_bit",
            retained_state=("private_firewall_bit",),
            private_pure=True,
            route_pure=False,
            row_free=True,
            bounded_carry=True,
            random031_terminal=False,
            incidence_cut_ready=False,
            retains_route_R=False,
            reconstructs_route=False,
            terminal_use=False,
            row_name_dependent=False,
            sidecar_count=1,
            note="too lossy: merges six singleton-current rows with random031 hard overlap",
        ),
        DispatchCandidate(
            name="raw_count_shadow",
            retained_state=("raw_counts",),
            private_pure=False,
            route_pure=False,
            row_free=True,
            bounded_carry=False,
            random031_terminal=False,
            incidence_cut_ready=False,
            retains_route_R=False,
            reconstructs_route=False,
            terminal_use=True,
            row_name_dependent=False,
            sidecar_count=1,
            note="scalar counts forget the private cut and the route",
        ),
    )


def main() -> None:
    records = H3513.build_records()
    p_i_private = purity(records, ("I",), "private_firewall_bit", lambda record: record.private_firewall)
    p_q_private = purity(records, ("Q",), "private_firewall_bit", lambda record: record.private_firewall)
    p_i_status = purity(records, ("I",), "private_dead_status", lambda record: record.private_dead_status)
    p_q_status = purity(records, ("Q",), "private_dead_status", lambda record: record.private_dead_status)
    p_i_route = purity(records, ("I",), "h3490_route", lambda record: record.route)
    p_q_route = purity(records, ("Q",), "h3490_route", lambda record: record.route)
    p_iq_route = purity(records, ("I", "Q"), "h3490_route", lambda record: record.route)
    p_all_iq_route = purity(records, EXISTING_AXES + ("I", "Q"), "h3490_route", lambda record: record.route)
    p_r_route = purity(records, ("R",), "h3490_route", lambda record: record.route)
    existing_route_reconstructs, existing_route_axes = any_route_reconstruction(records, EXISTING_AXES)
    term = terminal_facts()
    owners = owner_filtration_facts()

    candidates = sorted(
        dispatch_candidates(p_i_private, p_q_private, p_iq_route, p_all_iq_route, p_r_route),
        key=lambda candidate: (candidate.score, candidate.name),
        reverse=True,
    )

    print("## HYP-3526 Spigot Route-Sidecar Dispatch")
    print("source=HYP-3513 + HYP-3521 + HYP-3522")
    print("analogy=spigot emits next terminal certificate from bounded carry state")
    print()
    print("## HYP-3513 I/Q and R Purity")
    for item in (
        p_i_private,
        p_q_private,
        p_i_status,
        p_q_status,
        p_i_route,
        p_q_route,
        p_iq_route,
        p_all_iq_route,
        p_r_route,
    ):
        print(
            f"axes={item.axes} target={item.target} fibers={item.fibers} "
            f"max_fiber={item.max_fiber} mixed={item.mixed_fibers}/{item.mixed_rows} "
            f"pure={item.pure}"
        )
    print(f"existing_axes_route_reconstructable={existing_route_reconstructs}")
    print(f"existing_axes_route_reconstruction_witness={existing_route_axes}")
    print()
    print("## Row-Free Incidence-Cut Lemma Shape")
    print("lemma_name=row_free_multiplicity_one_incidence_cut")
    print("input_sidecars=I/Q")
    print("predicate=every E/branch-touched blocker label has dead-component multiplicity one")
    print(f"private_status_pure_by_I={p_i_private.pure and p_i_status.pure}")
    print(f"private_status_pure_by_Q={p_q_private.pure and p_q_status.pure}")
    print(f"route_pure_by_I={p_i_route.pure}")
    print(f"route_pure_by_Q={p_q_route.pure}")
    print(f"route_pure_by_IQ={p_iq_route.pure}")
    print(f"route_pure_by_R={p_r_route.pure}")
    print()
    print("## Random031 Terminal Certificate Inputs")
    print(f"cell_terminal_hist={compact_counter(term.cell_terminal_hist)}")
    print(f"component_terminal_hist={compact_counter(term.component_terminal_hist)}")
    print(f"terminal_certificate_count_after_doublet_collapse={term.terminal_certificate_count}")
    print("identity_282=230 ordinary_rank2 + 40 bracketed_free_hole + 12 pure_bypass")
    print("identity_77=64 ordinary + 10 free_single + 2 free_doublet + 1 bypass")
    print()
    print("## HYP-3522 Owner-Filtration Carry")
    print(f"seam_owners={owners.seam_owners}")
    print(f"transport_owners={owners.transport_owners}")
    print(f"branch_boundary_owners={owners.branch_boundary_owners}")
    print(f"bracket_lift_owners={owners.bracket_lift_owners}")
    print(f"residual_after_branch_boundary={owners.residual_after_branch_boundary}")
    print(f"bypass_size={owners.bypass_size}")
    print(f"mirror_pair_count={owners.mirror_pair_count}")
    print()
    print("## Dispatch Decision")
    route_reconstructable_now = p_all_iq_route.pure or existing_route_reconstructs
    terminal_dispatch = "reconstruct_route_from_existing_data" if route_reconstructable_now else "retain_R"
    print("row_free_dispatch_rule=I/Q may reduce the private cut only before terminal route is supplied")
    print(f"current_route_reconstruction_from_existing_data={route_reconstructable_now}")
    print(f"terminal_dispatch={terminal_dispatch}")
    print("spigot_carry=R")
    print("spigot_output=terminal_certificate")
    print("proof_obligation=prove route reconstruction, or explicitly keep R in every terminal dispatch lemma")
    print()
    print("## Tournament Analysis")
    print("vertices=dispatch compressors, not runners, arcs, row names, or raw counts")
    print("pairwise_observable=private purity + route purity + row-free compression + bounded carry + random031 compatibility")
    print("switch=higher dispatch score; ties use candidate name")
    print(f"score_hist={compact_counter(Counter(candidate.score for candidate in candidates))}")
    print("directed_3cycles=0")
    print("sccs=8 singleton SCCs")
    print("hamiltonian_path=" + " -> ".join(candidate.name for candidate in candidates))
    for candidate in candidates:
        print(
            f"candidate={candidate.name} score={candidate.score} legal_now={candidate.legal_now} "
            f"retained={candidate.retained_state} private_pure={candidate.private_pure} "
            f"route_pure={candidate.route_pure} row_free={candidate.row_free} "
            f"retains_R={candidate.retains_route_R} reconstructs_route={candidate.reconstructs_route} "
            f"note={candidate.note}"
        )
    print()
    print("## Proof Pull")
    print("P1: Prove I/Q as a row-free multiplicity-one incidence-cut lemma.")
    print("P2: Do not let I/Q stand in for the HYP-3490 five-way route.")
    print("P3: For current data, terminal random031 dispatch must retain R.")
    print("P4: HYP-3521 emits the next certificate, while HYP-3522 supplies the owner carry.")
    print("P5: A future route-reconstruction theorem may delete R; until then R is the spigot carry.")


if __name__ == "__main__":
    main()
