#!/usr/bin/env python3
"""Reproduce a six-crossing false TRUE_CERTIFIED verdict.

Spherogram 2.4.1 is installed under ./vendor.  Starting from a one-crossing
unknot, its reverse Reidemeister operations construct the hostile PD by
R1+, R2+, R2+, R3, R3.  Thus the hostile input is still the unknot.  The
current repo engine's greedy input check stalls, while changing crossing 4
makes its own reducer reach the empty diagram.  It therefore labels an
u=0 knot TRUE_CERTIFIED for u=1.
"""

import importlib.util
import json
import random
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "vendor"))

from spherogram import Link  # noqa: E402
from spherogram.links import simplify  # noqa: E402


EXPECTED_FINAL = [
    [1, 11, 2, 10],
    [6, 10, 7, 9],
    [3, 8, 4, 9],
    [11, 5, 12, 4],
    [7, 2, 8, 3],
    [5, 1, 6, 12],
]
EXPECTED_TYPES = [1, 2, 2, 3, 3]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def one_indexed_pd(link):
    return [[label + 1 for label in row] for row in link.PD_code()]


def construct():
    # This fixed run records a short path found by a hostile random search.
    # Calls to the repo decider consume no random state, so the construction
    # can be replayed without doing any earlier invariant computations.
    random.seed(31415926)
    chosen_path = None
    chosen_types = None
    for trial in range(3835):
        link = Link([(0, 0, 1, 1)])
        path = [one_indexed_pd(link)]
        types = []
        next_label = 0
        for _ in range(random.randint(2, 9)):
            move_type = random.choices([1, 2, 3], [.18, .32, .5])[0]
            next_label += move_type % 3
            before = one_indexed_pd(link)
            simplify.random_reverse_move(link, move_type, next_label)
            link._rebuild(same_components_and_orientations=True)
            after = one_indexed_pd(link)
            if after != before:
                types.append(move_type)
                path.append(after)
        if trial == 3834:
            chosen_path, chosen_types = path, types
            break
    require(chosen_types == EXPECTED_TYPES, "reverse-move types changed")
    require(chosen_path[-1] == EXPECTED_FINAL, "hostile PD changed")
    return chosen_path, chosen_types


def load_repo_engine():
    engine_path = HERE.parents[1] / "04-computation" / "unknot1_decider.py"
    spec = importlib.util.spec_from_file_location("unknot1_decider", engine_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    path, move_types = construct()
    engine = load_repo_engine()
    diagram = engine.parse_pd(EXPECTED_FINAL)
    reduced, initial_log, initial_success = engine.try_reduce(diagram)
    result = engine.decide(EXPECTED_FINAL)

    require(not initial_success, "input reducer unexpectedly recognized unknot")
    require(result["verdict"] == engine.VERDICT_TRUE,
            "engine no longer emits the hostile TRUE verdict")
    require(result["invariants"]["det"] == 1, "hostile determinant changed")
    require(result["invariants"]["sigma"] == 0, "hostile signature changed")

    print("construction move types =", move_types,
          "(R1+, R2+, R2+, R3, R3)")
    for number, pd in enumerate(path):
        print("path[%d] =" % number, json.dumps(pd))
    print("input greedy log =", initial_log,
          "; stalled crossings =", reduced.n)
    print("engine verdict =", result["verdict"])
    print("engine certificate =", result["certificate"])
    print("mathematical truth = input is the unknot, hence u(K)=0")


if __name__ == "__main__":
    main()
