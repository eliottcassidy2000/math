#!/usr/bin/env python3
"""
Projective-matroid annealer for the equal-suffix arithmetic-Kakeya game.

The search uses the normal form A=[B|D_rho B] and distinguished-coloop
peeling from ak_projective_matroid_audit.py.  Inner-loop ranks are computed
over two large prime fields; every record is replayed over Q by the
independent exact audit before it is printed.

This is a search instrument, not a proof.  Reproduction examples:

  python3 ak_matroid_anneal.py 2,2,2 20000 1
  python3 ak_matroid_anneal.py 2,2,2,2 50000 2
"""

from fractions import Fraction
from itertools import product
import importlib.util
import math
import random
import sys
import time
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


AUDIT_PATH = Path(__file__).with_name("ak_projective_matroid_audit.py")
SPEC = importlib.util.spec_from_file_location("ak_projective_audit", AUDIT_PATH)
AUDIT = importlib.util.module_from_spec(SPEC)
require(SPEC.loader is not None, "cannot load audit module")
SPEC.loader.exec_module(AUDIT)


RHO = (
    Fraction(-3), Fraction(-2), Fraction(-1), Fraction(-1, 2),
    Fraction(-1, 3), Fraction(0), Fraction(1, 3), Fraction(1, 2),
    Fraction(1), Fraction(2), Fraction(3),
)
PRIMES = (1000003, 1000033)


def rho_mod(rho, prime):
    return (rho.numerator * pow(rho.denominator, -1, prime)) % prime


def rho_label(rho):
    """Integral (a,b) with (a-b)/(a+b)=rho and a+b != 0."""
    p, q = rho.numerator, rho.denominator
    return (q + p, q - p)


def rref_mod(rows, ncols, prime):
    work = [list(row) for row in rows if any(x % prime for x in row)]
    rr = 0
    pivot_positions = []
    for col in range(ncols):
        pivot = next(
            (i for i in range(rr, len(work)) if work[i][col] % prime), None
        )
        if pivot is None:
            continue
        work[rr], work[pivot] = work[pivot], work[rr]
        inv = pow(work[rr][col] % prime, -1, prime)
        work[rr] = [(x * inv) % prime for x in work[rr]]
        for i in range(len(work)):
            if i == rr:
                continue
            scale = work[i][col] % prime
            if scale:
                work[i] = [
                    (x - scale * y) % prime
                    for x, y in zip(work[i], work[rr])
                ]
        pivot_positions.append((col, rr))
        rr += 1
        if rr == len(work):
            break
    return [(col, work[i]) for col, i in pivot_positions]


class Shape:
    def __init__(self, dims):
        self.dims = tuple(dims)
        self.vertices = tuple(product(*[range(1, d + 1) for d in dims]))
        self.vid = {v: i for i, v in enumerate(self.vertices)}
        self.n = len(self.vertices)
        self.groups = []
        for axis, d in enumerate(dims):
            domain = product(
                *([range(1, x + 1) for x in dims[:axis]]
                  + [range(1, d)])
            )
            suffixes = tuple(product(
                *[range(1, x + 1) for x in dims[axis + 1:]]
            ))
            for pre in domain:
                right = pre[:axis] + (pre[axis] + 1,)
                edges = tuple(
                    (self.vid[pre + suffix], self.vid[right + suffix])
                    for suffix in suffixes
                )
                self.groups.append((axis, pre, len(edges), edges))

    def initial_state(self):
        groups = tuple(-1 for _ in self.groups)
        seeds = frozenset(
            (v, r)
            for v in range(self.n)
            for r in (RHO.index(Fraction(-1)), RHO.index(Fraction(1)))
        )
        return groups, seeds, frozenset()

    def known_lift_state(self):
        """Lift the verified strict 13/7 [2,2,2] gadget over any suffix."""
        if len(self.dims) < 3 or self.dims[:3] != (2, 2, 2):
            return self.initial_state()
        choice = {
            (0, (1,)): RHO.index(Fraction(-1, 3)),
            (1, (1, 1)): RHO.index(Fraction(1)),
            (1, (2, 1)): RHO.index(Fraction(1)),
            (2, (2, 1, 1)): RHO.index(Fraction(1, 3)),
            (2, (2, 2, 1)): RHO.index(Fraction(-1)),
        }
        groups = tuple(
            choice.get((axis, pre), -1)
            for axis, pre, _, _ in self.groups
        )
        suffixes = tuple(product(
            *[range(1, d + 1) for d in self.dims[3:]]
        ))
        seeds = set()
        t0 = set()
        seed_template = (
            ((1, 2, 2), Fraction(-1)),
            ((2, 2, 1), Fraction(0)),
            ((1, 1, 1), Fraction(1)),
        )
        for suffix in suffixes:
            for pre, rho in seed_template:
                seeds.add((self.vid[pre + suffix], RHO.index(rho)))
            t0.add(self.vid[(1, 1, 2) + suffix])
        return groups, frozenset(seeds), frozenset(t0)

    def cost(self, state):
        groups, seeds, _ = state
        return (
            sum(self.groups[i][2] for i, r in enumerate(groups) if r >= 0)
            + len(seeds)
        )

    def rows(self, state, live, prime):
        groups, seeds, _ = state
        order = sorted(live)
        pos = {v: j for j, v in enumerate(order)}
        rmods = tuple(rho_mod(r, prime) for r in RHO)
        rows = []
        for i, ridx in enumerate(groups):
            if ridx < 0:
                continue
            rho = rmods[ridx]
            for u, v in self.groups[i][3]:
                row = [0] * (2 * len(order))
                if u in live:
                    j = pos[u]
                    row[2 * j], row[2 * j + 1] = 1, rho
                if v in live:
                    j = pos[v]
                    row[2 * j] = (row[2 * j] - 1) % prime
                    row[2 * j + 1] = (row[2 * j + 1] - rho) % prime
                if any(row):
                    rows.append(row)
        for v, ridx in seeds:
            if v not in live:
                continue
            row = [0] * (2 * len(order))
            j = pos[v]
            row[2 * j], row[2 * j + 1] = 1, rmods[ridx]
            rows.append(row)
        return order, rows

    def force_mod(self, state, prime):
        _, _, t0 = state
        T = set(t0)
        while len(T) < self.n:
            live = set(range(self.n)) - T
            order, rows = self.rows(state, live, prime)
            if not rows:
                break
            reduced = rref_mod(rows, 2 * len(order), prime)
            fired = []
            for pivot, row in reduced:
                if pivot % 2 != 1:
                    continue
                if sum(x % prime != 0 for x in row) == 1:
                    fired.append(order[pivot // 2])
            if not fired:
                break
            T.update(fired)
        return len(T), T

    def exact_instance(self, state):
        group_choices, seeds, t0 = state
        fs = [dict() for _ in self.dims]
        for i, ridx in enumerate(group_choices):
            if ridx < 0:
                continue
            axis, pre, _, _ = self.groups[i]
            fs[axis][pre] = rho_label(RHO[ridx])
        seed_data = [
            (self.vertices[v], rho_label(RHO[ridx]))
            for v, ridx in sorted(seeds)
        ]
        t_data = [self.vertices[v] for v in sorted(t0)]
        return AUDIT.Instance(self.dims, fs, seed_data, t_data)

    def exact_force(self, state):
        instance = self.exact_instance(state)
        return AUDIT.projective_force(instance, "strict")[0], instance


def mutate(shape, state, rng):
    groups, seeds, t0 = state
    groups = list(groups)
    seeds = set(seeds)
    t0 = set(t0)
    move = rng.random()
    if move < 0.48:
        i = rng.randrange(len(groups))
        old = groups[i]
        choices = tuple(range(-1, len(RHO)))
        groups[i] = rng.choice(tuple(x for x in choices if x != old))
    elif move < 0.92:
        if seeds and rng.random() < 0.46:
            seeds.remove(rng.choice(tuple(seeds)))
        else:
            item = (rng.randrange(shape.n), rng.randrange(len(RHO)))
            if item in seeds:
                seeds.remove(item)
            else:
                seeds.add(item)
    else:
        v = rng.randrange(shape.n)
        if v in t0:
            t0.remove(v)
        elif len(t0) < 2:
            t0.add(v)
    return tuple(groups), frozenset(seeds), frozenset(t0)


def anneal(dims, steps, seed):
    shape = Shape(dims)
    rng = random.Random(seed)
    state = shape.known_lift_state()
    cache = {}

    def evaluate(candidate):
        if candidate in cache:
            return cache[candidate]
        forced1, _ = shape.force_mod(candidate, PRIMES[0])
        cost = shape.cost(candidate)
        den = shape.n - len(candidate[2])
        # A missing vertex costs more than the full target's overhead, but
        # leaves finite-temperature routes through temporarily broken states.
        energy = Fraction(cost * 10 + 22 * (shape.n - forced1), 10 * den)
        cache[candidate] = (energy, forced1)
        return energy, forced1

    energy, forced = evaluate(state)
    initial_exact, initial_instance = shape.exact_force(state)
    require(initial_exact, "initial state must force exactly")
    best_score = initial_instance.score
    best_state = state
    exact_records = 0
    start = time.time()
    for step in range(steps):
        candidate = mutate(shape, state, rng)
        if rng.random() < 0.18:
            candidate = mutate(shape, candidate, rng)
        new_energy, new_forced = evaluate(candidate)
        temp = 0.25 * (1 - step / max(1, steps)) + 0.015
        delta = float(new_energy - energy)
        if delta <= 0 or rng.random() < math.exp(-delta / temp):
            state, energy, forced = candidate, new_energy, new_forced

        if new_forced != shape.n:
            continue
        forced2, _ = shape.force_mod(candidate, PRIMES[1])
        if forced2 != shape.n:
            continue
        cost = shape.cost(candidate)
        den = shape.n - len(candidate[2])
        score = Fraction(cost, den)
        if score >= best_score:
            continue
        exact_ok, instance = shape.exact_force(candidate)
        if not exact_ok:
            continue
        exact_records += 1
        best_score, best_state = score, candidate
        print(
            f"RECORD dims={dims} seed={seed} step={step} score={score} "
            f"m={instance.m} r={instance.r} n={instance.n} t={instance.t} "
            f"elapsed={time.time()-start:.2f}s",
            flush=True,
        )
        print(f"  fs={instance.fs}", flush=True)
        print(f"  seeds={instance.seeds} t0={instance.t0}", flush=True)
        if score <= Fraction(1675, 1000):
            print("TARGET_LE_1.675", flush=True)
    print(
        f"DONE dims={dims} seed={seed} steps={steps} best={best_score} "
        f"exact_records={exact_records} states={len(cache)} "
        f"elapsed={time.time()-start:.2f}s",
        flush=True,
    )
    return best_score, best_state


def main():
    dims = tuple(map(int, sys.argv[1].split(","))) if len(sys.argv) > 1 else (2, 2, 2)
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    anneal(dims, steps, seed)


if __name__ == "__main__":
    main()
