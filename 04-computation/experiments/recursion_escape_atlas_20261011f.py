"""Supplied any-distance equality -> a finite unit-distance obstruction.

The implementation uses the inherited exact Moser field.  It authenticates a
given equality seed, not its existence for six colours.  Coincidences between
copies are retained; no intermediate colourability assumption is made.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product
import json

import unit_equality_candidates_20261011e as field


proof = field.proof
CHECKS = 0


def need(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def edge(u, v):
    return (min(u, v), max(u, v))


@dataclass(frozen=True)
class EqualitySeed:
    points: tuple
    edges: frozenset
    ports: tuple
    colours: int
    colouring: tuple
    unequal_refutation: tuple


def verify_seed(seed):
    need(type(seed) is EqualitySeed, "exact equality-seed packet")
    points = field.validate_points(seed.points)
    n = len(points)
    field.graph_edges(n, seed.edges)
    need(seed.edges <= field.unit_edges(points), "all seed edges must be unit")
    need(type(seed.ports) is tuple and len(seed.ports) == 2
         and all(type(i) is int and 0 <= i < n for i in seed.ports)
         and seed.ports[0] != seed.ports[1], "two distinct named ports")
    need(type(seed.colours) is int and seed.colours >= 1, "exact positive palette")
    proof.verify_colouring(n, seed.edges, seed.colours, seed.colouring)
    # This extra edge is an ABSTRACT inequality, not a claimed unit segment.
    augmented = seed.edges | frozenset((edge(*seed.ports),))
    proof.verify_refutation(n, augmented, seed.colours, seed.unequal_refutation)
    return points


@dataclass(frozen=True)
class Amplification:
    seed: EqualitySeed
    copies: int
    rotation: tuple
    points: tuple
    edges: frozenset
    maps: tuple
    tips: tuple


def _geometry(seed, copies, rotation):
    points = verify_seed(seed)
    need(type(copies) is int and copies >= 1, "positive exact copy count")
    rotation = field.element(rotation)
    need(field.norm(rotation) == field.ONE, "actual unit rotation")
    u, v = seed.ports
    offsets = tuple(field.sub(p, points[u]) for p in points)
    h = offsets[v]
    tip = field.scale(copies, h)
    need(field.norm(field.times(field.sub(rotation, field.ONE), tip)) == field.ONE,
         "the two actual far tips must be unit separated")
    output, index, maps, edges = [], {}, [], set()
    for frame in (field.ONE, rotation):
        for j in range(copies):
            image = []
            for p in offsets:
                q = field.times(frame, field.add(p, field.scale(j, h)))
                if q not in index:
                    index[q] = len(output)
                    output.append(q)
                image.append(index[q])
            need(len(set(image)) == len(points), "a rigid copy remains injective")
            maps.append(tuple(image))
            edges.update(edge(image[a], image[b]) for a, b in seed.edges)
    tips = (maps[copies - 1][v], maps[2 * copies - 1][v])
    need(tips[0] != tips[1], "distinct tips")
    edges.add(edge(*tips))
    return rotation, tuple(output), frozenset(edges), tuple(maps), tips


def amplify(seed, copies, rotation):
    rotation, points, edges, maps, tips = _geometry(seed, copies, rotation)
    packet = Amplification(seed, copies, rotation, points, edges, maps, tips)
    verify_amplification(packet)
    return packet


def verify_amplification(packet):
    need(type(packet) is Amplification, "exact amplification packet")
    points = field.validate_points(packet.points)
    field.graph_edges(len(points), packet.edges)
    need(type(packet.maps) is tuple and all(type(m) is tuple
         and all(type(i) is int and 0 <= i < len(points) for i in m)
         for m in packet.maps), "exact source-to-output maps")
    need(type(packet.tips) is tuple and len(packet.tips) == 2
         and all(type(i) is int and 0 <= i < len(points) for i in packet.tips),
         "exact tip indices")
    expected = _geometry(packet.seed, packet.copies, packet.rotation)
    need((packet.rotation, points, packet.edges, packet.maps, packet.tips) == expected,
         "retain the complete canonical rigid-copy construction")
    need(packet.edges <= field.unit_edges(points), "authenticate every output edge")
    n, m = len(packet.seed.points), packet.copies
    need(len(points) <= 2 * m * (n - 1) + 1, "serial vertex upper bound")
    u, v = packet.seed.ports
    need(packet.maps[0][u] == packet.maps[m][u], "shared initial port")
    for start in (0, m):
        for j in range(m - 1):
            need(packet.maps[start + j][v] == packet.maps[start + j + 1][u],
                 "every consecutive equality seam is retained")
    # A colouring of the output would restrict to each authenticated seed,
    # equating all ports; the literal unit edge between the tips contradicts it.
    return True


def discover_seed(points, edges, ports, colours):
    """Explicit finite test helper; production verification never searches."""
    colouring, _ = proof.colour_or_refute(len(points), edges, colours)
    need(colouring is not None, "seed must actually be colourable")
    unequal, refutation = proof.colour_or_refute(
        len(points), edges | frozenset((edge(*ports),)), colours)
    need(unequal is None, "ports are not forced equal")
    seed = EqualitySeed(points, edges, ports, colours, colouring, refutation)
    verify_seed(seed)
    return seed


def critical_edge_seed(packet):
    """Finite extraction from a verified obstruction, not a new seed oracle."""
    verify_amplification(packet)
    n, k = len(packet.points), packet.seed.colours
    current = packet.edges
    no_colour, _ = proof.colour_or_refute(n, current, k)
    need(no_colour is None, "the output must reject the supplied palette")
    for e in sorted(current):
        trial = current - frozenset((e,))
        colouring, _ = proof.colour_or_refute(n, trial, k)
        if colouring is None:
            current = trial
    need(bool(current), "nonempty edge-critical obstruction")
    chosen = min(current)
    return discover_seed(packet.points, current - frozenset((chosen,)), chosen, k)


def main():
    zero, one, a, b = field.ZERO, field.ONE, field.A, field.B
    middle = field.add(field.scale(F(-7, 5), one), field.scale(F(9, 5), b))
    close = discover_seed((zero, middle, field.scale(F(1, 5), one)),
                          frozenset(((0, 1), (1, 2))), (0, 2), 2)
    turn = field.add(field.scale(F(-16, 9), one), field.scale(F(5, 3), b))
    diamond_points = (zero, one, a, field.add(one, a))
    diamond = discover_seed(diamond_points, field.unit_edges(diamond_points), (0, 3), 3)
    coincident = discover_seed((zero, a, one), frozenset(((0, 1), (1, 2))), (0, 2), 2)
    cases = (("distance_one_fifth", close, 3, turn),
             ("classical_diamond", diamond, 1, b),
             ("coincident_copies", coincident, 1, a))
    summaries = []
    for name, seed, m, rot in cases:
        packet = amplify(seed, m, rot)
        check(verify_amplification(packet), "authenticated amplification")
        colour, refutation = proof.colour_or_refute(len(packet.points), packet.edges, seed.colours)
        check(colour is None, "independent literal colouring search")
        nodes = proof.verify_refutation(len(packet.points), packet.edges, seed.colours, refutation)
        check(nodes > 0, "independent complete refutation")
        extracted = critical_edge_seed(packet)
        check(verify_seed(extracted) == packet.points, "critical edge retains source points")
        check(field.norm(field.sub(extracted.points[extracted.ports[0]],
                                   extracted.points[extracted.ports[1]])) == one,
              "extracted equality has actual unit-separated ports")
        summaries.append(dict(name=name, palette=seed.colours, copies=m,
                              vertices=len(packet.points), edges=len(packet.edges),
                              bound=2*m*(len(seed.points)-1)+1, proof_nodes=nodes))
        # Transport the same authenticated source proof in six orientations
        # and three translations.  Coincidence maps remain source-owned.
        for j, shift in product(range(6), (zero, one, field.add(a, b))):
            unit = field.power(a, j)
            pts = tuple(field.add(field.times(unit, p), shift) for p in seed.points)
            moved = replace(seed, points=pts)
            out = amplify(moved, m, rot)
            check(verify_amplification(out), "isometric seed transport")
            check(len(out.points) == len(packet.points) and out.edges == packet.edges,
                  "rigid transport preserves exact incidences and coincidence names")
    check(summaries[0]["vertices"] == 13, "three-copy short-distance control")
    check(summaries[1]["vertices"] == 7 and summaries[1]["edges"] == 11,
          "classical seven-point spindle")
    check(summaries[2]["vertices"] == 4 < summaries[2]["bound"], "real inter-copy coincidence")

    hostile_count = 0
    def reject(call):
        nonlocal hostile_count
        try:
            call()
        except ValueError:
            hostile_count += 1
            check(True, "hostile rejected")
        else:
            raise RuntimeError("hostile was accepted")

    packet = amplify(close, 3, turn)
    reject(lambda: amplify(close, True, turn))
    reject(lambda: amplify(close, 1, turn))
    reject(lambda: amplify(close, 3, field.scale(2, one)))
    reject(lambda: verify_seed(replace(close, colours=True)))
    reject(lambda: verify_seed(replace(close, ports=(False, 2))))
    reject(lambda: verify_seed(replace(close, colouring=(0, 0, 0))))
    reject(lambda: verify_seed(replace(close, unequal_refutation=("edge", 0, 2))))
    reject(lambda: verify_seed(replace(close, edges=close.edges | frozenset(((0, 2),)))))
    reject(lambda: verify_seed(replace(close, points=(zero, middle, (0.2, 0, 0, 0)))))
    reject(lambda: verify_amplification(replace(packet, tips=(True, packet.tips[1]))))
    wrong_maps = ((True,) + packet.maps[0][1:],) + packet.maps[1:]
    reject(lambda: verify_amplification(replace(packet, maps=wrong_maps)))
    reject(lambda: verify_amplification(replace(packet, edges=packet.edges - frozenset((edge(*packet.tips),)))))

    # Cheap interfaces from the atlas: equal scalar output counts or algebraic
    # closure are not unit-distance/colouring predicates.
    for x in range(-8, 9):
        check(2*(x+2) - (2*x+2) == 2, "translation/doubling order defect")
        check((x+2)**2 - (x*x+2) == 4*x+2, "translation/square carry")
    check(field.norm(field.scale(2, one)) != one, "dilation does not preserve unit edges")
    check(field.times(one, one) == field.times(field.scale(-1, one), field.scale(-1, one)),
          "squaring can identify distinct source points")
    left = ((7, 11), (8, 12))
    check([2*v-u for u, v in left] == [15, 16], "same difference loses raw-product data")
    for x, y, z in product(range(5), repeat=3):
        if x+y+z != 4:
            continue
        for a0, b0, c0 in product(range(5), repeat=3):
            if a0+b0+c0 != 4:
                continue
            delta = (a0-x, b0-y, c0-z)
            if sorted(delta) == [-1, 0, 1]:
                check((x-y) % 3 != (a0-b0) % 3, "edgewise triangular rank-three colouring")
    gx, gy = F(-4), F(28)
    slope = (3*gx*gx + 218*gx + 224)/(2*gy)
    dx = slope*slope - 109 - 2*gx
    dy = slope*(gx-dx) - gy
    check((dx, dy) == (F(676, 49), F(55796, 343)), "exact elliptic duplication")
    check(dy*dy == dx**3 + 109*dx*dx + 224*dx, "duplication stays on its source curve")

    print("recursion_escape_atlas_20261011f")
    print(json.dumps(summaries, sort_keys=True))
    print("Scope: three supplied k=2/3 seeds; 54 rigid transports; no supplied six-equality seed.")
    print("Coincidences retained; independent finite refutations and critical-edge extraction pass.")
    print(f"Hostiles: {hostile_count}; checks: {CHECKS}")


if __name__ == "__main__":
    main()
