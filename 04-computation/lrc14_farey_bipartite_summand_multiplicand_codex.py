"""Farey packets as summand, multiplicand, and bipartite graph carriers.

This scout is aimed at the LRC(14) proof search.  A reduced Farey term a/b
with 0 < a <= b is read three ways:

    a + b          summand/pinch denominator
    a * b          multiplicand node and |E(K_{a,b})|
    K_{a,b}        complete bipartite graph carrier

The new point checked here is that F_4 is the first Farey level where this
carrier contains the Kuratowski K_{3,3} obstruction: the reduced packet 3/4
is K_{3,4}, which contains K_{3,3} by deleting one vertex on the 4-side.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Callable


N = 14
PAIR_SUM_MODULUS = 2 * N - 1
APEX = N // 2


@dataclass(frozen=True, order=True)
class Packet:
    """Reduced Farey packet a/b with its three basic carriers."""

    a: int
    b: int

    @property
    def value(self) -> Fraction:
        return Fraction(self.a, self.b)

    @property
    def summand(self) -> int:
        return self.a + self.b

    @property
    def product(self) -> int:
        return self.a * self.b

    @property
    def planar_bipartite(self) -> bool:
        return min(self.a, self.b) <= 2

    @property
    def has_k33_minor(self) -> bool:
        return min(self.a, self.b) >= 3

    @property
    def shell27(self) -> int:
        r = self.product % PAIR_SUM_MODULUS
        return min(r, (-r) % PAIR_SUM_MODULUS)

    @property
    def sector7(self) -> int:
        return self.product % 7

    def label(self) -> str:
        return f"{self.a}/{self.b}"

    def carrier(self) -> str:
        graph = f"K_{{{self.a},{self.b}}}"
        state = "nonplanar" if self.has_k33_minor else "planar"
        return (
            f"{self.label():>4}  sum={self.summand:>2}  "
            f"prod={self.product:>3}  shell27={self.shell27:>2}  "
            f"{graph:<8} {state}"
        )


def farey_packets(level: int) -> list[Packet]:
    packets = [
        Packet(a, b)
        for b in range(1, level + 1)
        for a in range(1, b + 1)
        if gcd(a, b) == 1
    ]
    return sorted(packets, key=lambda p: p.value)


def divisors(n: int) -> set[int]:
    return {d for d in range(1, n + 1) if n % d == 0}


def relation_graph(
    packets: list[Packet], relation: Callable[[Packet, Packet], bool]
) -> dict[Packet, set[Packet]]:
    graph = {p: set() for p in packets}
    for i, p in enumerate(packets):
        for q in packets[i + 1 :]:
            if relation(p, q):
                graph[p].add(q)
                graph[q].add(p)
    return graph


def component_sizes(graph: dict[Packet, set[Packet]]) -> list[int]:
    seen: set[Packet] = set()
    sizes: list[int] = []
    for start in graph:
        if start in seen:
            continue
        q = deque([start])
        seen.add(start)
        size = 0
        while q:
            p = q.popleft()
            size += 1
            for nxt in graph[p]:
                if nxt not in seen:
                    seen.add(nxt)
                    q.append(nxt)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def edge_count(graph: dict[Packet, set[Packet]]) -> int:
    return sum(len(nbrs) for nbrs in graph.values()) // 2


def graph_stats(
    packets: list[Packet], name: str, relation: Callable[[Packet, Packet], bool]
) -> tuple[str, int, int, list[int]]:
    graph = relation_graph(packets, relation)
    sizes = component_sizes(graph)
    return name, edge_count(graph), len(sizes), sizes[:8]


def farey_det1(p: Packet, q: Packet) -> bool:
    return abs(p.a * q.b - q.a * p.b) == 1


def comparable_bipartite_subgraph(p: Packet, q: Packet) -> bool:
    return (p.a <= q.a and p.b <= q.b) or (q.a <= p.a and q.b <= p.b)


def consecutive_farey_relation(packets: list[Packet]) -> Callable[[Packet, Packet], bool]:
    neighbors = set()
    for p, q in zip(packets, packets[1:]):
        neighbors.add(frozenset((p, q)))

    def relation(p: Packet, q: Packet) -> bool:
        return frozenset((p, q)) in neighbors

    return relation


def first_nonplanar_by_level(limit: int) -> Packet | None:
    for level in range(1, limit + 1):
        candidates = [p for p in farey_packets(level) if p.has_k33_minor]
        if candidates:
            return candidates[0]
    return None


def packets_with_product(level: int, product: int) -> list[Packet]:
    return [p for p in farey_packets(level) if p.product == product]


def level_table(limit: int) -> None:
    print("Farey-bipartite packet table")
    print("level terms new distinct-product-sum product-set-is-divisors nonplanar new-nonplanar")
    prev: set[Packet] = set()
    for level in range(1, limit + 1):
        packets = farey_packets(level)
        packet_set = set(packets)
        new_packets = sorted(packet_set - prev, key=lambda p: p.value)
        products = {p.product for p in packets}
        nonplanar = [p for p in packets if p.has_k33_minor]
        new_nonplanar = [p.label() for p in new_packets if p.has_k33_minor]
        divisor_target = max(products) if products else 1
        is_divisors = products == divisors(divisor_target)
        print(
            f"F_{level:<2} {len(packets):>5} {len(new_packets):>3} "
            f"{sum(products):>15} {str(is_divisors):>22} "
            f"{len(nonplanar):>10} {','.join(new_nonplanar) or '-'}"
        )
        prev = packet_set


def relation_tournament(packets: list[Packet]) -> None:
    relations = [
        ("Farey consecutive", consecutive_farey_relation(packets)),
        ("Farey determinant-1", farey_det1),
        ("same summand a+b", lambda p, q: p.summand == q.summand),
        ("same product ab", lambda p, q: p.product == q.product),
        ("same product shell mod27", lambda p, q: p.shell27 == q.shell27),
        ("same product sector mod7", lambda p, q: p.sector7 == q.sector7),
        ("bipartite subgraph comparable", comparable_bipartite_subgraph),
    ]
    rows = [graph_stats(packets, name, rel) for name, rel in relations]
    print("\nBinary-relation tournament on nonzero F_14 packets")
    print("relation                         edges components largest-components")
    for name, edges, comps, sizes in rows:
        print(f"{name:<32} {edges:>5} {comps:>10} {sizes}")


def main() -> None:
    level_table(N)

    f3 = farey_packets(3)
    f4 = farey_packets(4)
    p3 = {p.product for p in f3}
    p4 = {p.product for p in f4}
    new_f4 = sorted(set(f4) - set(f3), key=lambda p: p.value)

    print("\nF_3/F_4 product checks")
    print(f"F_3 packets: {[p.label() for p in f3]}")
    print(f"F_3 products: {sorted(p3)}; sum={sum(p3)}; divisors(6)={sorted(divisors(6))}")
    print(f"new F_4 packets: {[p.label() for p in new_f4]}")
    print(f"F_4 products: {sorted(p4)}; sum={sum(p4)}; divisors(12)={sorted(divisors(12))}")

    first_np = first_nonplanar_by_level(N)
    print("\nFirst complete-bipartite Kuratowski packet")
    if first_np is None:
        print("none through this limit")
    else:
        print(first_np.carrier())
        print(
            "checks: a+b=7 is the LRC14 apex; ab=12 is the old seed denominator; "
            "2(a+b)=14 and ab+2=14"
        )
        print(
            "K_{3,3} itself would be 3/3, but gcd(3,3)>1, so the reduced Farey "
            "frontier first sees it as K_{3,4}."
        )

    print("\nK_5 edge-count aliases inside F_14")
    for p in packets_with_product(N, 10):
        print(f"{p.carrier()}  edge-count alias for |E(K_5)|=10")
    print("These aliases remain planar; the ab model sees the K_{3,3} side first.")

    print("\nPackets sharing the product-12 / shell-12 neighborhood")
    product12 = packets_with_product(N, 12)
    shell12 = [p for p in farey_packets(N) if p.shell27 == 12]
    print(f"product=12: {[p.label() for p in product12]}")
    print(f"shell27=12: {[p.label() for p in shell12]}")

    relation_tournament(farey_packets(N))

    print("\nProof-search checksum")
    print("triple carrier: a/b -> (a+b, ab, K_{a,b})")
    print("3/4 -> (7, 12, K_{3,4}); this is the first reduced packet containing K_{3,3}.")
    print("Previous scout found product transform shell-octahedron witness (3,4,6,9,12,13).")
    print("So product 12 is simultaneously: F_3 product-sum, F_4 new obstruction, and a shell-octahedron vertex.")


if __name__ == "__main__":
    main()
