"""Exact one-vertex tournament Hamiltonian response as a weighted directed cut.

Standard library only. Vertices are 0,...,n-1; out[v] is an outgoing bitmask.
For S encoded by a bitmask, the new vertex points to precisely S.
"""

from collections import deque
from dataclasses import dataclass
from math import gcd


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(out):
    n = len(out)
    require(all(isinstance(x, int) and 0 <= x < 1 << n for x in out),
            "outgoing masks must be n-bit nonnegative integers")
    for u in range(n):
        require(not (out[u] >> u & 1), "loops are not allowed")
        for v in range(u + 1, n):
            require(((out[u] >> v) ^ (out[v] >> u)) & 1,
                    "each distinct pair must have exactly one orientation")


def endpoint_tables(out):
    """Subset DP ending and starting counts, including zeros for empty sets."""
    n = len(out)
    end = [[0] * n for _ in range(1 << n)]
    start = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        end[1 << v][v] = start[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            rest = mask ^ (1 << v)
            if rest:
                end[mask][v] = sum(end[rest][u] for u in range(n)
                                   if out[u] >> v & 1)
                start[mask][v] = sum(start[rest][u] for u in range(n)
                                     if out[v] >> u & 1)
    return end, start


@dataclass
class Response:
    out: tuple
    h: int
    starts: list
    ends: list
    covers: list
    weights: list
    degrees: list

    def capacity(self, cut):
        n = len(self.out)
        require(isinstance(cut, int) and 0 <= cut < 1 << n, "invalid cut")
        return sum(self.weights[u][v] for u in range(n) if cut >> u & 1
                   for v in range(n)
                   if self.out[u] >> v & 1 and not (cut >> v & 1))

    def value(self, cut):
        return self.h + 2 * self.capacity(cut)

    def all_values(self):
        values = [self.h] * (1 << len(self.out))
        for mask in range(1, len(values)):
            bit = mask & -mask
            v = bit.bit_length() - 1
            rest = mask ^ bit
            values[mask] = values[rest] + 2 * self.degrees[v] - 2 * sum(
                self.weights[u][v] for u in range(len(self.out)) if rest >> u & 1)
        return values

    def gain_content(self):
        """GCD of all extension gains, without enumerating the cuts."""
        divisor = 0
        for u in range(len(self.out)):
            for v in range(u + 1, len(self.out)):
                divisor = gcd(divisor, self.weights[u][v])
        return 2 * divisor

    def minimum_positive(self):
        """Minimum strictly increasing extension, or None if no such cut exists."""
        answer = None
        for u in range(len(self.out)):
            for v in range(len(self.out)):
                if self.out[u] >> v & 1:
                    candidate = self.minimum(1 << u, 1 << v)
                    require(candidate[0] > self.h, "positive-cut certificate failed")
                    if answer is None or candidate < answer:
                        answer = candidate
        return answer

    def minimum(self, inside=0, outside=0):
        """Minimum extension count subject to forced cut inclusions/exclusions.

        Returns (count, cut), using integral Edmonds-Karp on the original
        weighted tournament plus source/sink constraint arcs.
        """
        n = len(self.out)
        require(0 <= inside < 1 << n and 0 <= outside < 1 << n,
                "invalid forced masks")
        require(not (inside & outside), "contradictory forced memberships")
        source, sink = n, n + 1
        residual = [[0] * (n + 2) for _ in range(n + 2)]
        total = 0
        for u in range(n):
            for v in range(n):
                if self.out[u] >> v & 1:
                    residual[u][v] = self.weights[u][v]
                    total += self.weights[u][v]
        infinite = total + 1
        for u in range(n):
            if inside >> u & 1:
                residual[source][u] = infinite
            if outside >> u & 1:
                residual[u][sink] = infinite
        flow = 0
        while True:
            predecessor = [-1] * (n + 2)
            predecessor[source] = source
            queue = deque([source])
            while queue and predecessor[sink] < 0:
                u = queue.popleft()
                for v, capacity in enumerate(residual[u]):
                    if capacity and predecessor[v] < 0:
                        predecessor[v] = u
                        queue.append(v)
            if predecessor[sink] < 0:
                break
            delta = infinite
            v = sink
            while v != source:
                u = predecessor[v]
                delta = min(delta, residual[u][v])
                v = u
            v = sink
            while v != source:
                u = predecessor[v]
                residual[u][v] -= delta
                residual[v][u] += delta
                v = u
            flow += delta
        cut = sum(1 << u for u in range(n) if predecessor[u] >= 0)
        require((cut & inside) == inside and not (cut & outside),
                "flow failed to respect forced memberships")
        require(self.capacity(cut) == flow, "flow-cut certificate mismatch")
        return self.h + 2 * flow, cut


def compile_response(out):
    out = tuple(out)
    validate(out)
    n = len(out)
    full = (1 << n) - 1
    end, start = endpoint_tables(out)
    h = sum(end[full]) if n else 1
    covers = [[0] * n for _ in range(n)]
    for mask in range(1, full):
        rest = full ^ mask
        for u in range(n):
            if end[mask][u]:
                for v in range(n):
                    covers[u][v] += end[mask][u] * start[rest][v]
    weights = [[0] * n for _ in range(n)]
    for u in range(n):
        for v in range(u + 1, n):
            twice = covers[u][v] + covers[v][u]
            require(twice % 2 == 0, "nonintegral canonical weight")
            weights[u][v] = weights[v][u] = twice // 2
    degrees = [sum(weights[u][v] for v in range(n) if out[u] >> v & 1)
               for u in range(n)]
    return Response(out, h, start[full], end[full], covers, weights, degrees)


def extension(out, cut):
    n = len(out)
    require(0 <= cut < 1 << n, "invalid cut")
    return tuple(out[v] | (0 if cut >> v & 1 else 1 << n) for v in range(n)) + (cut,)


def from_code(n, code):
    require(0 <= code < 1 << (n * (n - 1) // 2), "invalid tournament code")
    out = [0] * n
    bit = 0
    for u in range(n):
        for v in range(u + 1, n):
            if code >> bit & 1:
                out[u] |= 1 << v
            else:
                out[v] |= 1 << u
            bit += 1
    return tuple(out)
