"""Finite conditional root kernels; literal proof transport, no assumed seeds.

Run normally or with -O. --write saves deterministic output. The graph is a
finite observation object, never a symbolic claim of universal coverage.
"""
from collections import deque
from dataclasses import dataclass
from pathlib import Path
import argparse
import json

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def odd(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'positive odd integer')


def step(n):
    odd(n)
    m, a = 3 * n + 1, 0
    while m % 2 == 0:
        m //= 2
        a += 1
    return m, a


def replay(n, word):
    odd(n)
    need(type(word) is tuple, 'valuation tuple')
    for a in word:
        need(n != 1, 'no first-hit root padding')
        need(type(a) is int and a > 0, 'positive exact valuation')
        n, actual = step(n)
        need(a == actual, 'literal valuation')
    return n


def discover(n, cap):
    """A bounded experiment, not a promised terminating certificate oracle."""
    odd(n)
    need(type(cap) is int and cap >= 0, 'finite discovery cap')
    word = ()
    for _ in range(cap):
        if n == 1:
            return word
        n, a = step(n)
        word += (a,)
    return word if n == 1 else None


@dataclass(frozen=True)
class Join:
    left: int
    right: int
    left_word: tuple
    right_word: tuple

    def validate(self):
        need(replay(self.left, self.left_word) == replay(self.right, self.right_word),
             'actual common endpoint')

    def transport(self, source, certificate):
        self.validate()
        need(source in (self.left, self.right), 'incident rooted source')
        need(replay(source, certificate) == 1, 'supplied seed certificate')
        if source == self.left:
            own, other, target = self.left_word, self.right_word, self.right
        else:
            own, other, target = self.right_word, self.left_word, self.left
        need(certificate[:len(own)] == own, 'deterministic shared prefix')
        result = other + certificate[len(own):]
        need(replay(target, result) == 1, 'transported literal root certificate')
        return target, result


class Kernel:
    def __init__(self, targets):
        self.targets = tuple(targets)
        for n in self.targets:
            odd(n)
        self.adjacency = {n: [] for n in (1, *self.targets)}
        self.joins = set()

    def insert(self, join):
        need(type(join) is Join, 'typed join')
        join.validate()
        if join in self.joins:
            return
        self.joins.add(join)
        for n in (join.left, join.right):
            self.adjacency.setdefault(n, []).append(join)

    def components(self):
        seen, result = set(), []
        for n in sorted(self.adjacency):
            if n in seen:
                continue
            todo, component = [n], set()
            while todo:
                u = todo.pop()
                if u in component:
                    continue
                component.add(u)
                for edge in self.adjacency[u]:
                    todo.extend((edge.left, edge.right))
            seen |= component
            result.append(component)
        return result

    def unresolved_kernel(self, grounded=(1,)):
        ground = set(grounded)
        targets = set(self.targets)
        return sorted(min(c) for c in self.components() if c & targets and not c & ground)

    def discharge(self, supplied):
        proven = {1: ()}
        for n, certificate in supplied.items():
            need(n in self.adjacency, 'seed belongs to recorded graph')
            need(replay(n, certificate) == 1, 'seed must carry an actual proof')
            proven[n] = certificate
        todo = deque(proven)
        while todo:
            u = todo.popleft()
            for edge in self.adjacency[u]:
                v = edge.right if u == edge.left else edge.left
                if v not in proven:
                    v, word = edge.transport(u, proven[u])
                    proven[v] = word
                    todo.append(v)
        return proven


def observed(bound, depth):
    graph = Kernel(range(1, bound + 1, 2))
    for n in graph.targets:
        for _ in range(depth):
            if n == 1:
                break
            m, a = step(n)
            graph.insert(Join(n, m, (a,), ()))
            n = m
    return graph


def reject(call):
    try:
        call()
    except ValueError:
        return
    raise ArithmeticError('invalid evidence was accepted')


def run():
    rows = []
    for bound in (31, 127, 255):
        graph = observed(bound, 4)
        seeds = graph.unresolved_kernel()
        need(all(s % 4 == 3 for s in seeds), 'minimal ungrounded representatives have valuation1')
        rooted_before = graph.discharge({})
        proofs = {s: discover(s, 10000) for s in seeds}
        need(all(w is not None for w in proofs.values()), 'declared finite seeds certified within cap')
        rooted = graph.discharge(proofs)
        need(all(n in rooted for n in graph.targets), 'all finite targets discharged')
        # Independent forward exploration checks the exported exact words.
        need(all(rooted[n] == discover(n, 10000) for n in graph.targets),
             'transported certificates equal independent first-hit routes')
        rows.append({'target_bound': bound, 'target_count': len(graph.targets),
                     'observed_edges': len(graph.joins), 'components': len(graph.components()),
                     'unrooted_targets_before': sum(n not in rooted_before for n in graph.targets),
                     'minimal_point_assumptions_for_fixed_graph': seeds,
                     'seed_odd_lengths': [len(proofs[s]) for s in seeds],
                     'certified_targets_after': len(graph.targets)})
    need(rows[0]['minimal_point_assumptions_for_fixed_graph'] == [27], 'small transparent kernel')
    need(rows[2]['minimal_point_assumptions_for_fixed_graph'] == [27, 111, 127, 159, 223, 231, 255],
         'frozen seed manifest')
    # New certificates reveal new edges and can merge formerly distinct components.
    graph = observed(255, 4)
    calls, new_edges = [], 0
    while graph.unresolved_kernel():
        seed = graph.unresolved_kernel()[0]
        word = discover(seed, 10000)
        need(word is not None, 'bounded adaptive seed certificate')
        calls.append(seed)
        n = seed
        before = len(graph.joins)
        for a in word:
            m, actual = step(n)
            need(a == actual, 'discovered route replay')
            graph.insert(Join(n, m, (a,), ()))
            n = m
        new_edges += len(graph.joins) - before
    rooted_final = graph.discharge({})
    need(all(n in rooted_final for n in graph.targets), 'adaptive graph grounded at1')
    circular = Kernel((3, 5))
    circular.insert(Join(3, 5, (1,), ()))
    circular.insert(Join(5, 3, (), (1,)))
    need(circular.unresolved_kernel() == [3] and set(circular.discharge({})) == {1},
         'an ungrounded equivalence cycle is not a proof')
    need(discover(27, 1) is None, 'finite discovery failure remains unresolved')
    reject(lambda: circular.discharge({3: ()}))
    reject(lambda: Join(3, 7, (1,), ()).validate())
    reject(lambda: replay(1, (2,)))
    reject(lambda: replay(True, ()))
    # A true descending schema can still lack family closure.
    need(9 % 4 == 1 and step(9)[0] == 7 and 7 % 4 != 1,
         'one-mod4 descent leaves its claimed family')
    return {'status': 'PASS', 'checks': CHECKS, 'finite_experiments': rows,
            'adaptive_seed_calls_for_bound255': calls, 'adaptive_new_edges': new_edges,
            'adaptive_final_edges': len(graph.joins),
            'scope': 'Exact finite assumption kernels and transported point certificates; infinite family laws require separate universal proofs. No universal Collatz coverage.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, sort_keys=True) + '\n'
    if args.write:
        root = Path(__file__).resolve().parents[2]
        (root / '05-knowledge/results/finite_seed_kernel_20261005.out').write_text(
            output, encoding='utf-8', newline='\n')
    print(output, end='')
