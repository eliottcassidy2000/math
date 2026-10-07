"""Marked completion plans consume authenticated ROOT words; they never find them.

A local paid common-future implication is not itself a grounded proof.  This
module retains OPEN and CYCLE defects through exact plan substitutions.
"""
from dataclasses import dataclass, replace
from itertools import product
from collections import Counter


def require(ok, message):
    if not ok:
        raise ValueError(message)


def odd(n):
    require(type(n) is int and n > 0 and n % 2 == 1, "exact positive odd source")


def word_type(word):
    require(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
            "exact tuple of positive valuations")


def replay(source, word):
    """Strict first-hit, constant-state-memory verifier; no orbit discovery."""
    odd(source)
    word_type(word)
    n = source
    for a in word:
        require(n != 1, "ROOT padding is not a first-hit certificate")
        z = 3*n+1
        actual = (z & -z).bit_length()-1
        require(a == actual, "supplied valuation does not match actual source")
        n = z >> a
    return n


@dataclass(frozen=True)
class Root:
    source: int
    word: tuple


@dataclass(frozen=True)
class Hole:
    source: int
    label: str


@dataclass(frozen=True)
class Transport:
    source: int
    child: int
    left: tuple
    right: tuple
    endpoint: int
    kind: str = "join"


def audit(node):
    require(type(node) in (Root, Hole, Transport), "exact proof-node type")
    odd(node.source)
    if type(node) is Root:
        require(replay(node.source, node.word) == 1, "supplied ROOT word required")
    elif type(node) is Hole:
        require(type(node.label) is str and bool(node.label), "nonempty defect label")
        require(node.source != 1, "ROOT cannot be relabelled as an assumption")
    else:
        odd(node.child)
        odd(node.endpoint)
        require(node.endpoint > 1, "an endpoint at ROOT is a direct Root certificate")
        require(type(node.kind) is str and node.kind in ("join", "forward", "paid"),
                "declared transport kind")
        require(replay(node.source, node.left) == node.endpoint
                and replay(node.child, node.right) == node.endpoint,
                "both actual paths must reach the same endpoint")
        require(bool(node.left or node.right), "nontrivial transport required")
        if node.kind == "paid":
            require(node.child < node.source, "payment compares the immutable source")
        if node.kind == "forward":
            require(node.right == () and node.child == node.endpoint and bool(node.left),
                    "forward rule has the exact observed terminal port")
    return node


def paid(source, child, left, right, endpoint):
    """Convert independently checked receipt fields, absorbing any ROOT endpoint."""
    odd(source)
    odd(child)
    odd(endpoint)
    require(child < source, "strict source payment")
    require(replay(child, right) == endpoint, "actual child side")
    require(replay(source, left) == endpoint, "actual source side")
    if endpoint == 1:
        return audit(Root(source, left))
    return audit(Transport(source, child, left, right, endpoint, "paid"))


def discharge(link, child_root):
    """Splice only the same child's authenticated suffix."""
    require(type(link) is Transport and type(child_root) is Root, "typed discharge")
    audit(link)
    audit(child_root)
    require(child_root.source == link.child, "child identity cannot be substituted")
    cut = len(link.right)
    require(child_root.word[:cut] == link.right, "same actual child prefix")
    return audit(Root(link.source, link.left+child_root.word[cut:]))


def compose(first, second):
    """Compress two actual implications while retaining the terminal obligation."""
    require(type(first) is Transport and type(second) is Transport, "typed transports")
    audit(first)
    audit(second)
    require(first.child == second.source, "composition keeps the exact shared port")
    b, c = first.right, second.left
    if len(b) <= len(c):
        require(c[:len(b)] == b, "shared deterministic prefix")
        left, right, endpoint = first.left+c[len(b):], second.right, second.endpoint
    else:
        require(b[:len(c)] == c, "shared deterministic prefix")
        left, right, endpoint = first.left, second.right+b[len(c):], first.endpoint
    # A two-edge implication may become a vacuous equality at the same source.
    # It is deliberately not promoted to a new useful transport or ROOT proof.
    if not left and not right:
        require(first.source == second.child, "empty composite source identity")
        return None
    return audit(Transport(first.source, second.child, left, right, endpoint))


@dataclass(frozen=True)
class Plan:
    nodes: tuple


def plan(nodes):
    require(type(nodes) is tuple, "tuple of marked nodes")
    for node in nodes:
        audit(node)
    require(len({x.source for x in nodes}) == len(nodes), "one selected rule per source")
    table = {x.source: x for x in nodes}
    require(table.get(1) == Root(1, ()), "literal ROOT boundary is retained")
    for node in nodes:
        if type(node) is Transport:
            require(node.child in table, "every dependency has an explicit port")
    return Plan(tuple(sorted(nodes, key=lambda x: x.source)))


def audit_plan(packet):
    require(type(packet) is Plan, "exact Plan")
    require(plan(packet.nodes) == packet, "canonical marked plan")
    return {node.source: node for node in packet.nodes}


@dataclass(frozen=True)
class Trace:
    source: int
    status: str
    path: tuple
    terminal: int
    cycle: tuple = ()


def trace(packet, source):
    table = audit_plan(packet)
    odd(source)
    require(source in table, "original requested source must have a port")
    path, seen = [], {}
    current = source
    while current not in seen:
        seen[current] = len(path)
        path.append(current)
        node = table[current]
        if type(node) is Root:
            return Trace(source, "GROUNDED", tuple(path), current)
        if type(node) is Hole:
            return Trace(source, "OPEN", tuple(path), current)
        current = node.child
    return Trace(source, "CYCLE", tuple(path), current, tuple(path[seen[current]:]))


def close(packet, source):
    """A total bounded consumer: return Root or an explicit Trace defect."""
    result = trace(packet, source)
    if result.status != "GROUNDED":
        return result
    table = {node.source: node for node in packet.nodes}
    certificate = table[result.terminal]
    for n in reversed(result.path[:-1]):
        certificate = discharge(table[n], certificate)
    require(certificate.source == source, "export preserves immutable source")
    return certificate


def rooted_clock(packet, source):
    """Exact odd first-hit length from an authenticated leaf and signed edge costs."""
    result = trace(packet, source)
    require(result.status == "GROUNDED", "an open clock is not a ROOT deadline")
    table = {node.source: node for node in packet.nodes}
    return len(table[result.terminal].word)+sum(
        len(table[n].left)-len(table[n].right) for n in result.path[:-1])


def conditional_link(packet, source):
    """The checked source-to-hole implication, never a ROOT certificate."""
    result = trace(packet, source)
    require(result.status == "OPEN", "only an open path has an unresolved leaf")
    table = {node.source: node for node in packet.nodes}
    combined = None
    for n in result.path[:-1]:
        combined = table[n] if combined is None else compose(combined, table[n])
    return combined, table[result.terminal]


def refine(packet, source, replacement_nodes):
    """Replace exactly one hole; keep colliding ports equal or reject.

    Newly referring to an old port is allowed, so refinement can expose a
    cycle.  The returned packet must be traced again before ROOT export.
    """
    table = audit_plan(packet)
    odd(source)
    require(source in table and type(table[source]) is Hole, "replace a declared open port")
    require(type(replacement_nodes) is tuple, "tuple of replacement nodes")
    incoming = {}
    for node in replacement_nodes:
        audit(node)
        require(node.source not in incoming, "no duplicate replacement ports")
        incoming[node.source] = node
    require(source in incoming, "replacement retains the exact hole source")
    del table[source]
    for n, node in incoming.items():
        require(n not in table or table[n] == node, "cannot overwrite another retained port")
        table[n] = node
    return plan(tuple(table.values()))


def provide(packet, source, word):
    odd(source)
    return refine(packet, source, (audit(Root(source, word)),))


def supply_root(packet, source, word):
    """Install a checked anchor even at a cyclic port; the input plan is immutable."""
    table = audit_plan(packet)
    odd(source)
    require(source in table, "a supplied anchor keeps an existing exact source port")
    certificate = audit(Root(source, word))
    table[source] = certificate
    return plan(tuple(table.values()))


def ranks(packet):
    """Least selected-proof heights, with no value for open/cyclic ports."""
    table = audit_plan(packet)
    result = {n: 0 for n, node in table.items() if type(node) is Root}
    changed = True
    while changed:
        changed = False
        for n, node in table.items():
            if n not in result and type(node) is Transport and node.child in result:
                result[n] = result[node.child]+1
                changed = True
    return result


def verify_ranks(packet, rank):
    table = audit_plan(packet)
    require(type(rank) is dict and all(type(n) is int for n in rank)
            and set(rank) == set(table), "complete exact labelled rank")
    require(all(type(r) is int and r >= 0 for r in rank.values()), "natural proof ranks")
    for n, node in table.items():
        require(type(node) is not Hole, "an assigned rank does not authenticate a hole")
        if type(node) is Transport:
            require(rank[node.child] < rank[n], "strict dependency-rank decrease")
    return True


def defect_cut(packet):
    """A selected-bank non-grounding witness; never a divergence certificate."""
    table = audit_plan(packet)
    grounded = ranks(packet)
    return frozenset(set(table)-set(grounded))


def verify_defect_cut(packet, cut):
    table = audit_plan(packet)
    require(type(cut) is frozenset and all(type(n) is int for n in cut)
            and cut <= set(table), "exact selected-bank cut")
    for n in cut:
        node = table[n]
        require(type(node) is not Root, "a cut cannot hide an authenticated ROOT seed")
        if type(node) is Transport:
            require(node.child in cut, "all selected dependencies remain behind the cut")
    return True


def main():
    checks = Counter()

    def check(ok, key):
        checks[key] += 1
        if not ok:
            raise ArithmeticError(f"{key}: {checks[key]}")

    def rejects(fn):
        try:
            fn()
        except (ValueError, TypeError):
            check(True, "malformed")
        else:
            check(False, "malformed")

    root_words = {3: (1, 4), 5: (4,), 7: (1, 1, 2, 3, 4)}
    to_five = {3: (1,), 5: (), 7: (1, 1, 2, 3)}
    universes = []
    for n in (3, 5, 7):
        choices = [Hole(n, f"ROOT({n})"), Root(n, root_words[n])]
        choices += [Transport(n, m, to_five[n], to_five[m], 5,
                              "paid" if m < n else "join") for m in (3, 5, 7) if m != n]
        universes.append(choices)
    totals = Counter()
    for selected in product(*universes):
        packet = plan((Root(1, ()),)+selected)
        r = ranks(packet)
        cut = defect_cut(packet)
        check(verify_defect_cut(packet, cut), "plans")
        for source in (1, 3, 5, 7):
            result = trace(packet, source)
            totals[result.status] += 1
            check((source in r) == (result.status == "GROUNDED"), "plans")
            certificate = close(packet, source)
            if source in r:
                check(type(certificate) is Root and certificate.source == source
                      and replay(source, certificate.word) == 1, "exports")
                check(rooted_clock(packet, source) == len(certificate.word), "clocks")
            else:
                check(type(certificate) is Trace, "exports")
            if result.status == "OPEN":
                link, hole = conditional_link(packet, source)
                check(hole.source == result.terminal, "conditional")
                if link is not None:
                    check(link.source == source and link.child == hole.source, "conditional")
                    grounded = discharge(link, Root(hole.source, root_words[hole.source]))
                    check(replay(source, grounded.word) == 1, "conditional")
        if not cut:
            check(verify_ranks(packet, r), "plans")
        for node in selected:
            if type(node) is Hole:
                patched = provide(packet, node.source, root_words[node.source])
                check(set(r) <= set(ranks(patched)), "monotone")
                for source in r:
                    check(close(packet, source) == close(patched, source), "monotone")

    initial = plan((Root(1, ()), Transport(3, 5, (1,), (), 5, "forward"),
                    Hole(5, "unproved completion edge to ROOT")))
    check(trace(initial, 3).status == "OPEN", "cycle_hostile")
    cyclic = refine(initial, 5, (paid(5, 3, (), (1,), 5),))
    check(trace(cyclic, 3).cycle == (3, 5), "cycle_hostile")
    check(type(close(cyclic, 3)) is Trace and defect_cut(cyclic) == frozenset((3, 5)),
          "cycle_hostile")
    check(compose(cyclic.nodes[1], cyclic.nodes[2]) is not None
          and compose(cyclic.nodes[2], cyclic.nodes[1]) is None, "cycle_hostile")
    rejects(lambda: rooted_clock(cyclic, 3))
    rejects(lambda: verify_ranks(cyclic, {1: 0, 3: 1, 5: 2}))
    complete = provide(initial, 5, (4,))
    check(close(complete, 3) == Root(3, (1, 4)), "cycle_hostile")
    check(verify_ranks(complete, ranks(complete)), "cycle_hostile")
    repaired_cycle = supply_root(cyclic, 5, (4,))
    check(close(repaired_cycle, 3) == Root(3, (1, 4))
          and trace(cyclic, 3).status == "CYCLE", "cycle_hostile")

    # The bounded head compiler contributes arithmetic, not its missing child proof.
    import collatz_twoanchor_head_decoder_20261007b as heads
    receipt = heads.native_receipt(heads.compile_partner((10,)), 0)
    link = paid(receipt.source, receipt.child, receipt.source_word,
                receipt.child_word, receipt.endpoint)
    head_plan = plan((Root(1, ()), link, Hole(receipt.child, "supplied native child ROOT")))
    check(trace(head_plan, 63829).status == "OPEN", "native_head")
    supplied = (3,1,1,3,3,2,1,3,1,1,3,4,1,3,1,2,3,4)
    fixed = provide(head_plan, 2365, supplied)
    actual = close(fixed, 63829)
    check(type(actual) is Root and len(actual.word) == 15, "native_head")
    check(actual.word == heads.discharge(heads.compile_partner((10,)), 0, supplied), "native_head")
    rejects(lambda: provide(head_plan, 63829, actual.word))
    rejects(lambda: discharge(link, Root(3, (1, 4))))

    # One newly supplied J2 child is consumed; discovery is deliberately not called.
    import collatz_two_twos_grounding_20261007c as supplied_j2
    large_receipt = supplied_j2.parent_receipt(6129, (1, 10), (1, 1, 2, 2, 3), 1)
    large_link = paid(large_receipt.source, large_receipt.child,
                      large_receipt.source_word, large_receipt.child_word,
                      large_receipt.endpoint)
    large_open = plan((Root(1, ()), large_link,
                       Hole(large_receipt.child, "frozen J2 child M6125")))
    check(trace(large_open, large_receipt.source).status == "OPEN", "J2_integration")
    child_word = supplied_j2.root_word(6125)
    large_closed = provide(large_open, large_receipt.child, child_word)
    parent_root = close(large_closed, large_receipt.source)
    check(type(parent_root) is Root and parent_root.source == (1 << 6129)-1,
          "J2_integration")
    check((len(child_word), sum(child_word)) == (28762, 51712), "J2_integration")
    check((len(parent_root.word), sum(parent_root.word)) == (28762, 51716),
          "J2_integration")
    check(parent_root.word == supplied_j2.discharge(large_receipt, child_word),
          "J2_integration")
    check(rooted_clock(large_closed, large_receipt.source) == len(parent_root.word),
          "J2_integration")
    check(trace(large_open, large_receipt.source).status == "OPEN", "J2_integration")

    # A direct ROOT endpoint is absorbed, so a known terminal is never hidden.
    check(paid(5, 1, (4,), (), 1) == Root(5, (4,)), "root_boundary")
    rejects(lambda: audit(Transport(5, 1, (4,), (), 1, "paid")))
    rejects(lambda: audit(Root(1, (2,))))
    rejects(lambda: provide(initial, 5, (4, 2)))
    for bad in (True, 3.0, -3, 2):
        rejects(lambda bad=bad: trace(initial, bad))
        rejects(lambda bad=bad: audit(Hole(bad, "bad")))
    rejects(lambda: plan((Root(1, ()), Hole(3, "a"), Hole(3, "b"))))
    rejects(lambda: plan((Root(1, ()), Transport(3, 5, (1,), (), 5, "forward"))))
    rejects(lambda: refine(initial, 5, (Root(7, root_words[7]),)))
    rejects(lambda: refine(initial, 5, (Root(5, (4,)), Root(3, (1, 4)))))
    rejects(lambda: verify_defect_cut(cyclic, frozenset((True, 3, 5))))
    rejects(lambda: verify_ranks(complete, {1: False, 3: 1, 5: 0}))
    rejects(lambda: paid(3, 5, (1,), (), 5))
    rejects(lambda: supply_root(cyclic, 5, root_words[3]))
    for bad_endpoint in (True, 1.0):
        rejects(lambda bad_endpoint=bad_endpoint: paid(5, 1, (4,), (), bad_endpoint))

    print("Selected-plan universe: 64 choices on ports 3, 5, 7, plus literal ROOT")
    print("Port traces:", dict(sorted(totals.items())))
    print("Local-paid substitution hostile: 3 -> 5 -> 3 remains CYCLE")
    print("Supplied ROOT(5) repairs the open plan; exported ROOT(3) word=(1, 4)")
    print("Native head: source63829 / child2365 / join281; supplied child proof exports15 edges")
    print("Frozen J2 integration: M6125 ROOT rank28762/cost51712 closes M6129 rank28762/cost51716")
    print("The immutable earlier M6125 port remains OPEN in its original plan; no orbit discovery")
    print("Checks:", dict(sorted(checks.items())))
    print("Total checks:", sum(checks.values()))


if __name__ == "__main__":
    main()
