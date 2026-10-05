"""Fair supplied-source Collatz extension; guarded words are not free literal edges.

Only ROOT1 is initially grounded. Search never calls a home oracle. Post-search
certificate audits are separate, counted validation work. Explicit checks survive -O.
"""
from collections import Counter, deque
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import argparse
import json

import inverse_ray_ternary_addresses_20261004 as codec
import debt_word_join_search_20261004 as debt
import collatz_binary_ternary_guard_fusion_20261004 as fusion


def need(ok, message):
    if not ok:
        raise ValueError(message)


def source_guard(n):
    need(type(n) is int and n > 0 and n % 2, 'positive exact odd supplied source')


def step(n):
    source_guard(n)
    raw = 3*n+1
    a = (raw & -raw).bit_length()-1
    return raw >> a, a


def literal(n):
    source_guard(n)
    z, a = 3*n+1, 0
    while z % 2 == 0:
        z //= 2
        a += 1
    return z, a


@dataclass(frozen=True)
class Arc:
    source: int
    target: int
    word: tuple
    origin: str

    def validate(self):
        source_guard(self.source)
        source_guard(self.target)
        need(self.source != 1, 'no root self-arc')
        need(self.word and all(type(a) is int and a > 0 for a in self.word), 'positive exact word')
        P, Q, B = debt.data(self.word)
        need(P*self.source+B == Q*self.target, 'retained source and endpoint affine identity')
        need((P*self.source+B) % (2*Q) == Q, 'exact final oddness source guard')


class ProofGraph:
    """Rooted least closure of guarded actual paths, in either direction."""
    def __init__(self):
        self.arcs = []
        self.keys = set()
        self.adjacency = {1: []}
        self.outgoing = {}
        self.parent = {1: 1}
        self.size = {1: 1}

    def vertex(self, n):
        if n not in self.parent:
            self.parent[n] = n
            self.size[n] = 1
            self.adjacency[n] = []

    def find(self, n):
        self.vertex(n)
        while self.parent[n] != n:
            self.parent[n] = self.parent[self.parent[n]]
            n = self.parent[n]
        return n

    def grounded(self, n):
        return self.find(n) == self.find(1)

    def add(self, arc):
        arc.validate()
        key = (arc.source, arc.target, arc.word)
        if key in self.keys:
            return False
        self.keys.add(key)
        index = len(self.arcs)
        self.arcs.append(arc)
        for n in (arc.source, arc.target):
            self.vertex(n)
            self.adjacency[n].append(index)
        self.outgoing.setdefault(arc.source, []).append(index)
        a, b = self.find(arc.source), self.find(arc.target)
        if a != b:
            if self.size[a] < self.size[b]:
                a, b = b, a
            self.parent[b] = a
            self.size[a] += self.size[b]
        return True

    def route_arc(self, n):
        indices = self.outgoing.get(n, ())
        if not indices:
            return None
        return self.arcs[min(indices, key=lambda i: (-len(self.arcs[i].word), i))]


LONG_WORD = (1,1,1,2,1,1,1,1,2,1,1,2,1,1,3,2,3,1,1,2,1,1,1,5,1,1,1,1,4,2,6)
COARSE_WORDS = ((2,), LONG_WORD, (1,1,2,1,2,1,4), (1,1,1,1,2,1,1,5))


def load_bank():
    root = Path(__file__).resolve().parents[2]
    record = json.loads((root/'05-knowledge/results/collatz_binary_ternary_guard_fusion_20261004.json').read_text())
    bank = {(row['s'], row['e']): (tuple(row['source_word']), tuple(row['child_word']))
            for row in record['debt_rows']}
    need(len(bank) == 16, 'frozen sixteen-row debt bank')
    for (s, e), (u, v) in bank.items():
        P, Q, B = debt.data(u)
        p, q, b = debt.data(v)
        need(p == P*3**e and Q == q*2**s and P+B == b*2**s, 'bank affine theorem premise')
    return bank


def coarse_match(n, word):
    P, Q, B = debt.data(word)
    if (P*n+B) % Q:
        return None
    nominal = (P*n+B)//Q
    h = (nominal & -nominal).bit_length()-1
    return n, nominal >> h, tuple(word[:-1])+(word[-1]+h,)


def rule_candidate(n, index, bank):
    """One finite arithmetic guard; returns actual arcs and an optional child.

    No exact edge replay or home search occurs here. The caller validates every
    returned source/word/endpoint triple before it can influence closure.
    """
    if 2 <= index < 6:
        result = coarse_match(n, COARSE_WORDS[index-2])
        return ([] if result is None else [result]), None
    if index == 1:
        r = ((n+1) & -(n+1)).bit_length()-2
        if r < 1:
            return [], None
        z = 3**r*(n+1)//2**r-1
        endpoint, a = step(z)
        if a < 3:
            return [], None
        m = (n-1)//2
        return [(n, endpoint, (1,)*r+(a,)),
                (m, endpoint, (1,)*(r-1)+(2,a-2))], m
    if index == 0:
        hit = fusion.select_debt(n, bank)
        if hit is None:
            return [], None
        return [(n, hit['endpoint'], tuple(hit['parent_word'])),
                (hit['child'], hit['endpoint'], tuple(hit['child_word']))], hit['child']
    k = index-6
    row = fusion.sibling_row(k)
    if n % row['modulus'] != row['residue']:
        return [], None
    S = 4**k*n+(4**k-1)//3
    m = 2**row['ell']*((S+1)//row['modulus'])-1
    need(0 < m < n and m % 2, 'guarded sibling smaller child')
    endpoint, a = step(n)
    return [(n, endpoint, (a,)), (m, endpoint, (1,)*row['ell']+(a+2*k,))], m


class FairSearch:
    """Each round reserves an original-request turn before one symbolic action."""
    def __init__(self, requested, rules=True, busy=False):
        self.requested = tuple(requested)
        need(bool(self.requested) and len(set(self.requested)) == len(self.requested), 'nonempty distinct finite requests')
        for n in self.requested:
            source_guard(n)
        need(type(rules) is bool and type(busy) is bool, 'exact policy booleans')
        self.rules = rules
        self.bank = load_bank() if rules else {}
        self.graph = ProofGraph()
        self.cursors = list(self.requested)
        self.next_request = 0
        self.queue = deque()
        self.probed = set()
        self.children = set()
        self.stats = Counter()
        self.discovery_trace = []
        for n in self.requested:
            self.register(n)
        if busy:
            self.queue.append(('busy', 0, 0))

    def register(self, n):
        self.graph.vertex(n)
        if self.rules and n != 1 and n not in self.probed:
            self.probed.add(n)
            self.queue.append(('probe', n, 0))

    def done(self):
        return all(self.graph.grounded(n) for n in self.requested)

    def add(self, source, target, word, origin):
        if source == 1:
            need(target == 1 and all(a == 2 for a in word), 'formal root tail is checked but omitted')
            self.stats['omitted_root_tails'] += 1
            return
        new = self.graph.add(Arc(source, target, tuple(word), origin))
        self.stats[origin+'_arc_submissions'] += 1
        if new:
            self.stats[origin+'_new_arcs'] += 1
            self.stats[origin+'_stored_word_letters'] += len(word)
        self.register(source)
        self.register(target)

    def advance_cursor(self, current, origin):
        arc = self.graph.route_arc(current)
        if arc is not None:
            self.stats[origin+'_cached_arcs'] += 1
            self.stats[origin+'_cached_word_letters'] += len(arc.word)
            return arc.target
        target, a = step(current)
        need(literal(current) == (target, a), 'independent literal observation validation')
        self.stats[origin+'_literal_queries'] += 1
        self.stats['query_validation_steps'] += 1
        self.discovery_trace.append((current, target, a, origin))
        self.add(current, target, (a,), 'literal')
        return target

    def original_turn(self):
        N = len(self.requested)
        for _ in range(N):
            i = self.next_request
            self.next_request = (i+1) % N
            if self.graph.grounded(self.requested[i]):
                continue
            need(not self.graph.grounded(self.cursors[i]), 'request cursor path remains retained')
            self.cursors[i] = self.advance_cursor(self.cursors[i], 'original')
            self.stats['original_turns'] += 1
            return

    def symbolic_turn(self):
        if not self.queue:
            return
        self.stats['symbolic_turns'] += 1
        # Reserve alternating auxiliary actions for already-guarded dependencies.
        # Other actions retain FIFO progress, so even this priority cannot starve probes.
        priority = next((i for i, job in enumerate(self.queue) if job[0] == 'walk'), None)
        if self.stats['symbolic_turns'] % 2 == 0 and priority is not None:
            kind, n, index = self.queue[priority]
            del self.queue[priority]
        else:
            kind, n, index = self.queue.popleft()
        if kind == 'busy':
            self.stats['unfinished_decoder_steps'] += 1
            self.queue.append((kind, n, index+1))
        elif kind == 'walk':
            if not self.graph.grounded(n):
                target = self.advance_cursor(n, 'auxiliary')
                self.queue.append(('walk', target, 0))
        else:
            need(kind == 'probe', 'known resumable job kind')
            if self.graph.grounded(n):
                self.stats['grounded_probe_skips'] += 1
                return
            self.stats['guard_tests'] += 1
            arcs, child = rule_candidate(n, index, self.bank)
            if arcs:
                self.stats['guard_hits'] += 1
                for source, target, word in arcs:
                    self.add(source, target, word, 'guard')
                if child is not None and child != 1 and child not in self.children:
                    self.children.add(child)
                    self.queue.append(('walk', child, 0))
            # One source-aware sibling index per action. Positivity gives a finite cap.
            if index < 6+(n+1).bit_length()//2:
                self.queue.append(('probe', n, index+1))

    def run(self, rounds):
        need(type(rounds) is int and rounds >= 0, 'nonnegative exact round budget')
        for _ in range(rounds):
            if self.done():
                break
            self.original_turn()
            self.symbolic_turn()
            self.stats['rounds'] += 1
        return 'COMPLETE' if self.done() else 'PENDING'

    def snapshot(self):
        return dict(requested=self.requested, rules=self.rules, cursors=self.cursors,
                    next_request=self.next_request, queue=list(self.queue), probed=sorted(self.probed),
                    children=sorted(self.children), stats=dict(self.stats), trace=self.discovery_trace,
                    arcs=[(a.source,a.target,a.word,a.origin) for a in self.graph.arcs])

    @classmethod
    def restore(cls, state):
        result = cls(state['requested'], state['rules'])
        result.graph = ProofGraph()
        for source, target, word, origin in state['arcs']:
            result.graph.add(Arc(source,target,tuple(word),origin))
        result.cursors = list(state['cursors'])
        result.next_request = state['next_request']
        result.queue = deque(tuple(x) for x in state['queue'])
        result.probed = set(state['probed'])
        result.children = set(state['children'])
        result.stats = Counter(state['stats'])
        result.discovery_trace = [tuple(x) for x in state['trace']]
        for n in result.requested:
            result.graph.vertex(n)
        need(len(result.cursors) == len(result.requested) and type(result.next_request) is int
             and 0 <= result.next_request < len(result.requested), 'retained fixed request schedule')
        for source, cursor in zip(result.requested, result.cursors):
            source_guard(cursor)
            seen, todo = {source}, [source]
            while todo and cursor not in seen:
                for index in result.graph.outgoing.get(todo.pop(), ()):
                    target = result.graph.arcs[index].target
                    if target not in seen:
                        seen.add(target)
                        todo.append(target)
            need(cursor in seen, 'restored cursor has an actual retained path from its original source')
        for kind, n, index in result.queue:
            need(kind in ('probe','walk','busy') and type(index) is int and index >= 0,
                 'resumable finite job descriptor')
            if kind != 'busy':
                source_guard(n)
            if kind == 'probe':
                need(result.rules and index <= 6+(n+1).bit_length()//2, 'source-aware probe stage')
        need(all(type(value) is int and value >= 0 for value in result.stats.values()), 'nonnegative retained counters')
        return result

    def report(self):
        grounded = sum(self.graph.grounded(n) for n in self.requested)
        return dict(status='COMPLETE' if grounded == len(self.requested) else 'PENDING',
                    requested=len(self.requested), grounded=grounded, arcs=len(self.graph.arcs),
                    vertices=len(self.graph.parent), pending_jobs=len(self.queue), **dict(sorted(self.stats.items())))


def certificates(graph):
    """Construct canonical first-hit ASTs by a ROOT-rooted spanning tree.

    Backward traversal prepends inverses; forward traversal cuts an existing
    certificate. Root padding is checked then discarded, never stored.
    """
    certs, queue, count = {1: codec.ROOT}, deque([1]), 0
    while queue:
        known = queue.popleft()
        for index in graph.adjacency[known]:
            arc = graph.arcs[index]
            other = arc.target if known == arc.source else arc.source
            if other in certs:
                continue
            cert = certs[known]
            if known == arc.source:
                for a in arc.word:
                    count += 1
                    if cert == codec.ROOT:
                        need(a == 2, 'only formal root padding after first hit')
                    else:
                        need(codec.exponent(cert) == a, 'forward cut agrees with supplied first-hit route')
                        cert = cert.parent
            else:
                value = known
                for a in reversed(arc.word):
                    count += 1
                    numerator = (value << a)-1
                    need(numerator % 3 == 0, 'integral checked inverse')
                    source = numerator//3
                    if source == 1:
                        need(value == 1 and a == 2 and cert == codec.ROOT, 'only root self-return is discarded')
                    else:
                        least = codec.kappa(value % 9, source % 3)
                        need(a >= least and (a-least) % 6 == 0, 'inverse address agrees with exact edge')
                        cert = codec.extend(cert, source % 3, (a-least)//6)
                    value = source
                need(value == other, 'retained original source label')
            need(codec.expand(cert,bit_cap=codec.bit_bounds(cert)[1]) == other,
                 'spanning-tree certificate source identity with the codec conservative bound')
            certs[other] = cert
            queue.append(other)
    need(set(certs) == {n for n in graph.parent if graph.grounded(n)}, 'least grounded component only')
    return certs, count


def audit_search(search):
    letters = 0
    for arc in search.graph.arcs:
        x = arc.source
        for a in arc.word:
            x, actual = literal(x)
            need(actual == a, 'independent post-search literal macro validation')
            letters += 1
        need(x == arc.target, 'literal macro endpoint validation')
    certs, construction = certificates(search.graph)
    odd_edges = 0
    for n in search.requested:
        if n in certs:
            independent = codec.encode_source(n, step_cap=10000)
            need(certs[n] == independent, 'independent canonical first-hit certificate')
            odd_edges += codec.ranks(independent)[0]
    return dict(post_search_arc_letter_replays=letters, certificate_splice_letters=construction,
                independent_requested_route_edges=odd_edges)


def generic_root_component(pairs):
    adjacency = {}
    for a, b in pairs:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    seen, todo = {1}, [1]
    while todo:
        for target in adjacency.get(todo.pop(), ()):
            if target not in seen:
                seen.add(target)
                todo.append(target)
    return seen


def export_certificates(search, destination):
    """Export retained proofs only after search; reload and strictly replay them."""
    need(search.done() and search.rules, 'completed guarded search supplies the export')
    digest = sha256(json.dumps(search.requested,separators=(',',':')).encode()).hexdigest()
    need(len(search.requested) == 239 and digest ==
         'c1a01b7d890f9bcf6205572a65f112b4038c290b67c10ef86e5d1eba40df31b7',
         'export exactly the frozen239 requested sources')
    certs, _ = certificates(search.graph)
    rows = []
    for source in search.requested:
        word = tuple(codec.exponent(node) for node in codec.chain(certs[source]))
        rows.append(dict(source=source,valuation_word=word,odd_rank=len(word),
                         halving_cost=sum(word),ordinary_rank=len(word)+sum(word)))
    payload = dict(status='VERIFIED finite first-hit certificates; universal Collatz coverage OPEN',
                   universe='Exactly the ordered frozen239 inherited seed sources; ROOT1-only guarded search',
                   source_count=239,source_sha256=digest,root=1,
                   proof_provenance='Output of completed fair search; never supplied as search inputs',
                   search_counters=search.report(),certificates=rows)
    text = json.dumps(payload,indent=2)+'\n'
    # Validate the serialized representation, not just the in-memory ASTs.
    decoded = json.loads(text)
    need(tuple(row['source'] for row in decoded['certificates']) == search.requested,
         'export preserves exactly the requested source universe and order')
    for row in decoded['certificates']:
        value = row['source']
        word = row['valuation_word']
        need(row['odd_rank'] == len(word) and row['halving_cost'] == sum(word) and
             row['ordinary_rank'] == len(word)+sum(word), 'exported ranks')
        for exponent in word:
            need(value != 1, 'export has no padded root self-edge')
            value, actual = literal(value)
            need(actual == exponent, 'serialized valuation is exact')
        need(value == 1, 'serialized first-hit route terminates at ROOT')
    destination = Path(destination)
    destination.write_text(text,encoding='utf-8')
    need(json.loads(destination.read_text(encoding='utf-8')) == decoded, 'written certificate round-trip')


def main(export_path=None):
    root = Path(__file__).resolve().parents[2]
    record = json.loads((root/'05-knowledge/results/checked_switch_phase19_20261004.json').read_text())
    seeds = tuple(record['compiler'][-1]['seed_sources'])
    digest = sha256(json.dumps(seeds,separators=(',',':')).encode()).hexdigest()
    need(len(seeds) == len(set(seeds)) == 239 and digest ==
         'c1a01b7d890f9bcf6205572a65f112b4038c290b67c10ef86e5d1eba40df31b7', 'frozen 239 supplied requests')
    print('FAIR FRONTIER EXTENSION: supplied-state fairness; ROOT1 only; universal termination OPEN')
    print('Frozen239 seed hash:', digest)
    print('Each round: at most one original cursor action, then one resumable symbolic action.')
    print('Literal queries, source guards, stored macro letters and post-search validation are separate costs.')
    runs = {}
    for rules in (False, True):
        search = FairSearch(seeds, rules)
        previous = 0
        for cap in (239, 1000, 2000, 4000, 10000):
            search.run(cap-previous)
            previous = cap
            print('CHECKPOINT', 'guarded' if rules else 'literal', cap, json.dumps(search.report(),sort_keys=True))
            if search.done():
                break
        need(search.done(), 'declared finite benchmark cap10000')
        print('AUDIT', 'guarded' if rules else 'literal', json.dumps(audit_search(search),sort_keys=True))
        runs[rules] = search
    uninterrupted = FairSearch(seeds, True)
    uninterrupted.run(10000)
    resumed = FairSearch(seeds, True)
    need(resumed.run(0) == 'PENDING' and not resumed.graph.arcs, 'zero budget leaves evidence unchanged')
    resumed.run(137)
    frozen = json.dumps(resumed.snapshot(),sort_keys=True)
    resumed = FairSearch.restore(json.loads(frozen))
    resumed.run(863)
    resumed = FairSearch.restore(json.loads(json.dumps(resumed.snapshot())))
    resumed.run(9000)
    need(resumed.snapshot() == uninterrupted.snapshot(), 'budget partitions and JSON resumption preserve exact state')
    print('RESUME exact equality:0+137+863+9000 rounds versus10000; snapshot bytes',len(frozen))
    bank_record = json.loads((root/'05-knowledge/results/collatz_binary_ternary_guard_fusion_20261004.json').read_text())
    heldout = tuple(row['first_tested_example']['source']+
                    2**(sum(row['first_tested_example']['parent_word'])+1)
                    for row in bank_record['debt_rows'])
    need(len(set(heldout)) == 16 and not set(heldout) & set(seeds), 'separate prespecified positive-lift controls')
    print('HELDOUT16_SOURCES',json.dumps(heldout))
    for rules in (False, True):
        search = FairSearch(heldout,rules)
        need(search.run(10000) == 'COMPLETE', 'declared held-out cap10000')
        print('HELDOUT16','guarded' if rules else 'literal',json.dumps(search.report(),sort_keys=True))
        print('HELDOUT_AUDIT','guarded' if rules else 'literal',json.dumps(audit_search(search),sort_keys=True))
    hostile = FairSearch((27,703), rules=False, busy=True)
    hostile.run(1000)
    need(hostile.done() and hostile.stats['unfinished_decoder_steps'] == hostile.stats['rounds'],
         'endlessly unfinished symbolic branch cannot starve supplied literal requests')
    print('UNFINISHED_DECODER',json.dumps(hostile.report(),sort_keys=True))
    singleton = FairSearch((27,),False)
    need(singleton.run(1000) == 'COMPLETE', 'singleton27 completion regression')
    need(certificates(singleton.graph)[0][27] == codec.encode_source(27),
         'certificate export uses codec bound rather than endpoint bit length')
    # Root-free cycles are tested generically, without asserting an unknown positive Collatz cycle.
    need(generic_root_component(((1,9),(3,5),(5,3))) == {1,9}, 'a formal unsupported cycle cannot certify itself')
    need(generic_root_component((a.source,a.target) for a in runs[True].graph.arcs) ==
         {n for n in runs[True].graph.parent if runs[True].graph.grounded(n)}, 'independent least-component audit')
    disconnected = ProofGraph()
    disconnected.add(Arc(5,1,(4,),'control'))
    disconnected.add(Arc(7,11,(1,),'control'))
    need(disconnected.grounded(5) and not disconnected.grounded(7), 'ungrounded observations are not proofs')
    need(not FairSearch((7,),False).graph.grounded(7), 'unobserved supplied source is not a seed')
    for invalid in (True, 1.0, 0, -1, 2):
        try:
            FairSearch((invalid,),False)
        except ValueError:
            pass
        else:
            raise ValueError('malformed supplied state accepted')
    try:
        Arc(27,1,(1,),'hostile').validate()
    except ValueError:
        pass
    else:
        raise ValueError('wrong endpoint certified from a matching-looking word')
    try:
        Arc(1,1,(2,),'hostile').validate()
    except ValueError:
        pass
    else:
        raise ValueError('root self-arc accepted')
    print('HOSTILES:unproved supplied source; disconnected edge; wrong endpoint; root arc; five malformed source types.')
    print('PASS: conditional completeness for finite home routes; no assertion that every supplied source terminates.')
    if export_path is not None:
        export_certificates(runs[True],export_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export-certificates',metavar='PATH',
                        help='write and round-trip verify exactly the frozen239 guarded-run first-hit routes')
    main(parser.parse_args().export_certificates)
