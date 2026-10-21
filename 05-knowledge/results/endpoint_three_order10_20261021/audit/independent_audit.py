"""Independent proof audit. Standard library only; no producer code imports.

Uses sets of variable labels for GF(2) polynomials, independent row reduction,
and directed-arc sets. Does not enumerate the 86,016 source assignments.
All checks use explicit exceptions and remain active with python -O.
"""
import argparse
import ast
import copy
import hashlib
import itertools
import json
from pathlib import Path
import re

V = tuple(range(10))
PAIRS = tuple(itertools.combinations(V, 2))
CUTS = ((1, 8), (2, 8), (2, 7))
ONE = frozenset((-1,))
ZERO = frozenset()
REPS = tuple(sorted({min((u, v), (9-v, 9-u)) for u, v in itertools.combinations(range(1, 9), 2) if v > u+1}))

def demand(ok, why):
    if not ok:
        raise RuntimeError(why)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def value(poly, mask):
    return sum(1 if k == -1 else (mask >> k) & 1 for k in poly) % 2

def reduction(equations):
    """Return coefficient rank, inconsistency and an independent reduced basis."""
    basis = {}
    inconsistent = False
    for polynomial in equations:
        row = polynomial
        while row - ONE:
            pivot = min(row - ONE)
            if pivot in basis:
                row = row ^ basis[pivot]
            else:
                basis[pivot] = row
                break
        if row == ONE:
            inconsistent = True
    return len(basis), inconsistent, basis

def symbolic_source(incoming):
    matrix = {}
    for u, v in PAIRS:
        if u == 0:
            expression = ONE if v < 9 and v not in incoming else ZERO
        elif v == 9:
            expression = ZERO if 9-u in incoming else ONE
        elif v == u+1:
            expression = ONE
        else:
            expression = frozenset((REPS.index(min((u, v), (9-v, 9-u))),))
        matrix[u, v], matrix[v, u] = expression, expression ^ ONE
    return matrix

def verify_parameter_space(incoming, matrix):
    """Start from all 45 independent tournament arcs, not the 12-bit claim."""
    raw = {}
    for index, (u, v) in enumerate(PAIRS):
        raw[u, v] = frozenset((index,))
        raw[v, u] = frozenset((index, -1))
    constraints = [raw[u, v] ^ raw[9-v, 9-u] for u, v in PAIRS]
    constraints += [raw[u, u+1] ^ ONE for u in range(9)]
    constraints += [raw[9, 0] ^ ONE]
    constraints += [raw[0, v] ^ (ONE if v not in incoming else ZERO) for v in range(1, 9)]
    rank, bad, _ = reduction(constraints)
    demand(not bad and rank == 33, 'full 45-arc source constraints must have dimension 12')
    for polynomial in constraints:
        substituted = ZERO
        for term in polynomial:
            substituted ^= ONE if term == -1 else matrix[PAIRS[term]]
        demand(substituted == ZERO, '12-bit parametrization violates a source constraint')
    demand(tuple(matrix[p] for p in REPS) == tuple(frozenset((k,)) for k in range(12)), 'free coordinates not independent')
    demand(all(matrix[u, v] == matrix[9-v, 9-u] for u, v in PAIRS), 'source reflection')
    demand(sum(matrix[0, v] == ONE for v in range(1, 10)) == 6, 'first degree')
    demand(sum(matrix[9, v] == ONE for v in range(9)) == 3, 'last degree')

def path_edges(word):
    demand(len(word) == 10 and set(word) == set(V), 'word is not a permutation')
    return {frozenset((u, v)) for u, v in zip(word, word[1:])}

def move(matrix, word):
    kept = path_edges(word)
    return {arc: polynomial ^ (ZERO if frozenset(arc) in kept else ONE) for arc, polynomial in matrix.items()}

def realize(matrix, mask):
    return {arc for arc, polynomial in matrix.items() if value(polynomial, mask)}

def arc_rows(arcs):
    return [sum(1 << v for u, v in arcs if u == row) for row in V]

def from_rows(rows):
    demand(len(rows) == 10 and all(type(row) is int and 0 <= row < 1024 for row in rows), 'row shape')
    arcs = {(u, v) for u in V for v in V if (rows[u] >> v) & 1}
    demand(all((u, u) not in arcs for u in V), 'diagonal arc')
    demand(all(((u, v) in arcs) != ((v, u) in arcs) for u, v in PAIRS), 'not a tournament')
    return arcs

def colour_failures(arcs, word):
    reversal = dict(zip(word, reversed(word)))
    return [(u, v, reversal[v], reversal[u]) for u, v in PAIRS
            if ((u, v) in arcs) != ((reversal[v], reversal[u]) in arcs)]

def literal_move(arcs, word):
    kept = set(zip(word, word[1:]))
    path_edges(word)
    demand(kept <= arcs, 'literal word is not directed in its actual source')
    return kept | {(v, u) for u, v in arcs - kept}

def replay(rows, words):
    now = from_rows(rows)
    states, witnesses, score_sums = [], [], []
    for word in words:
        states.append(now)
        demand(colour_failures(now, word), 'literal word is blue')
        witness = colour_failures(now, word)[0]
        witnesses.append(list(witness) + [int(tuple(witness[:2]) in now), int(tuple(witness[2:]) in now)])
        score_sums.append(sum(u in (word[0], word[-1]) for u, v in now))
        now = literal_move(now, word)
    demand(now == literal_move(from_rows(rows), V), 'literal labelled target mismatch')
    demand(len(states) == 3, 'not three literal steps')
    return states, witnesses, score_sums, now

def prove_xor(equations, indices):
    demand(type(indices) is list and all(type(i) is int and 0 <= i < len(equations) for i in indices), 'proof index range')
    demand(len(indices) == len(set(indices)), 'repeated proof index')
    result = ZERO
    for i in indices:
        result ^= equations[i]
    return result

def certificate_audit(data):
    demand(data['order'] == 10 and data['dimension'] == 12 and data['cut_pairs'] == [list(x) for x in CUTS], 'certificate header')
    demand(data['variables'] == [sorted([list(p) for p in {(u, v), (9-v, 9-u)}]) for u, v in REPS], 'certificate variable order')
    expected = {(ports, bits) for ports in itertools.combinations(range(2, 9), 2) for bits in itertools.product((0, 1), repeat=3)}
    seen, residuals, system_ranks, branches = set(), [], {}, {}
    joins = contradictory = units = 0
    for record in data['cells']:
        ports, bits = tuple(record['incoming']), tuple(record['cell'])
        key = ports, bits
        demand(key in expected and key not in seen, 'extra, duplicate or invalid cell')
        seen.add(key)
        source = symbolic_source(ports)
        verify_parameter_space(ports, source)
        guards = [source[p] ^ (ONE if bit else ZERO) for p, bit in zip(CUTS, bits)]
        rank, bad, _ = reduction(guards)
        demand(rank == 3 and not bad, 'cell must have dimension nine')
        fixed = {next(iter(row - ONE)): (-1 in row) for row in guards}
        demand(len(fixed) == 3 and all(len(row - ONE) == 1 for row in guards), 'cuts must be independent source coordinates')
        demand(len(record['words']) == len(record['proofs']) == 3, 'wrong step count')
        now = source
        edges = set()
        branches[record['branch']] = branches.get(record['branch'], 0) + 1
        for step, (word, proof) in enumerate(zip(record['words'], record['proofs']), 1):
            edges ^= path_edges(word)
            for arc in zip(word, word[1:]):
                expr = now[arc] ^ ONE
                for variable, bit in fixed.items():
                    if variable in expr:
                        expr ^= frozenset((variable,)) ^ (ONE if bit else ZERO)
                demand(not expr, 'join not identically directed in cell')
                joins += 1
            reversal = dict(zip(word, reversed(word)))
            equations = guards + [now[u, v] ^ now[reversal[v], reversal[u]] for u, v in PAIRS]
            rank, bad, _ = reduction(equations)
            system_ranks[str((rank, bad))] = system_ranks.get(str((rank, bad)), 0) + 1
            if set(proof) == {'contradiction'}:
                demand(prove_xor(equations, proof['contradiction']) == ONE, 'invalid contradiction certificate')
                demand(bad, 'independent elimination fails to confirm contradiction')
                contradictory += 1
            else:
                demand(set(proof) == {'assignment', 'units'}, 'invalid certificate kind')
                assignment = proof['assignment']
                demand(type(assignment) is int and 0 <= assignment < 4096, 'assignment range')
                demand(not bad and rank == 12, 'residual must be exactly one point')
                demand(not any(value(row, assignment) for row in equations), 'residual does not solve all equations')
                demand(len(proof['units']) == 12, 'missing unit proof')
                for k, indices in enumerate(proof['units']):
                    demand(prove_xor(equations, indices) == frozenset((k,)) ^ (ONE if (assignment >> k) & 1 else ZERO), 'invalid unit proof')
                    units += 1
                residuals.append({'incoming': list(ports), 'cell': list(bits), 'branch': record['branch'], 'step': step, 'mask': assignment})
            now = move(now, word)
        demand(edges == path_edges(V), 'literal labelled undirected XOR fails')
        demand(now == move(source, V), 'symbolic literal target fails')
    demand(seen == expected, 'incomplete source partition')
    demand(contradictory == 501 and len(residuals) == 3 and joins == 4536 and units == 36, 'summary counts')
    return residuals, branches, system_ranks

def check_witnesses(states, words, witnesses, require_arc=False):
    demand(len(witnesses) == 3, 'witness step count')
    for arcs, word, w in zip(states, words, witnesses):
        demand(len(w) == 4 and all(type(x) is int and x in V for x in w), 'witness shape')
        u, v, a, b = w
        reversal = dict(zip(word, reversed(word)))
        demand(u != v and (a, b) == (reversal[v], reversal[u]), 'not actual positional anti-image')
        demand(((u, v) in arcs) != ((a, b) in arcs), 'witness not an anti-reversal failure in actual source')
        if require_arc:
            demand((u, v) in arcs, 'producer witness must begin with actual source arc')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--packet', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    packet = args.packet
    pinned = {}
    for line in (packet/'MANIFEST.sha256').read_text().splitlines():
        expected, name = line.split()
        demand(name not in pinned and digest(packet/name) == expected, 'packet hash mismatch: ' + name)
        pinned[name] = expected
    demand(len(pinned) == 10, 'packet manifest size')
    pins = json.loads((packet/'source_pins.json').read_text())
    if args.repo:
        for name, expected in pins['sources'].items():
            demand(digest(args.repo/name) == expected, 'inherited source hash mismatch: ' + name)
    supplemental_path = Path(__file__).with_name('supplemental_dependencies.json')
    supplemental = json.loads(supplemental_path.read_text())
    if args.repo:
        for name, expected in supplemental['sources'].items():
            demand(digest(args.repo/name) == expected, 'supplemental dependency hash mismatch: ' + name)
        def ledger_ast(path):
            return ast.dump(next(node for node in ast.parse(path.read_text()).body
                                 if isinstance(node, ast.FunctionDef) and node.name == 'ledger'))
        demand(ledger_ast(packet/'check.py') == ledger_ast(args.repo/'05-knowledge/results/endpoint_three_repairs_20261020h/root_check.py'), 'inherited selector formulas changed')
    text = (packet/'README.md').read_text()
    printed_reps = tuple(ast.literal_eval('[' + re.search(r'following order:\s*\n\s*((?:\([0-9, ]+\),?\s*)+)', text).group(1).strip().rstrip(',') + ']'))
    demand(printed_reps == REPS, 'printed free-orbit order is wrong')
    data = json.loads((packet/'certificates.json').read_text())
    residuals, branches, system_ranks = certificate_audit(data)
    table = []
    for ports, bits, branch, step, mask in re.findall(r'^\| \{([0-9,]+)\} \| \(([0-9,]+)\) \| ([a-z_]+) \| ([123]) \| ([0-9]+) \|$', text, re.M):
        table.append({'incoming': list(ast.literal_eval('('+ports+',)')), 'cell': list(ast.literal_eval('('+bits+')')), 'branch': branch, 'step': int(step), 'mask': int(mask)})
    demand(table == residuals, 'printed residual table disagrees with exact systems')
    exceptions = json.loads((packet/'exceptions.json').read_text())
    demand(len(exceptions) == 3, 'exception count')
    demand({(tuple(r['incoming']), r['mask']) for r in residuals} == {(tuple(f['incoming']), f['mask']) for f in exceptions}, 'exception coverage')
    printed = re.findall(r'### I=\{([0-9,]+)\}, J=\{([0-9,]+)\}, mask ([0-9]+)\s+T=(\([^\n]+\))\.\s+Q1=(\([^\n]+\))\s+Q2=(\([^\n]+\))\s+Q3=(\([^\n]+\))', text)
    demand(len(printed) == 3, 'missing literal printed examples')
    receipts = [json.loads((packet/name).read_text()) for name in ('receipt.normal.json', 'independent.normal.json')]
    replay_receipts = []
    for index, (fix, display) in enumerate(zip(exceptions, printed)):
        ports, out, mask, rows, q1, q2, q3 = display
        demand(list(ast.literal_eval('('+ports+',)')) == fix['incoming'], 'printed incoming ports')
        demand(list(ast.literal_eval('('+out+',)')) == fix['outgoing'], 'printed outgoing ports')
        demand(int(mask) == fix['mask'] and list(ast.literal_eval(rows)) == fix['rows'], 'printed source mask or rows')
        demand([list(ast.literal_eval(q)) for q in (q1, q2, q3)] == fix['words'], 'printed words differ from data')
        source = symbolic_source(fix['incoming'])
        demand(arc_rows(realize(source, fix['mask'])) == fix['rows'], 'source does not match coordinate point')
        demand(sorted(9-u for u in fix['incoming']) == fix['outgoing'], 'incorrect reflected outgoing ports')
        states, witnesses, scores, target = replay(fix['rows'], fix['words'])
        demand(arc_rows(target) == fix['partner_rows'], 'stored partner rows wrong')
        check_witnesses(states, fix['words'], receipts[0]['literal_exception_replays'][index]['anti_failures'], True)
        check_witnesses(states, fix['words'], receipts[1]['numeric_replays'][index]['anti_failures'])
        replay_receipts.append({'incoming': fix['incoming'], 'mask': fix['mask'], 'actual_source_anti_failures_with_values': witnesses, 'endpoint_score_sums': scores, 'literal_target_rows': arc_rows(target)})
    demand(replay_receipts[0]['endpoint_score_sums'] == [9, 9, 9], 'printed endpoint-score warning')
    hostile_match = re.search(r'sharp cyclic-cut hostile is also retained:\s*(\([^\n]+\))\.', text)
    demand(hostile_match is not None, 'missing printed hostile')
    hostile_rows = list(ast.literal_eval(hostile_match.group(1)))
    hostile = from_rows(hostile_rows)
    demand(not colour_failures(hostile, V) and (9, 0) in hostile, 'hostile source must be closed blue')
    demand(sum(u == 0 for u, v in hostile) == 6 and sum(u == 9 for u, v in hostile) == 3, 'hostile endpoint scores')
    cut = tuple(range(1, 10)) + (0,)
    literal_move(hostile, cut)
    demand(not colour_failures(hostile, cut), 'hostile cyclic cut must be blue')
    states, witnesses, _, target = replay(hostile_rows, exceptions[2]['words'])
    check_witnesses(states, exceptions[2]['words'], receipts[0]['cyclic_cut_hostile_repaired'], True)
    # Negative controls exercise independent proof checking, rather than trusting assertions.
    controls = {}
    for name in ('omit_cell', 'duplicate_cell', 'corrupt_contradiction', 'corrupt_singleton', 'invalid_word'):
        changed = copy.deepcopy(data)
        if name == 'omit_cell':
            changed['cells'].pop()
        elif name == 'duplicate_cell':
            changed['cells'].append(changed['cells'][0])
        elif name == 'corrupt_contradiction':
            p = next(p for c in changed['cells'] for p in c['proofs'] if 'contradiction' in p)
            p['contradiction'] = []
        elif name == 'corrupt_singleton':
            p = next(p for c in changed['cells'] for p in c['proofs'] if 'assignment' in p)
            p['assignment'] ^= 1
        else:
            changed['cells'][0]['words'][0][0] = changed['cells'][0]['words'][0][1]
        try:
            certificate_audit(changed)
        except RuntimeError as error:
            controls[name] = str(error)
        else:
            raise RuntimeError('negative control was accepted: ' + name)
    demand((packet/'receipt.normal.json').read_bytes() == (packet/'receipt.optimized.json').read_bytes(), 'producer optimized receipts differ')
    demand((packet/'independent.normal.json').read_bytes() == (packet/'independent.optimized.json').read_bytes(), 'producer verifier optimized receipts differ')
    result = {
        'status': 'PASS: independent mathematical and literal-text audit',
        'scope': 'FINITE-EXACT n=10; inherited structural n>=11 theorem is not reproved here',
        'source_universe': {'raw_arc_variables': 45, 'constraint_rank': 33, 'free_coordinates': 12, 'port_pairs': 21, 'cells': 168, 'cell_dimension': 9, 'marked_presentations': 86016, 'assignment_census_performed': False, 'filters': []},
        'affine_verification': {'systems': 504, 'inconsistent': 501, 'singleton': 3, 'unit_certificates': 36, 'joins': 4536, 'literal_XOR_targets': 168, 'rank_distribution': system_ranks, 'branch_cells': branches},
        'residuals': residuals, 'printed_repair_checks': replay_receipts,
        'hostile': {'rows': hostile_rows, 'cyclic_cut': list(cut), 'cyclic_cut_is_blue': True, 'repair_actual_source_anti_failures_with_values': witnesses, 'literal_target_rows': arc_rows(target)},
        'stored_actual_source_witnesses_verified': 21,
        'negative_controls': controls,
        'source_packet_manifest_sha256': digest(packet/'MANIFEST.sha256'),
        'source_packet_hashes': pinned,
        'inherited_repository_commit': pins['repository_commit'],
        'inherited_pins_verified': len(pins['sources']) if args.repo else 0,
        'supplemental_dependencies_sha256': digest(supplemental_path),
        'supplemental_dependencies_verified': len(supplemental['sources']) if args.repo else 0,
        'inherited_selector_ast_equality': True if args.repo else 'not checked without --repo',
        'checker_sha256': digest(Path(__file__))
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
