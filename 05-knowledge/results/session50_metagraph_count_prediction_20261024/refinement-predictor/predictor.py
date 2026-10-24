#!/usr/bin/env python3
"""Blind order-eight count: attachment-subset orbits and directed refinement.

Reads only the frozen parent atlas and this source. The structural theorem
sets the strong order-eight refinement error to zero. All other errors are
computed from input parent actions. No result receipt or generation module is
read or imported. The executable deliberately rejects target order nine.
"""
import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

ATLAS_SHA256 = 'df16881db1060e0423ed73658bdcb5923266f839681b74ac5b1bddd5a232f3ad'
OPERATIONS = Counter()
CHECKS = Counter()


def require(ok, message):
    CHECKS[message] += 1
    if not ok:
        raise RuntimeError(message)


def enforce_source_boundary(source):
    """A narrow, inspected call/import surface, independent of Python asserts."""
    tree = ast.parse(source)
    allowed_imports = {'argparse', 'ast', 'collections', 'hashlib', 'json', 'pathlib', 'sys'}
    forbidden_calls = {'canonical', 'canonicalize', 'canonical_label', 'deletion',
                       'delete_vertex', 'rooted_key', 'lookup', 'reroot',
                       'automorphisms', 'permutations', 'permute', 'eval', 'exec',
                       'compile', '__import__', 'getattr', 'setattr', 'open',
                       'load_module', 'exec_module', 'run', 'Popen', 'system'}
    # These are output/control constants, not legitimate blind inputs.
    forbidden_numeric_inputs = {6880, 55040, 55068, 54256, 28, 14}
    imports = []
    calls = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for item in node.names:
                require(item.name in allowed_imports, 'only audited standard-library imports')
                imports.append(item.name)
        elif isinstance(node, ast.ImportFrom):
            require(node.module in allowed_imports and node.level == 0,
                    'only audited standard-library imports')
            imports.append(node.module)
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                name = node.func.id
            elif isinstance(node.func, ast.Attribute):
                name = node.func.attr
            else:
                raise RuntimeError('indirect calls are forbidden')
            require(name not in forbidden_calls, 'no forbidden isomorphism or dynamic calls')
            calls.add(name)
    # The six integer literals above are allowed only inside this boundary's
    # diagnostic set, never in prediction, input handling, or main.
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == 'enforce_source_boundary':
            continue
        for descendant in ast.walk(node):
            if isinstance(descendant, ast.Constant) and type(descendant.value) is int:
                require(descendant.value not in forbidden_numeric_inputs,
                        'no answer or trace constants in prediction')
    return {'imports': sorted(set(imports)), 'call_names': sorted(calls),
            'forbidden_isomorphism_calls_in_source': 0,
            'target_or_trace_constants_outside_guard': 0}


def restrict_io(source_path, atlas_path, output_path):
    readable = {str(source_path), str(atlas_path)}
    writable = str(output_path)

    def audit_hook(event, arguments):
        if event == 'open':
            target, mode, flags = arguments
            require(isinstance(target, (str, bytes)), 'no descriptor-based extra input')
            path = str(Path(target).resolve())
            writing = mode is not None and any(c in mode for c in 'wax+')
            require(path == writable if writing else path in readable,
                    'blind prediction file access boundary')
        require(not event.startswith(('subprocess.', 'socket.', 'ctypes.')),
                'no external execution or network during prediction')
        require(event not in {'os.system', 'os.exec', 'os.spawn'},
                'no external execution or network during prediction')

    sys.addaudithook(audit_hook)


def stable_colors(rows):
    n = len(rows)
    colors = tuple(row.bit_count() for row in rows)
    while True:
        OPERATIONS['refinement_rounds'] += 1
        labels = sorted(set(colors))
        signatures = tuple((colors[v], tuple(sum((rows[v] >> w) & 1
                           for w in range(n) if colors[w] == color)
                           for color in labels)) for v in range(n))
        names = {signature: index for index, signature in enumerate(sorted(set(signatures)))}
        refined = tuple(names[signature] for signature in signatures)
        if len(set(refined)) == len(set(colors)):
            return refined
        colors = refined


def check_parent(rows, actions):
    n = len(rows)
    require(all(type(row) is int and 0 <= row < (1 << n) for row in rows),
            'parent bitsets have valid width')
    require(all(not ((rows[v] >> v) & 1) for v in range(n)), 'parent has no loops')
    require(all(((rows[v] >> w) & 1) + ((rows[w] >> v) & 1) == 1
                for v in range(n) for w in range(v)), 'parent is a tournament')
    group = set(actions)
    require(len(group) == len(actions) and tuple(range(n)) in group,
            'parent action set distinct and contains identity')
    require(all(tuple(sorted(g)) == tuple(range(n)) for g in group),
            'parent actions are permutations')
    for g in group:
        OPERATIONS['parent_action_arc_checks'] += n * n
        require(all(((rows[v] >> w) & 1) == ((rows[g[v]] >> g[w]) & 1)
                    for v in range(n) for w in range(n)),
                'supplied parent actions preserve arcs')
        for h in group:
            require(tuple(g[h[v]] for v in range(n)) in group,
                    'supplied parent actions form a group')
    return group


def vertex_partition(actions, n):
    remaining = set(range(n))
    blocks = []
    while remaining:
        v = min(remaining)
        block = {g[v] for g in actions}
        require(block <= remaining, 'parent action orbits partition vertices')
        blocks.append(tuple(sorted(block)))
        remaining -= block
    return blocks


def is_strong(rows):
    n = len(rows)
    if not n:
        return False
    for start in range(n):
        reached = 1 << start
        while True:
            following = reached
            for v in range(n):
                if (reached >> v) & 1:
                    following |= rows[v]
            if following == reached:
                break
            reached = following
        if reached != (1 << n) - 1:
            return False
    return True


def subset_image(mask, action):
    OPERATIONS['parent_subset_action_evaluations'] += 1
    return sum(1 << action[v] for v in range(len(action)) if (mask >> v) & 1)


def attach_vertex(rows, incoming):
    n = len(rows)
    child = list(rows) + [0]
    for v in range(n):
        if (incoming >> v) & 1:
            child[v] |= 1 << n
        else:
            child[n] |= 1 << v
    OPERATIONS['child_extensions'] += 1
    return tuple(child)


def predict(payload):
    require(payload['schema'] == 'tournament-parent-atlas-actions-v1', 'input atlas schema')
    require(payload['maximum_parent_order'] == 7, 'input parent order is seven')
    levels = payload['levels']
    require([level['order'] for level in levels] == list(range(8)),
            'input has exactly orders zero through seven')
    counts = [len(level['parents']) for level in levels]
    require(counts[0] == counts[1] == 1, 'empty and singleton base classes')
    strong_errors = [0] * 9
    total_parent_errors = [0] * 8
    positive_parents = []
    final_parents = []
    for level in levels:
        n = level['order']
        parent_rows = set()
        for record in level['parents']:
            require(set(record) == {'rows', 'actions'}, 'parent record contains only rows and actions')
            rows = tuple(record['rows'])
            actions = tuple(tuple(g) for g in record['actions'])
            require(len(rows) == n and rows not in parent_rows, 'distinct labelled parent records')
            parent_rows.add(rows)
            group = check_parent(rows, actions)
            if n:
                OPERATIONS['parent_refinement_calls'] += 1
                colors = stable_colors(rows)
                orbits = vertex_partition(group, n)
                error = 0
                for orbit in orbits:
                    require(len({colors[v] for v in orbit}) == 1,
                            'parent orbit contained in stable cell')
                    error += colors.count(colors[orbit[0]]) - len(orbit)
                total_parent_errors[n] += error
                strong = is_strong(rows)
                if strong:
                    strong_errors[n] += error
                if error:
                    positive_parents.append({'order': n, 'rows': list(rows),
                        'action_order': len(group), 'vertex_orbit_sizes': sorted(map(len, orbits)),
                        'stable_cell_sizes': sorted(Counter(colors).values()),
                        'strong': strong, 'error': error})
            if n == 7:
                final_parents.append((rows, group))
    require(all(error == 0 for error in strong_errors[:7]),
            'computed lower strong errors agree with structural vanishing')
    # THEOREM.md of the separate holonomy packet proves this zero without
    # enumerating, comparing, or canonically labeling any order-eight child.
    strong_errors[8] = 0
    for n in range(1, 8):
        inherited = sum(strong_errors[k] * sum(counts[i] * counts[n-k-i]
                        for i in range(n-k+1)) for k in range(1, n+1))
        require(inherited == total_parent_errors[n], 'SCC error identity on known parent levels')
    correction = sum(strong_errors[k] * sum(counts[i] * counts[8-k-i]
                     for i in range(8-k+1)) for k in range(1, 9))
    upper_trace = 0
    states = 0
    histogram = Counter()
    stream = hashlib.sha256()
    for rows, actions in final_parents:
        seen = set()
        for mask in range(1 << len(rows)):
            if mask in seen:
                continue
            orbit = {subset_image(mask, g) for g in actions}
            require(mask == min(orbit) and not orbit & seen, 'attachment orbit represented exactly once')
            seen.update(orbit)
            child = attach_vertex(rows, mask)
            OPERATIONS['child_refinement_calls'] += 1
            colors = stable_colors(child)
            cell_size = colors.count(colors[-1])
            histogram[cell_size] += 1
            upper_trace += cell_size
            states += 1
            stream.update((json.dumps([rows, mask, cell_size], separators=(',', ':')) + '\n').encode())
        require(seen == set(range(1 << len(rows))), 'attachment orbits cover all subsets')
    numerator = upper_trace - correction
    require(numerator > 0 and numerator % 8 == 0, 'corrected trace is a positive multiple of eight')
    require(OPERATIONS['child_extensions'] == OPERATIONS['child_refinement_calls'] == states,
            'one refinement and extension per pointed state')
    return {'status': 'PASS', 'target_order': 8, 'input_parent_class_counts': counts,
            'pointed_extension_states': states, 'stable_root_cell_histogram': sorted(histogram.items()),
            'upper_refinement_trace': upper_trace, 'computed_parent_total_errors': total_parent_errors,
            'computed_strong_errors_through_seven': strong_errors[:8],
            'theorem_strong_error_order_eight': strong_errors[8], 'SCC_total_error_order_eight': correction,
            'corrected_unpointing_trace': numerator, 'predicted_A000568_order_eight': numerator // 8,
            'positive_error_input_parents': positive_parents, 'pointed_refinement_stream_sha256': stream.hexdigest(),
            'operation_counts': dict(sorted(OPERATIONS.items())),
            'child_canonicalization_calls': 0, 'deletion_or_reroot_isomorphism_calls': 0,
            'zero_call_basis': 'source boundary, standard-library-only imports, and restricted runtime input',
            'scope': 'finite order eight, conditional on complete parent atlas and the proved strong-eight vanishing theorem'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--parents', default=str(Path(__file__).resolve().parent / 'parents-through-seven.json'))
    parser.add_argument('--receipt', default='receipt.json')
    parser.add_argument('--target-order', type=int, choices=[8], default=8)
    args = parser.parse_args()
    source = Path(__file__).resolve()
    atlas = Path(args.parents).resolve()
    output = Path(args.receipt).resolve()
    require(output not in {source, atlas}, 'output cannot overwrite source or input')
    restrict_io(source, atlas, output)
    boundary = enforce_source_boundary(source.read_text())
    raw = atlas.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == ATLAS_SHA256, 'frozen input atlas hash')
    result = predict(json.loads(raw))
    result['input_atlas_sha256'] = ATLAS_SHA256
    result['source_boundary'] = boundary
    # Check counts are not part of the deterministic receipt: the interpreter
    # may emit different audit event counts under -O.
    output.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
