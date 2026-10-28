#!/usr/bin/env python3
"""Audit live S325/S327 codecs, without executing historical explorations.

Defaults locate 04-computation/ beneath a repository ancestor of this file or
of the working directory. Use --repo-root or per-source overrides for staging.
Requires NumPy. Every check remains active under python -O.
"""
import argparse
import ast
import copy
import hashlib
import json
from pathlib import Path

import numpy as np

checks = 0


def require(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)


def assigned_name(node):
    if isinstance(node, ast.Assign) and len(node.targets) == 1:
        target = node.targets[0]
        if isinstance(target, ast.Name):
            return target.id
    return None


def loop_named(nodes, name):
    found = [node for node in nodes if isinstance(node, ast.For)
             and isinstance(node.target, ast.Name) and node.target.id == name]
    require(len(found) == 1, ('unique loop', name))
    return found[0]


def extract_s325(tree):
    found = [node for node in tree.body if isinstance(node, ast.FunctionDef)
             and node.name == 'build_square']
    require(len(found) == 1, 'build_square not unique')
    return copy.deepcopy(found[0])


def extract_s327(tree):
    # Select the unique census loop by its local rank_dist initialization.
    loops = [node for node in tree.body if isinstance(node, ast.For)
             and isinstance(node.target, ast.Name) and node.target.id == 'n'
             and any(assigned_name(s) == 'rank_dist' for s in node.body)]
    require(len(loops) == 1, 'S327 matrix census loop not unique')
    outer = loops[0]
    cutoff = next(i for i, s in enumerate(outer.body)
                  if assigned_name(s) == 'rank_dist')
    setup = outer.body[:cutoff]
    # Fail closed rather than silently selecting unrelated exploratory code.
    require([assigned_name(s) for s in setup]
            == ['side', 'tiles_n', None, 'tiles_n1', None, 'mn', 'mn1'],
            'S327 setup shape changed')
    for position, name in ((2, 'tiles_n'), (4, 'tiles_n1')):
        statement = setup[position]
        require(isinstance(statement, ast.Expr)
                and isinstance(statement.value, ast.Call)
                and isinstance(statement.value.func, ast.Attribute)
                and isinstance(statement.value.func.value, ast.Name)
                and statement.value.func.value.id == name
                and statement.value.func.attr == 'sort'
                and not statement.value.args and not statement.value.keywords,
                ('S327 expected tile sort', name))
    inner = loop_named(loop_named(outer.body, 'nb').body, 'nb1')
    require(len(inner.body) == 5 and assigned_name(inner.body[0]) == 'M'
            and all(isinstance(s, ast.For) for s in inner.body[1:3])
            and all(isinstance(s, ast.AugAssign) for s in inner.body[3:]),
            'S327 inline matrix block shape changed')
    for statement, name in zip(inner.body[3:], ('rank_dist', 'det_dist')):
        require(isinstance(statement.target, ast.Subscript)
                and isinstance(statement.target.value, ast.Name)
                and statement.target.value.id == name,
                ('S327 excluded census operation', name))
    function = ast.parse('def build_square(nb, nb1, n):\n    pass\n').body[0]
    # Only the source's tile setup and matrix initialization/filling are used.
    # Neither Counter operation nor rank/determinant evaluation is executed.
    function.body = copy.deepcopy(setup + inner.body[:3]) + [ast.Return(ast.Name('M', ast.Load()))]
    return ast.fix_missing_locations(function)


def load_codec(path, extractor):
    data = path.read_bytes()
    tree = ast.parse(data, filename=str(path))
    function = extractor(tree)
    namespace = {'np': np}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['build_square'], hashlib.sha256(data).hexdigest()


def input_tiles(n):
    large = [(i, j) for i in range(n) for j in range(i + 2, n)]
    small = [(i, j) for i in range(n - 1) for j in range(i + 2, n - 1)]
    return large, small


def oracle(a, b, n):
    m = n - 2
    large, small = input_tiles(n)
    ai = {pair: k for k, pair in enumerate(large)}
    bi = {pair: k for k, pair in enumerate(small)}
    # Independent inverse-cell specification, not the producer's loop direction.
    return tuple((a >> ai[(c, r + 2)]) & 1 if r >= c
                 else (b >> bi[(r, c + 1)]) & 1
                 for r in range(m) for c in range(m))


def audit(encode, label):
    start_checks = checks
    rows = []
    for n in range(3, 7):
        m = n - 2
        large, small = input_tiles(n)
        seen = set()
        for a in range(1 << len(large)):
            for b in range(1 << len(small)):
                matrix = encode(a, b, n)
                require(matrix.shape == (m, m), (label, 'shape', n, a, b))
                word = tuple(int(v) for v in matrix.flat)
                require(word == oracle(a, b, n), (label, 'inverse-cell oracle', n, a, b))
                require(word not in seen, (label, 'injectivity', n, a, b))
                seen.add(word)
        require(len(seen) == 1 << (m * m), (label, 'surjectivity', n))
        rows.append({'tiling_sizes': [n, n - 1], 'matrix_side': m,
                     'all_input_pairs': len(seen)})
    # Larger-side coordinate checks supplement, not replace, the all-size proof.
    for m in range(1, 13):
        n = m + 2
        large, small = input_tiles(n)
        cases = [(0, 0)] + [(1 << k, 0) for k in range(len(large))]
        cases += [(0, 1 << k) for k in range(len(small))]
        cases += [((1 << len(large)) - 1, (1 << len(small)) - 1)]
        for a, b in cases:
            matrix = encode(a, b, n)
            require(matrix.shape == (m, m), (label, 'basis shape', n, a, b))
            require(tuple(int(v) for v in matrix.flat) == oracle(a, b, n),
                    (label, 'coordinate basis', n, a, b))
    require(not np.array_equal(encode(0, 0, 5), encode(0, 5, 5)),
            (label, 'legacy 000/101 collision survives'))
    return {'always_active_checks': checks - start_checks, 'universes': rows,
            'coordinate_basis_sides': '1..12 (zero, each basis bit, all ones)'}


def resolve_sources(args):
    names = ('tiling_square_s325.py', 'new_sequences_s327.py')
    explicit = (args.s325_source, args.s327_source)
    root = args.repo_root
    if root is None and not all(explicit):
        candidates = list(Path(__file__).resolve().parents)
        candidates += [Path.cwd(), *Path.cwd().parents]
        root = next((candidate for candidate in candidates
                     if all(override is not None or (candidate / '04-computation' / name).is_file()
                            for name, override in zip(names, explicit))), None)
        if root is None:
            raise RuntimeError('Cannot locate live 04-computation sources. '
                               'Pass --repo-root or both --s325-source and --s327-source.')
    paths = [override if override is not None else root / '04-computation' / name
             for name, override in zip(names, explicit)]
    for path in paths:
        if not path.is_file():
            raise RuntimeError('Source file does not exist: ' + str(path))
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path)
    parser.add_argument('--s325-source', type=Path)
    parser.add_argument('--s327-source', type=Path)
    args = parser.parse_args()
    paths = resolve_sources(args)
    results = []
    for path, extractor, name, description in zip(
            paths, (extract_s325, extract_s327),
            ('tiling_square_s325.py', 'new_sequences_s327.py'),
            ('named build_square', 'inline tile setup and matrix block')):
        encode, digest = load_codec(path, extractor)
        result = {'repository_path': '04-computation/' + name,
                  'tested_source_sha256': digest, 'AST_extraction': description}
        result.update(audit(encode, name))
        results.append(result)
    print(json.dumps({'status': 'PASS', 'always_active_checks': checks,
                      'codecs': results,
                      'historical_narrative_or_other_explorations_executed': False,
                      'limitations': ['Finite tests supplement the all-size coordinate proof.',
                                      'No historical rank/determinant outputs regenerated.',
                                      'No unrelated exploratory or OEIS claims certified.']}, indent=2))


if __name__ == '__main__':
    main()
