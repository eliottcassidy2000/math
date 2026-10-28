"""Replay the two historical bugs without executing their benchmark blocks."""
import ast
import importlib.util
import json
import subprocess
from fractions import Fraction
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = 'b8672af00bcc6d440a9a00a6915325faf39fcdff'
ORACLE = HERE.parent/'session86_relations_attachment_carriers_20261028/recreate_a000568.py'
spec = importlib.util.spec_from_file_location('cyclic_pair_oracle', ORACLE)
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)

def functions_only(source, name):
    tree = ast.parse(source, filename=name)
    tree.body = [node for node in tree.body if isinstance(
        node, (ast.Import, ast.ImportFrom, ast.FunctionDef))]
    scope = {}
    exec(compile(tree, name, 'exec'), scope)
    return scope

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def main():
    reference = []
    for n in range(31):
        value = sum((Fraction(1 << oracle.cyclic_burnside_exponent(part),
                              oracle.centralizer(part))
                     for part in oracle.partitions(n)), Fraction())
        need(value.denominator == 1, 'independent pair-action integrality')
        reference.append(int(value))
    sources = []
    for path, function, old_expected in (
        ('04-computation/a000568_speedup.py', 'a000568_davis', [1,1,1,2,4,10,56]),
        ('04-computation/a000568_v3_s90dm.py', 'T_v3', [1,1,2,4,11,34,156]),
    ):
        original = subprocess.check_output(
            ['git', 'show', BASE+':'+path], cwd=ROOT).decode('utf-8')
        old = functions_only(original, path+'@'+BASE)[function]
        fixed = functions_only((ROOT/path).read_text(encoding='utf-8'), path)[function]
        before = [old(n) for n in range(7)]
        after = [fixed(n) for n in range(31)]
        need(before == old_expected, 'historical witness preserved')
        need(after == reference, 'repaired function vs independent pair action')
        sources.append({'path':path, 'old_counts_0_6':before,
                        'repaired_counts_0_30':after})
    # Literal pair-orbit tracing is a third path at the cheap small orders.
    for n in range(9):
        value = sum((Fraction(oracle.direct_fixed_orientations(part),
                              oracle.centralizer(part))
                     for part in oracle.partitions(n, odd=False)), Fraction())
        need(value == reference[n], 'unfiltered ordered-pair consistency audit')
    ratios = {str(n):str(Fraction(n*reference[n], (1 << (n-1))*reference[n-1]))
              for n in (10,20,30)}
    report = {'status':'PASS', 'base':BASE, 'sources':sources,
              'independent_cyclic_pair_action_verified_through':30,
              'unfiltered_direct_pair_action_verified_through':8,
              'exact_successive_asymptotic_ratios':ratios,
              'single_three_cycle_relative_at10':str(Fraction(2*comb(10,3),1 << 16))}
    (HERE/'legacy_repairs.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('PASS: legacy witnesses and repairs, independent pair actions through30.')

if __name__ == '__main__':
    main()
