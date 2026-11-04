"""Frozen checks and optional complete replay; expectations are never rewritten."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ENGINES = (
    ('core7', 'cores/run.py', ['--order', '7'], 'cores/order7.json', 'cores/order7.optimized.json'),
    ('core8', 'cores/run.py', ['--order', '8'], 'cores/order8.json', 'cores/order8.optimized.json'),
    ('core8_independent', 'cores/independent.py', [], 'cores/independent.json', 'cores/independent.optimized.json'),
    ('coefficients', 'recovery/check.py', [], 'recovery/normal.json', 'recovery/optimized.json'),
    ('incoming', 'recovery/incoming_check.py', [], 'recovery/incoming.json', 'recovery/incoming.optimized.json'),
    ('geometry', 'geometry/verify.py', [], 'geometry/receipt.json', 'geometry/receipt.optimized.json'),
    ('geometry_review', 'recovery/geometry_review.py', [], 'recovery/geometry_review.json', 'recovery/geometry_review.optimized.json'),
    ('scout', 'algebra/run.py', [], 'algebra/normal.json', 'algebra/optimized.json'),
    ('modules', 'algebra/module_check.py', [], 'algebra/modules.json', 'algebra/modules.optimized.json'),
    ('response', 'response/run.py', [], 'response/normal.json', 'response/optimized.json'),
    ('algebra_review', 'recovery/algebra_review.py', [], 'recovery/algebra_review.json', 'recovery/algebra_review.optimized.json'),
)


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads((HERE / path).read_text(encoding='utf8'))


def pins():
    files = load('MANIFEST.json')['files']
    actual = {p.relative_to(HERE).as_posix() for p in HERE.rglob('*')
              if p.is_file() and p.name != 'MANIFEST.json'
              and '__pycache__' not in p.parts and p.suffix != '.pyc'}
    need(set(files) == actual, 'artifact inventory changed')
    for name, digest in files.items():
        need(sha(HERE / name) == digest, 'artifact changed: ' + name)
    sources = load('SOURCE_PINS.json')['files']
    for name, digest in sources.items():
        need(sha(ROOT / name) == digest, 'inherited source changed: ' + name)
    return len(files), len(sources)


def main():
    need(sys.argv[1:] in ([], ['--replay']), 'usage: verify.py [--replay]')
    before = pins()
    expected = {}
    for key, script, args, normal, optimized in ENGINES:
        need(key not in expected, 'duplicate engine key')
        expected[key] = load(normal)
        need(expected[key] == load(optimized), 'normal/optimized mismatch: ' + key)
    core = expected['core8']
    need((core['all_graphs'], core['min_degree3_graphs'], core['metric_filter_graphs'], core['classes'])
         == (2**28, 53470067, 724080, 38), 'complete core scope')
    independent = expected['core8_independent']
    need(independent['universe'] == 2**28 and independent['eligible_graphs'] == 724080
         and independent['orbit_union_equals_complete_survivor_set'], 'independent complete core scope')
    coeff = expected['coefficients']
    need(coeff['order8_certificates'] == 30 and coeff['order8_complete_orientations'] == 548864,
         'independent coefficient universe')
    geo = expected['geometry']
    need((geo['moser_exterior_exclusions'], geo['circle_fan_exclusions'], geo['all_actual_pairs_checked'],
          geo['actual_unit_pairs'], geo['control_chromatic_number']) == (6, 2, 105, 25, 4), 'geometry scope')
    scout = expected['scout']
    need(scout['total_candidates'] == 311872, 'specified scout scope')
    need(sum(r['endpoint_parity_balanced'] for r in scout['families'].values()) == 146
         and sum(r['joint_survivors'] for r in scout['families'].values()) == 0,
         'specified scout necessary conditions')
    modules = expected['modules']
    need(modules['six_module_rooted_expansions'] == 336 and modules['six_module_joint_survivors'] == 0,
         'forced module scope')
    need(len(modules['prime_controls']) == 2 and all(r['tested_subsets'] == 4082
         and r['proper_nontrivial_modules'] == 0 for r in modules['prime_controls']), 'prime controls')
    response = expected['response']
    need(response['distinguishing_battery_entries'] == 6240, 'actual response separation')
    need(len(response['controls']) == 2, 'specified full tensor universe')
    need(sum(r['queries'] for r in response['controls']) == 12500, 'direct query inventory')
    for r in response['controls']:
        for name, power in [('S', 1), ('E', 1), ('M', 2), ('C', 2), ('U', 3), ('V', 3), ('D', 4)]:
            need(len(r[name]) == 12**power, 'full named tensor: ' + name)
        need(len(r['battery']) == 6241, 'complete interpolation battery')
        need(sum(r['S']) == sum(r['E']) == sum(r['M']), 'endpoint sums')
        need(sum(r['C']) == sum(r['U']) == sum(r['V']), 'paired-path sums')
        need(sum(r['D']) % 6 == 0, 'ordered three-path factor')
    witness = response['first_actual_query_witness']
    need((witness['H_T_plus_ports'], witness['H_converse_plus_same_ports']) == (627597, 1139039),
         'actual distinguishing query')
    replayed = 0
    if sys.argv[1:] == ['--replay']:
        for key, script, args, _, _ in ENGINES:
            for flags in ([], ['-O']):
                result = subprocess.run([sys.executable, '-B'] + flags + [str(HERE / script)] + args,
                                        cwd=(HERE / script).parent, check=True,
                                        stdout=subprocess.PIPE, text=True, encoding='utf8')
                need(json.loads(result.stdout) == expected[key], 'frozen replay differs: ' + key)
                replayed += 1
                print('PASS ' + key + (' -O' if flags else ''), flush=True)
    need(pins() == before, 'files changed during verification')
    print(json.dumps({'status': 'PASS', 'artifact_files': before[0], 'inherited_source_files': before[1],
                      'matching_receipt_pairs': len(ENGINES), 'replayed_engines': replayed,
                      'complete_eight_graph_universe': 2**28, 'twelve_scout_records': 311872,
                      'universal_two_port': 'OPEN', 'six_chromatic_witness': False,
                      'expectations_rewritten': False}, sort_keys=True))


if __name__ == '__main__':
    main()
