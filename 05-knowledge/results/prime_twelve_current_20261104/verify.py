"""Hash-pinned finite certificates; replay never overwrites retained evidence."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BANK = ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path, normalize=False):
    data = path.read_bytes()
    return sha256(data.replace(b'\r\n', b'\n') if normalize else data).hexdigest()


def read(name):
    return json.loads((HERE/name).read_text(encoding='utf8'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    manifest = read('manifest.json')
    inventory = {p.name for p in HERE.iterdir() if p.is_file() and p.name != 'manifest.json'}
    need(inventory == set(manifest['files']), 'packet inventory')
    for name, value in manifest['files'].items():
        need(digest(HERE/name) == value, 'packet hash: '+name)
    for name, value in manifest['dependencies_lf_sha256'].items():
        need(digest(ROOT/name, True) == value, 'dependency hash: '+name)
    for stem in ('current', 'switch', 'response', 'audit'):
        need((HERE/(stem+'.json')).read_bytes() ==
             (HERE/(stem+'.optimized.json')).read_bytes(), 'normal/-O: '+stem)
        pins = read(stem+'.json')['source_sha256']
        if isinstance(pins, str):
            pins = {stem+'.py': pins}
        for name, value in pins.items():
            need(digest(HERE/name) == value, 'producer source: '+name)
    current, orbit, response, audit = [read(s+'.json') for s in ('current', 'switch', 'response', 'audit')]
    need(digest(BANK) == current['bank_sha256'] == manifest['bank_sha256'], 'complete input bank')
    need(sum(v['classes'] for v in current['counts'].values()) == 198947, 'bank universe')
    need(len(current['records']) == 7 and not any(r['all_deletion_currents_zero'] for r in current['records']), 'current census')
    records = orbit['records']
    need([r['mask'] for r in records] == list(range(2048)), 'full orbit')
    need(all(r['S_parity'] ^ r['E_parity'] == 4095 for r in records), 'endpoint obstruction')
    need(sum(r['prime'] and r['strong'] for r in records) == 2024, 'strong prime presentations')
    need(dict(Counter(r['H'] for r in records)) == {95095: 24, 196815: 264, 236447: 1320, 256023: 440}, 'H strata')
    need([len({r['simple_paths_by_vertices'][k] for r in records}) for k in range(12)] == [1]*5+[4]*7, 'visited-set separation')
    need(digest(HERE/'switch.json') == audit['switch_receipt_sha256'], 'audit orbit input')
    need(audit['full_orbit_parity_checks'] == audit['full_orbit_closure_checks'] == 2048, 'independent orbit checks')
    need(audit['small_labelled_lattice_inputs'] == 1099 and audit['literal_one_port_current_checks'] == 1309, 'literal controls')
    need((HERE/'response_tables.json').read_bytes() == (HERE/'response_tables.optimized.json').read_bytes(), 'all response tables normal/-O')
    need(digest(HERE/'response_tables.json') == response['tables_sha256'], 'response table input')
    need(response['actual_compiled_queries'] == 49976 and response['exact_shared_cubic_checks'] == 1056, 'full named query counts')
    need(len(read('response_tables.json')['cases']) == len(response['cases']) == 4, 'prime control inventory')
    for n in (8, 12):
        scout = read('sat'+str(n)+'.json')
        need(scout['status'].startswith('UNKNOWN') and scout['rows'] is None, 'SAT has no conclusion')
        need(digest(HERE/'sat_scout.py') == scout['source_sha256'], 'scout source')
    if args.replay:
        with tempfile.TemporaryDirectory(prefix='prime12-replay-') as folder:
            for optimized in (False, True):
                prefix = [sys.executable]+(['-O'] if optimized else [])
                for stem in ('current', 'switch', 'response', 'audit'):
                    out = Path(folder)/(stem+'.json')
                    command = prefix+[str(HERE/(stem+'.py')), '--output', str(out)]
                    if stem == 'response':
                        tables = Path(folder)/'tables.json'
                        command += ['--tables', str(tables)]
                    subprocess.run(command, check=True)
                    need(out.read_bytes() == (HERE/(stem+'.json')).read_bytes(), 'replay: '+stem)
                    if stem == 'response':
                        need(tables.read_bytes() == (HERE/'response_tables.json').read_bytes(), 'full table replay')
    print(json.dumps(dict(status='PASS', files=len(inventory),
                          dependencies=len(manifest['dependencies_lf_sha256']),
                          orbit=2048, strong_prime_presentations=2024,
                          complete_bank=198947, full_named_queries=49976,
                          exact_shared_cubics=1056, receipt_pairs=4,
                          replayed=args.replay,
                          scope='One complete orbit excluded; universal prime12 and missing-spectrum goals OPEN'),
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
