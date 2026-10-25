"""Hash-pinned replay of the five declared exact research checkers."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CHECKERS = ('anchor_check.py', 'endpoint_check.py', 'gap_check.py',
            'gap_owner_audit.py', 'f1_check.py')
DEPENDENCIES = (
    '05-knowledge/results/metagraph_growth_20261020/enumeration.py',
    '05-knowledge/results/zero_heavy_ports_20261020g/zero_heavy.md',
    '05-knowledge/results/zero_heavy_ports_20261020g/microcosm.md',
    '05-knowledge/results/charge_palettes_boundary_20261020g/zero_heavy.md',
    '05-knowledge/results/zero_heavy_next_20261020h/README.md',
    '05-knowledge/results/defect_boundary_20261020/sparse_theory.md',
    '05-knowledge/results/endpoint_three_repairs_20261020h/README.md',
    '05-knowledge/results/endpoint_three_order10_20261021/README.md',
    '05-knowledge/results/endpoint_four_frontier_20261021/README.md',
    '05-knowledge/results/endpoint_four_repairs_20261021/README.md',
    '05-knowledge/results/session58_hourglass_complement_geometry_20261025/absolute-geometry/ABSOLUTE_GEOMETRY_BRIDGES.md',
    '05-knowledge/results/session58_hourglass_complement_geometry_20261025/projective-support-bridge/PROJECTIVE_SUPPORT_BRIDGE.md',
    '01-canon/theorems/THM-4734-endpoint-three-layer-closed-at-order-ten-and-the-endpoint-four-layer-at-orders-nine-and-ten.md',
    '01-canon/theorems/THM-4736-the-delta-equals-k-zero-heavy-layer-exact-counts-missing-gaps-and-the-adjacent-matched-pair-collapse.md',
)


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def norm(data):
    return data.replace(b'\r\n', b'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', action='store_true')
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    paths = [p for p in HERE.iterdir() if p.is_file() and
             p.suffix in ('.py', '.md', '.out', '.json') and
             p.name not in ('manifest.json', 'verify.out')]
    paths += [ROOT / p for p in DEPENDENCIES]
    pins = {p.relative_to(ROOT).as_posix(): hashlib.sha256(norm(p.read_bytes())).hexdigest()
            for p in sorted(set(paths))}
    manifest = HERE / 'manifest.json'
    if args.record:
        manifest.write_text(json.dumps({'normalization': 'CRLF to LF', 'files': pins},
                                       indent=2, sort_keys=True)+'\n', encoding='utf-8')
    need(json.loads(manifest.read_text(encoding='utf-8'))['files'] == pins, 'changed pins')
    print('Pinned source/input/receipt files:', len(pins), 'PASS', flush=True)
    for name in CHECKERS:
        tree = ast.parse((HERE / name).read_text(encoding='utf-8'))
        need(not any(isinstance(n, ast.Assert) for n in ast.walk(tree)), 'removable assertion '+name)
    print('Required checkers without removable assertions:', len(CHECKERS), flush=True)
    links = 0
    for path in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)', path.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):
                continue
            need((path.parent / target.split('#')[0]).resolve().exists(), 'broken link '+target)
            links += 1
    print('Local Markdown targets:', links, 'PASS', flush=True)
    if args.replay:
        for name in CHECKERS:
            path = HERE / name
            expected = norm(path.with_suffix('.out').read_bytes())
            for flags in ([], ['-O']):
                result = subprocess.run([sys.executable, *flags, str(path)], cwd=ROOT,
                                        capture_output=True, check=False)
                need(result.returncode == 0, name+' failed: '+result.stderr.decode(errors='replace'))
                need(norm(result.stdout) == expected, name+' receipt mismatch '+str(flags))
            print('Normal/optimized replay:', name, 'PASS', flush=True)
    print('PASS', flush=True)


if __name__ == '__main__':
    main()
