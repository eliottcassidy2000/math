"""Closed inventory, source pins, paired checks; --freeze is an explicit author action."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
PAIRS=[('check_moments.py','moments'),('check_lifts.py','lifts'),
       ('check_resolved.py','resolved'),('audit_independent.py','audit')]
SOURCE_PATHS=[
 '05-knowledge/results/normseven_f1_20261025/README.md',
 '05-knowledge/results/normseven_f1_20261025/KUMMER_AND_F1.md',
 '05-knowledge/results/normseven_f1_20261025/TWIST_AVERAGE.md',
 '05-knowledge/results/normseven_f1_20261025/PORTS_AND_CHARACTERS.md',
 '05-knowledge/results/normseven_f1_20261025/INTEGRAL_WRAP.md',
 '05-knowledge/results/normseven_boundary_20261025/BOUNDARY.md',
 '05-knowledge/results/normseven_boundary_20261025/STACK_GLUE.md',
 '05-knowledge/results/six_seven_geometry_20261025/README.md',
 '05-knowledge/results/six_seven_gale_cube_20261025/README.md',
 '05-knowledge/results/twist_conductor_reverse_20261025/FIXED_TWISTS.md',
 '05-knowledge/results/session63_integral_observers_root_arrangements_20261025/twists/RESULT.md',
 '05-knowledge/results/session63_integral_observers_root_arrangements_20261025/twists/ZETA_AND_EFFECTIVITY.md',
 '01-canon/theorems/THM-4765-quadratic-twists-of-a-square-root-cover-average-to-the-branch-complement.md',
 '01-canon/theorems/THM-4773-fixed-twist-rational-lifts-on-the-hourglass-boundary-central-conic-outer-genus-five.md',
 '05-knowledge/results/normseven_arithmetic_f1_20261026/README.md',
 '05-knowledge/results/normseven_arithmetic_f1_20261026/boundary/RESULT.md',
 '05-knowledge/results/normseven_arithmetic_f1_20261026/reduction/RESULT.md',
 '05-knowledge/results/normseven_arithmetic_f1_20261026/f1/RESULT.md',
 '05-knowledge/results/normseven_arithmetic_f1_20261026/gale/RESULT.md',
]


def need(ok,label):
    if not ok: raise RuntimeError(label)


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory():
    return sorted(p for p in BASE.iterdir() if p.is_file() and p.name!='MANIFEST.json')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--replay',action='store_true')
    ap.add_argument('--freeze',action='store_true')
    args=ap.parse_args()
    if args.freeze:
        # Pin only the inspected inherited revisions. Do not use this to hide
        # a source change after a rebase without reading the mathematical diff.
        sources={'initial_revision':'983768efde8483643fc3c6ba488eee04622ebd2d',
                 'incoming_revision':'2daf5e00af',
                 'files':{p:sha(ROOT/p) for p in SOURCE_PATHS},
                 'external_conventions_checked_2026_10_26':[
                    'https://stacks.math.columbia.edu/tag/04UV',
                    'https://arxiv.org/html/math/0404185v7']}
        (BASE/'SOURCES.json').write_text(json.dumps(sources,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
        for path in inventory():
            path.write_text(path.read_text(encoding='utf-8'),encoding='utf-8',newline='\n')
        manifest={'files':{p.name:sha(p) for p in inventory()}}
        (BASE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    manifest=json.loads((BASE/'MANIFEST.json').read_text())
    need(set(manifest['files'])=={p.name for p in inventory()},'closed packet inventory')
    for name,digest in manifest['files'].items(): need(sha(BASE/name)==digest,'artifact hash: '+name)
    sources=json.loads((BASE/'SOURCES.json').read_text())
    for path,digest in sources['files'].items(): need(sha(ROOT/path)==digest,'source hash: '+path)
    for script,stem in PAIRS:
        normal=(BASE/(stem+'.normal.json')).read_bytes()
        need(normal==(BASE/(stem+'.optimized.json')).read_bytes(),'paired bytes: '+stem)
        json.loads(normal)
        if args.replay:
            for flags in [[],['-O']]:
                proc=subprocess.run([sys.executable,'-B']+flags+[str(BASE/script)],cwd=BASE,
                                    capture_output=True,text=True,encoding='utf-8',check=True)
                # Universal-newline decoding makes the saved LF receipts portable.
                need(proc.stdout.encode('utf-8')==normal,'replay mismatch: '+script+str(flags))
                need(not proc.stderr,'unexpected stderr: '+script)
    print(json.dumps({'status':'PASS','packet_files':len(manifest['files']),
                      'source_pins':len(sources['files']),'paired_receipts':len(PAIRS),
                      'normal_and_optimized_replays':2*len(PAIRS) if args.replay else 0},sort_keys=True))


if __name__=='__main__': main()
