"""Hash the owned packet and its explicitly cited/read code dependencies."""
from pathlib import Path
import argparse, hashlib, json, re
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    manifest=HERE/'manifest.json'
    if args.check:
        data=json.loads(manifest.read_text())
        for rel,expected in data['files'].items():
            if digest(ROOT/rel)!=expected:raise RuntimeError(('hash mismatch',rel))
        print('Manifest PASS:',len(data['files']),'files')
        return
    files={p for p in HERE.iterdir() if p.is_file() and p!=manifest}
    urls=set()
    for md in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',md.read_text(encoding='utf-8')):
            if target.startswith(('http:','https:')):
                urls.add(target);continue
            target=target.split('#')[0]
            if not target:continue
            p=(md.parent/target).resolve()
            if not p.is_file():raise RuntimeError(('missing local reference',md.name,target))
            files.add(p)
    for rel in (
        '05-knowledge/results/port_residual_integral_20261025/precision/check.py',
        '05-knowledge/results/primegame_port_spectrum_20261103/port_scout.py',
        '05-knowledge/results/primegame_port_spectrum_20261103/verify.py',
        '05-knowledge/results/primegame_port_spectrum_20261103/growth_identity.py',
        '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
        '04-computation/almost_fixed_frame_s533.py',
    ):files.add(ROOT/rel)
    data={'starting_revision':'e30c89e9fcdc',
          'incoming_reviewed':['090fcaa1e7','490da01fbb','9b75ef18bf','a28303347b','947b6c3852'],
          'scope':'Owned artifacts and direct local proof/code references. Self-audit, not external review.',
          'external_sources':sorted(urls),
          'files':{p.relative_to(ROOT).as_posix():digest(p) for p in sorted(files)}}
    manifest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print('Manifest sealed:',len(data['files']),'files')
if __name__=='__main__':main()

