#!/usr/bin/env python3
"""Pin referenced repo sources and verify this packet's artifact hashes."""
from pathlib import Path
import argparse,hashlib,json,re,subprocess
HERE=Path(__file__).resolve().parent
ROOT=Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=HERE,text=True).strip()).resolve()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)
def digest(data):return hashlib.sha256(data).hexdigest()
def files():
    return sorted(p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json')
def main(write):
    if write:
        revision=git('rev-parse','HEAD').decode().strip();sources=set()
        for doc in HERE.glob('*.md'):
            for target in re.findall(r'\]\(([^)]+)\)',doc.read_text()):
                if '://' in target or target.startswith('#'):continue
                path=(doc.parent/target.split('#')[0]).resolve()
                if path==HERE or HERE in path.parents:continue
                sources.add(path.relative_to(ROOT).as_posix())
        sources.add('05-knowledge/results/square_sum_archive_stability_20261102/verify_reduction.py')
        rows=[]
        for path in sorted(sources):
            data=git('show',revision+':'+path)
            rows.append({'path':path,'git_blob':git('rev-parse',revision+':'+path).decode().strip(),'sha256':digest(data)})
        pins={'source_revision':revision,'scope':'Referenced proof inputs and recovered context; pinning is not an assertion that every historical theorem was independently reaudited.',
              'correction_note':'THM-461 is pinned after its repaired source-count and real-root statements; MISTAKES.md records their lineage. The honeycomb subtraction correction is also explicit there.',
              'sources':rows}
        (HERE/'SOURCE_PINS.json').write_text(json.dumps(pins,indent=2,sort_keys=True)+'\n')
        manifest={'scope':'Integrity only, not a proof of mathematical correctness. Shared navigation/correction edits are versioned by Git, not frozen here.',
                  'files':{p.relative_to(HERE).as_posix():digest(p.read_bytes()) for p in files()}}
        (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    manifest=json.loads((HERE/'MANIFEST.json').read_text())
    actual={p.relative_to(HERE).as_posix():digest(p.read_bytes()) for p in files()}
    if actual!=manifest['files']:raise RuntimeError('packet file inventory/hash mismatch')
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for row in pins['sources']:
        if digest(git('show',pins['source_revision']+':'+row['path']))!=row['sha256']:raise RuntimeError('source pin mismatch '+row['path'])
    print(json.dumps({'status':'passed','packet_files':len(actual),'pinned_sources':len(pins['sources']),'source_revision':pins['source_revision']},sort_keys=True))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');main(ap.parse_args().write)
