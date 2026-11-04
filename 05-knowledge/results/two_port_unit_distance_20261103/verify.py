"""Check frozen artifacts; optional replay never overwrites expectations."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ENGINES=(
 ('reconstruction/check.py','reconstruction/normal.json','reconstruction/optimized.json'),
 ('algebra/check.py','algebra/receipt.json','algebra/receipt.optimized.json'),
 ('geometry/verify.py','geometry/receipt.json','geometry/receipt.optimized.json'),
 ('geometry/audit_colour.py','geometry/audit_colour.json','geometry/audit_colour.optimized.json'),
 ('census/check.py','census/python8.json','census/python8.optimized.json'),
 ('census/run_cpp.py','census/cpp9.json','census/cpp9.optimized.json'),
 ('census/audit.py','census/audit.json','census/audit.optimized.json'),
)


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def pins():
    manifest=json.loads((HERE/'MANIFEST.json').read_text(encoding='utf8'))['files']
    actual={p.relative_to(HERE).as_posix() for p in HERE.rglob('*') if p.is_file()
            and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='MANIFEST.json'}
    require(set(manifest)==actual,'artifact inventory mismatch')
    for path,digest in manifest.items():require(sha(HERE/path)==digest,'artifact changed: '+path)
    sources=json.loads((HERE/'SOURCE_PINS.json').read_text(encoding='utf8'))['files']
    for path,digest in sources.items():require(sha(ROOT/path)==digest,'inherited source changed: '+path)
    return len(manifest),len(sources)


def main():
    require(sys.argv[1:] in ([],['--replay']),'usage: verify.py [--replay]')
    counts=pins()
    # Freeze every expectation before any producer runs.
    expected={script:(HERE/normal).read_text(encoding='utf8') for script,normal,opt in ENGINES}
    for script,normal,opt in ENGINES:
        require(expected[script]==(HERE/opt).read_text(encoding='utf8'),'normal/optimized receipt mismatch: '+script)
        json.loads(expected[script])
    a=json.loads(expected['census/check.py']);b=json.loads(expected['census/run_cpp.py'])
    require(a['orders']==b['orders'][:8],'Python/C++ census disagreement')
    require(not a['collisions'] and b['collision_count']==0,'unexpected collision')
    require(b['orders'][-1]['classes']==b['orders'][-1]['sidecar_classes']==191536,'order9 incomplete')
    require(all(r['automorphism_group_disagreements']==0 for r in b['orders']),'extra response automorphisms')
    replayed=0
    if sys.argv[1:]==['--replay']:
        for script,normal,opt in ENGINES:
            path=HERE/script
            for flag in ([],['-O']):
                result=subprocess.run([sys.executable,'-B']+flag+[str(path)],cwd=path.parent,
                                      capture_output=True,text=True,encoding='utf8',check=True)
                require(result.stdout==expected[script],'frozen replay mismatch: '+script)
                # census/check.py emits bounded progress to stderr by design.
                require(not result.stderr or script=='census/check.py','unexpected stderr: '+script)
                replayed+=1
                print('PASS '+script+(' -O' if flag else ''),flush=True)
    require(pins()==counts,'pins changed during verification')
    print(json.dumps({'status':'PASS','artifact_files':counts[0],'inherited_source_files':counts[1],
                      'matching_receipt_pairs':len(ENGINES),'replayed_engines':replayed,
                      'expected_receipts_modified':False},sort_keys=True))


if __name__=='__main__':main()
