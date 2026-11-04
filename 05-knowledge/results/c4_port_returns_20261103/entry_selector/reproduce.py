"""Reproduce and pin only the new entry-aware selector packet."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def need(ok,message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def main():
    outputs=[];receipts=[]
    for optimized in (False,True):
        command=[sys.executable]+(['-O'] if optimized else [])+[str(HERE/'check_selector.py')]
        result=subprocess.run(command,check=True,capture_output=True)
        need(not result.stderr,'unexpected selector stderr')
        outputs.append(result.stdout.replace(b'\r\n',b'\n'))
        receipts.append([(HERE/name).read_bytes() for name in ('receipt.json','profiles.json')])
    need(outputs[0]==outputs[1] and receipts[0]==receipts[1],'normal/-O selector mismatch')
    for name,data in zip(('normal.out','optimized.out'),outputs):
        (HERE/name).write_bytes(data)
    files=['PROOF.md','check_selector.py','reproduce.py','receipt.json','profiles.json','normal.out','optimized.out']
    old=HERE.parent/'reader'
    old_manifest=json.loads((old/'manifest.json').read_text())
    need(all(digest(old/name)==expected for name,expected in old_manifest['files'].items()),'pinned old reader packet changed')
    report=dict(normal_optimized_identical=True,previous_reader_packet_unchanged=True,
                files={name:digest(HERE/name) for name in files},
                scope='Three finite regenerated positive entries with actual native repairs/direct C2 children; six scoped nonneutral triangle restorations.')
    (HERE/'manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: entry-aware programs normal/-O byte-identical; old reader pins unchanged')


if __name__=='__main__':
    main()
