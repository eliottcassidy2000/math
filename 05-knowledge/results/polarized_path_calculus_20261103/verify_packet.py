"""Replay both implementations under normal and optimized Python."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    record={'checks':{},'source_hashes':{},'inherited_hashes':{}}
    for script,output in (('explore.py','explore.out'),('counts.py','counts.json'),('recognize.py','recognition.json')):
        values=[]
        for optimize in (False,True):
            args=[sys.executable,'-B','-X','utf8']+(['-O'] if optimize else [])+[str(HERE/script)]
            proc=subprocess.run(args,check=True,capture_output=True,text=True,encoding='utf-8')
            values.append(proc.stdout.replace('\r\n','\n'))
        if values[0]!=values[1]:
            raise RuntimeError(('optimized mismatch',script))
        (HERE/output).write_text(values[0],encoding='utf-8',newline='\n')
        record['checks'][script]={'normal_optimized_identical':True,'output_sha256':digest(HERE/output)}
    for name in ('explore.py','counts.py','recognize.py','verify_packet.py','README.md','FUNNELS.md','CALCULUS.md','AUDIT.md'):
        record['source_hashes'][name]=digest(HERE/name)
    for rel in ('golden_inverse_defects_20261103/inverse.py','golden_seam_recurrence_20261103/atom_atlas.json'):
        record['inherited_hashes'][rel]=digest(HERE.parent/rel)
    (HERE/'verification.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('All three programs agree under normal and optimized Python; receipts and hashes saved.')


if __name__=='__main__':
    main()
