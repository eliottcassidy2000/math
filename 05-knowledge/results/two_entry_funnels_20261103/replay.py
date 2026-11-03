"""Freeze normal/-O receipts and source pins without removable assertions."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def main():
    result=dict(status='PASS',checks={},local_hashes={},inherited_hashes={},
                hash_convention='SHA256 after CRLF-to-LF normalization')
    for script,output in (('verify.py','verification.json'),('carrier.py','carrier.json'),('boundary.py','boundary.json'),
                          ('enumeration.py','enumeration.json'),('research.py','research.out')):
        runs=[];search_receipts=[]
        for optimized in (False,True):
            args=[sys.executable,'-B','-X','utf8']+(['-O'] if optimized else [])+[str(HERE/script)]
            p=subprocess.run(args,check=True,capture_output=True,text=True,encoding='utf-8')
            runs.append(p.stdout.replace('\r\n','\n'))
            if script=='research.py':search_receipts.append((HERE/'template_search.json').read_text(encoding='utf-8'))
        if runs[0]!=runs[1] or (search_receipts and search_receipts[0]!=search_receipts[1]):
            raise RuntimeError(('normal/optimized disagreement',script))
        (HERE/output).write_text(runs[0],encoding='utf-8',newline='\n')
        result['checks'][script]=dict(normal_optimized_identical=True,output=output,output_sha256=digest(HERE/output))
    for name in ('README.md','TWO_ENTRY.md','READY.md','TRANSPORT.md','TILINGS.md','AUDIT.md',
                 'research.py','verify.py','carrier.py','boundary.py','enumeration.py','replay.py'):
        result['local_hashes'][name]=digest(HERE/name)
    for rel in ('polarized_path_calculus_20261103/explore.py','polarized_path_calculus_20261103/recognize.py',
                'source_selection_closure_20261103/graded_carrier.py','golden_inverse_defects_20261103/inverse.py',
                'golden_seam_recurrence_20261103/atom_atlas.json','adaptive_metagraph_reformulation_20261103/metagraph_probe.py'):
        result['inherited_hashes'][rel]=digest(HERE.parent/rel)
    result['checks']['template_search.json']=dict(output_sha256=digest(HERE/'template_search.json'))
    (HERE/'manifest.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('PASS: five normal/-O jobs, explicit joint/native/phase controls, source hashes saved.')


if __name__=='__main__':main()
