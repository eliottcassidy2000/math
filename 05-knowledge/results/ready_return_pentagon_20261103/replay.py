"""Reproduce the constructive atlas and source returns under normal/-O Python."""
from pathlib import Path
import hashlib,json,subprocess,sys

HERE=Path(__file__).resolve().parent


def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def main():
    result=dict(status='PASS',checks={},local_hashes={},inherited_hashes={},
                hash_convention='SHA256 after CRLF-to-LF normalization')
    for script,output in (('explore_returns.py','languages.out'),('obstruction.py','obstruction.json'),
                          ('return_program.py','returns.json'),('ready_inverse.py','inverse.json')):
        runs=[];languages=[]
        for optimized in (False,True):
            args=[sys.executable,'-B','-X','utf8']+(['-O'] if optimized else [])+[str(HERE/script)]
            p=subprocess.run(args,check=True,capture_output=True,text=True,encoding='utf-8')
            runs.append(p.stdout.replace('\r\n','\n'))
            if script=='explore_returns.py':languages.append((HERE/'languages.json').read_text(encoding='utf-8'))
        if runs[0]!=runs[1] or (languages and languages[0]!=languages[1]):raise RuntimeError(('normal/-O mismatch',script))
        (HERE/output).write_text(runs[0],encoding='utf-8',newline='\n')
        result['checks'][script]=dict(normal_optimized_identical=True,receipt=output,sha256=digest(HERE/output))
    for name in ('README.md','OBSTRUCTION.md','ATLAS.md','RETURN.md','PENTAGON.md','AUDIT.md',
                 'explore_returns.py','obstruction.py','return_program.py','ready_inverse.py','replay.py'):
        result['local_hashes'][name]=digest(HERE/name)
    for rel in ('two_entry_funnels_20261103/research.py','two_entry_funnels_20261103/verify.py',
                'two_entry_funnels_20261103/carrier.py','two_entry_funnels_20261103/READY.md',
                'two_entry_funnels_20261103/TWO_ENTRY.md','polarized_path_calculus_20261103/explore.py',
                'polarized_path_calculus_20261103/recognize.py','source_selection_closure_20261103/graded_carrier.py',
                'golden_inverse_defects_20261103/inverse.py','golden_seam_recurrence_20261103/atom_atlas.json',
                'adaptive_metagraph_reformulation_20261103/metagraph_probe.py'):
        result['inherited_hashes'][rel]=digest(HERE.parent/rel)
    result['checks']['languages.json']=dict(sha256=digest(HERE/'languages.json'))
    (HERE/'manifest.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('PASS: four normal/-O jobs, exact ready inverse, full joint returns, finite atlas, phase and precision controls.')


if __name__=='__main__':main()
