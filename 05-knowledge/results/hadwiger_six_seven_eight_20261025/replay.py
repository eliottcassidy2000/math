"""Replay four producers and four independent audits against their receipts."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
CASES=[
 ('metric_spindle/check.py','metric_spindle/normal.json'),
 ('projective_eight/check.py','projective_eight/normal.json'),
 ('fano_eight/check.py','fano_eight/receipt.json'),
 ('cube_boundary/check.py','cube_boundary/normal.json'),
 ('audits/check_metric_independent.py','audits/metric.normal.json'),
 ('audits/check_projective_independent.py','audits/projective_independent.json'),
 ('audits/check_palette_independent.py','audits/palette-normal.json'),
 ('audits/check_boundary_independent.py','audits/boundary_independent.json'),
]

def replay(case):
    script,saved,expected=case
    outputs=[]
    for flags in ([],['-O']):
        proc=subprocess.run([sys.executable,'-B',*flags,str(HERE/script)],
                            cwd=HERE,capture_output=True,timeout=300)
        if proc.returncode:
            raise RuntimeError(f'{script} {flags}: {proc.stderr.decode(errors="replace")}')
        outputs.append(proc.stdout)
    if outputs[0]!=outputs[1]:
        raise RuntimeError(f'normal/optimized byte mismatch: {script}')
    if json.loads(outputs[0])!=expected:
        raise RuntimeError(f'frozen receipt mismatch: {script}')
    return {'script':script,'receipt':saved,'status':'PASS',
            'source_sha256':hashlib.sha256((HERE/script).read_bytes()).hexdigest(),
            'normal_optimized_identical':True,
            'receipt_data_identical':True}

if __name__=='__main__':
    # Snapshot before running: two inherited checkers regenerate a JSON
    # artifact. Comparing only after execution would be tautological.
    cases=[(s,r,json.loads((HERE/r).read_text(encoding='utf-8-sig'))) for s,r in CASES]
    artifacts=['cube_boundary/receipt.json','audits/projective_independent.json']
    frozen={p:json.loads((HERE/p).read_text(encoding='utf-8-sig')) for p in artifacts}
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows=list(pool.map(replay,cases[:4]))
        # The boundary audit reads the producer's coordinate artifact.
        rows+=list(pool.map(replay,cases[4:]))
    for p,data in frozen.items():
        if json.loads((HERE/p).read_text(encoding='utf-8-sig'))!=data:
            raise RuntimeError(f'regenerated artifact changed: {p}')
    print(json.dumps({'status':'PASS','script_pairs':len(rows),'results':rows},
                     sort_keys=True,indent=2))
