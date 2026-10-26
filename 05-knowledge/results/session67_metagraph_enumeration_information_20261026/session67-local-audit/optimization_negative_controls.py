#!/usr/bin/env python3
"""Verify optimized Python rejects corrupted production node/line data."""
import argparse,csv,json,shutil,subprocess,sys,tempfile
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('production_output',type=Path);a=ap.parse_args()
here=Path(__file__).resolve().parent;receipts=[]
for kind in ('node_mass','line_multiplicity'):
    with tempfile.TemporaryDirectory(prefix='negative-',dir=here) as tmp:
        out=Path(tmp)
        for n in range(1,5):
            for prefix in ('nodes','lines'):shutil.copyfile(a.production_output/f'{prefix}_{n}.tsv',out/f'{prefix}_{n}.tsv')
        path=out/('nodes_4.tsv' if kind=='node_mass' else 'lines_4.tsv')
        rows=list(csv.DictReader(path.open(),delimiter='\t'));fields=list(rows[0])
        key='tilings_tc' if kind=='node_mass' else 'multiplicity'
        rows[0][key]=str(int(rows[0][key])+1)
        with path.open('w') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
        run=subprocess.run([sys.executable,'-O',str(here/'independent_check.py'),'--max-n','4','--compare-dir',str(out),'--output',str(out/'unexpected.json')],capture_output=True,text=True)
        if run.returncode==0 or 'RuntimeError' not in run.stderr:raise RuntimeError('corruption was not rejected: '+kind)
        receipts.append({'corruption':kind,'python_optimized':True,'return_code':run.returncode,'rejected':True,'last_error_line':run.stderr.strip().splitlines()[-1]})
(here/'optimization_negative_controls.json').write_text(json.dumps({'status':'PASS','controls':receipts},indent=2)+'\n')
print(json.dumps(receipts,indent=2))
