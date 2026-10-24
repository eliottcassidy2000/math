import json,shutil,subprocess,sys,tempfile
from pathlib import Path
src=Path(__file__).resolve().parent
results={}
with tempfile.TemporaryDirectory() as tmp:
    for kind in ['corrupt','extra','missing']:
        dst=Path(tmp)/kind
        shutil.copytree(src,dst)
        if kind=='corrupt':
            with (dst/'REPORT.md').open('a') as f: f.write('\ncorruption\n')
        if kind=='missing': (dst/'normal.json').unlink()
        if kind=='extra':
            (dst/'nested').mkdir()
            (dst/'nested'/'extra.txt').write_text('unexpected')
        run=subprocess.run([sys.executable,str(dst/'replay.py')],capture_output=True)
        if run.returncode==0: raise RuntimeError('hostile accepted '+kind)
        results[kind]='REJECTED'
print(json.dumps(results,sort_keys=True))
