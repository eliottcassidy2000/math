"""Reject weakened/incomplete certificates using this auditor's own validator."""
import contextlib, io, json, tempfile, sys
from pathlib import Path
from check import run
src=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent.parent/'lrc/triple_receipt.json'
original=json.loads(src.read_text())
checks=[]
for name in ('missing_core','missing_first_branch','missing_second_branch','false_second_width','false_cutoff'):
    d=json.loads(src.read_text())
    if name=='missing_core': d['families'].pop()
    elif name=='missing_first_branch': d['families'][0]['branches'].pop(0)
    elif name=='missing_second_branch': d['families'][0]['branches'][0]['finite_branches'].pop(0)
    elif name=='false_second_width': d['families'][0]['branches'][0]['finite_branches'][0][1]='1/100'
    else: d['families'][0]['branches'][0]['finite_branches'][0][2]+=1
    with tempfile.TemporaryDirectory() as folder:
        path=Path(folder)/'mutant.json'; path.write_text(json.dumps(d))
        try: run(path)
        except RuntimeError as error: checks.append({'mutation':name,'rejected_by':str(error)})
        else: raise RuntimeError('Mutation incorrectly accepted: '+name)
print(json.dumps({'status':'PASS','mutations':checks},sort_keys=True))
