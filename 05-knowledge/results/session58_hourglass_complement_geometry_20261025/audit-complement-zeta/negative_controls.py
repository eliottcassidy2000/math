from pathlib import Path
import subprocess,json,sys,ast
root=Path(__file__).resolve().parent
rows=[]
for optimized in (False,True):
    args=[sys.executable]+(['-O'] if optimized else [])
    r=subprocess.run(args+[str(root/'independent_audit.py'),'--negative-exception'],capture_output=True,text=True)
    if r.returncode==0 or 'injected negative exception' not in r.stderr: raise RuntimeError('independent negative control failed')
    rows.append({'target':'independent exception','optimized':optimized,'exit':r.returncode,'detected':True})
    src=(root/'replays/complement/verify.py').read_text().replace("require([r['tournaments'] for r in rows]==[1,1,2,4,12],'small tournament counts')","require([r['tournaments'] for r in rows]==[1,1,2,4,13],'small tournament counts')")
    target=root/'replays/complement/negative.py';target.write_text(src)
    r=subprocess.run(args+[str(target)],capture_output=True,text=True)
    if r.returncode==0 or 'small tournament counts' not in r.stderr: raise RuntimeError('producer negative control failed')
    rows.append({'target':'complement perturbed order-five count','optimized':optimized,'exit':r.returncode,'detected':True})
    src=(root/'replays/absolute/verify_bridges.py').read_text().replace('assert (N-f)%2==0','assert False, "injected broken assertion"')
    target=root/'replays/absolute/negative.py';target.write_text(src)
    r=subprocess.run(args+[str(target)],capture_output=True,text=True)
    if bool(r.returncode)==optimized: raise RuntimeError('unexpected original assertion behavior')
    rows.append({'target':'original absolute injected false assertion','optimized':optimized,'exit':r.returncode,'detected':bool(r.returncode),'known_defect':True})
for f in root.glob('replays/*/negative.py'):f.unlink()
(root/'negative_results.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
