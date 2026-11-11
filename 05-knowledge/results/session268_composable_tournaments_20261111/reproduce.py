#!/usr/bin/env python3
"""Replay only this session's new finite controls. Standard library only."""
import ast,base64,gzip,hashlib,json,os,pathlib,shutil,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent

def require(x,m):
    if not x: raise RuntimeError(m)

def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items() if k!='elapsed_seconds'}
    if isinstance(x,list):return [clean(v) for v in x]
    return x

def main():
    summary=json.loads((ROOT/'modular/receipt.json').read_text())
    raw=gzip.decompress(base64.b64decode((ROOT/'modular/full_receipt.json.gz.b64').read_bytes()))
    require(hashlib.sha256(raw).hexdigest()==summary['full_receipt_sha256'],'full receipt hash')
    full=json.loads(raw)
    source_checks=[]
    for f in sorted(ROOT.rglob('*.py')):
        parsed=ast.parse(f.read_text());require(not any(isinstance(n,ast.Assert) for n in ast.walk(parsed)),str(f)+' contains optimization-sensitive assert')
        source_checks.append(str(f.relative_to(ROOT)))
    cases=[('automata/probe.py',[], 'automata/receipt.json'),
      ('automata/interacting_ca.py',[],'automata/ca_receipt.json'),
      ('dynamics/replay.py',[],None),
      ('modular/independent_controls.py',[],'modular/controls.json'),
      ('modular/census.py',['--max-n','16','--out','modular-replay.json'],'modular-replay.json'),
      ('spectrum/replay.py',[],'spectrum/results.json'),
      ('spectrum/compile_target.py',['--measure','--output','compiler-replay.json'],'compiler-replay.json')]
    receipts=[]
    with tempfile.TemporaryDirectory(prefix='session268-replay-') as td:
        work=pathlib.Path(td)/'packet';shutil.copytree(ROOT,work,ignore=shutil.ignore_patterns('__pycache__'))
        for file,args,output in cases:
            modes=[]
            for opt in (False,True):
                cmd=[sys.executable,'-B']+(['-O'] if opt else [])+[file]+args
                proc=subprocess.run(cmd,cwd=work,text=True,capture_output=True,check=True)
                data=json.loads((work/output).read_text() if output else proc.stdout)
                modes.append(clean(data))
            require(modes[0]==modes[1],file+' normal/optimized mismatch')
            expected= full if file=='modular/census.py' else json.loads((ROOT/('spectrum/compiler_results.json' if file=='spectrum/compile_target.py' else 'dynamics/receipt.json' if output is None else output)).read_text())
            require(modes[0]==clean(expected),file+' recorded result mismatch')
            receipts.append({'program':file,'normal':'PASS','optimized':'PASS','recorded_result':'MATCH'})
    # Spatial partition recurrence: a direct matrix computation is independent of radicals.
    W=((215,179),(245,201));M=((1,0),(0,1));previous=[2,416]
    for n in range(1,21):
        M=tuple(tuple(sum(M[i][k]*W[k][j] for k in range(2)) for j in range(2)) for i in range(2))
        if n>=2:previous.append(416*previous[-1]+640*previous[-2])
        require(M[0][0]+M[1][1]==previous[n],'partition recurrence')
    print(json.dumps({'status':'PASS','scope':'Seven new programs normal and optimized; no inherited Lean or repository-wide audit','full_census_sha256':summary['full_receipt_sha256'],'programs':receipts,'assert_free_sources':source_checks,'partition_recurrence_controls':20,'spatial_Z2':previous[2]},indent=2))
if __name__=='__main__':main()
