"""Normal/optimized receipt equality; independent engines remain separate."""
from pathlib import Path
import sys,subprocess,gzip,json,hashlib,time
HERE=Path(__file__).resolve().parent
JOBS=[
 ('head_factor_search.py',[],'head_factor_search'),
 ('endpoint_profiles.py',[],'endpoint_profiles'),
 ('endpoint_profiles.py',['split'],'split_endpoint_profiles'),
 ('endpoint_profiles.py',['search','3'],'endpoint_search_3'),
 ('endpoint_profiles.py',['search','5'],'endpoint_search_5'),
 ('endpoint_profiles.py',['search','5','split'],'split_endpoint_search_5'),
 ('endpoint_periods.py',['3'],'endpoint_periods'),
 ('observer_rewrites.py',[],'observer_rewrites'),
 ('synchronized_descent.py',[],'synchronized_descent'),
 ('synchronized_repairs.py',[],'synchronized_repairs'),
 ('independent_descent.py',[],'independent_descent'),
 ('descent_grammar.py',[],'descent_grammar'),
 ('modular_bridge.py',[],'modular_bridge'),
 ('connection_checks.py',[],'connection_checks'),
 ('verify.py',[],'verification')]

def run():
    results=[];start=time.time();resume=sys.argv[2] if len(sys.argv)>2 and sys.argv[1]=='--from' else None;active=resume is None
    for script,args,name in JOBS:
        active=active or name==resume
        if active:
            for flags in [[],['-O']]:
                subprocess.run([sys.executable,'-X','utf8',*flags,str(HERE/script),*args],check=True)
        paths=[HERE/(name+tag+'.json.gz') for tag in ['', '.optimized']]
        if paths[0].read_bytes()!=paths[1].read_bytes():raise RuntimeError('normal/optimized mismatch: '+name)
        results.append({'receipt':name,'sha256':hashlib.sha256(paths[0].read_bytes()).hexdigest()})
        print('PAIRED',name,'elapsed',round(time.time()-start,1),flush=True)
    if not active:raise RuntimeError('unknown resume point')
    out={'status':'PASS','scope':'Fifteen producer/receipt pairs normal and optimized. Independent methods are those described in independent_descent and verify; same-code optimization agreement is not an independent proof.','pairs':results}
    (HERE/'reproduction.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':run()
