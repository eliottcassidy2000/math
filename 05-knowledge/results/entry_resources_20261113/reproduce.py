"""Recompute every producer normally and with optimization; compare bytes."""
from pathlib import Path
import sys,subprocess,gzip,json,hashlib,time
HERE=Path(__file__).resolve().parent
JOBS=[
 ('cut_gadgets.py',[],'cut_gadgets'),
 ('uniform_gadgets.py',['5'],'uniform_gadgets'),
 ('port_gadgets.py',[],'port_gadgets'),
 ('graceful_padding.py',[],'graceful_padding'),
 ('square_covers.py',[],'square_covers'),
 ('square_cut_audit.py',[],'square_cut_audit'),
 ('square_endpoint_certificate.py',[],'square_endpoint_certificate'),
 ('prime_phase_machine.py',[],'prime_phase_machine'),
 ('moving_root.py',[],'moving_root'),
 ('weighted_entry_search.py',[],'weighted_entry_search'),
 ('budget_complexes.py',['8'],'budget_complexes'),
 ('budget_hostiles.py',[],'budget_hostiles'),
 ('resource_checks.py',[],'resource_checks'),
 ('independent_checks.py',[],'independent_checks')]

def run():
    results=[];start=time.time()
    for script,args,name in JOBS:
        for flags in [[],['-O']]:subprocess.run([sys.executable,'-X','utf8',*flags,str(HERE/script),*args],check=True)
        a,b=[HERE/(name+tag+'.json.gz') for tag in ['', '.optimized']]
        if a.read_bytes()!=b.read_bytes():raise RuntimeError('normal/optimized mismatch: '+name)
        results.append({'receipt':name,'sha256':hashlib.sha256(a.read_bytes()).hexdigest()})
        print('PAIRED',name,'elapsed',round(time.time()-start,1),flush=True)
    (HERE/'reproduction.json').write_text(json.dumps({'status':'PASS','scope':'Fourteen normal/optimized producer pairs, byte-identical deterministic receipts. Independent mathematical engines are separately identified in AUDIT.md.','pairs':results},indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':run()
