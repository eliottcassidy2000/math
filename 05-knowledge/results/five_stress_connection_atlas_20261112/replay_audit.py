"""Replay exact checks with assertions disabled; preserve the larger censuses.

The normal-mode receipts must already exist. This is an execution audit,
not an independent mathematical implementation. Independent controls are
inside the individual checkers and described in AUDIT.md.
"""
from pathlib import Path
import subprocess,sys,json,gzip,hashlib,time
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
JOBS=[
 ('verify_stress.py',[],['stress_receipt.json']),
 ('c3_pair_family.py',[],['c3_pair_family.json']),
 ('carousel_clocks.py',[],['carousel_clocks.json','carousel_clock_certificates.json.gz']),
 ('seed_probes.py',[],['seed_probes.json']),
 ('close_heads.py',[],['heads.json']),
 ('classify_seed_census.py',[],['seed_classification.json']),
 ('paper_transfers.py',[],['paper_transfer_receipt.json']),
 ('cyclic_port_probe.py',[],['cyclic_port_receipt.json']),
 ('joint_port_extension.py',[],['joint_port_extension.json.gz','joint_port_forms.json']),
 ('verify_tournaments.py',['8'],['tournament_receipt.json']),
 ('extend_round.py',['14'],['round_extension.json']),
]
def read(p):
    if p.suffix=='.gz':return json.loads(gzip.decompress(p.read_bytes()))
    return json.loads(p.read_text(encoding='utf8'))
def need(x,msg):
    if not x:raise RuntimeError(msg)
def run():
    results=[]
    for script,args,outputs in JOBS:
        start=time.time()
        p=subprocess.run([sys.executable,'-X','utf8','-O',str(HERE/script),*args],cwd=ROOT,capture_output=True,text=True,encoding='utf8')
        need(p.returncode==0,script+' failed: '+p.stdout[-1500:]+p.stderr[-3000:])
        controls=[]
        for name in outputs:
            baseline=HERE/name
            optimized=HERE/name.replace('.json','.optimized.json')
            a,b=read(baseline),read(optimized)
            if script=='verify_tournaments.py':
                need(b['levels']==a['levels'][:8],'Python census prefix')
                need(b['seed_records']==[q for q in a['seed_records'] if q['n']<=8],'Python marked seeds prefix')
                need(b['full_labelled_controls']==a['full_labelled_controls'],'labelled controls')
                mode='All order1..8 levels, marked records and labelled controls match the saved order14 receipt.'
            elif script=='extend_round.py':
                need(b['levels']==a['levels'][:14],'compiled census prefix')
                need(not b['additional_seed_records'],'no extension records below original14')
                mode='All order1..14 levels match the saved order18 receipt; order15..18 not rerun under -O.'
            else:
                need(a==b,script+' receipt disagreement '+name)
                mode='Complete parsed receipt equality.'
            controls.append({'file':name,'mode':mode,'normal_sha256':hashlib.sha256(baseline.read_bytes()).hexdigest(),'optimized_sha256':hashlib.sha256(optimized.read_bytes()).hexdigest()})
            optimized.unlink()
        row={'script':script,'arguments':args,'seconds':round(time.time()-start,3),'checks':controls}
        results.append(row)
        print(json.dumps(row),flush=True)
    out={'scope':'Optimized-mode replay against previously produced normal-mode exact receipts. This does not replace independent controls. The expensive complete order18 census was run once; its order1..14 prefix was independently recomputed in Python and replayed here.','python':sys.version,'results':results}
    (HERE/'replay_receipt.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
if __name__=='__main__':run()
