"""No-reference certified completion, followed by separate exact controls.

Uses inherited Session277 modular engine and Session270 full counter. Timing
is explicitly not a universal or full-counter speedup claim.
"""
import argparse,json,sys,time,subprocess,tempfile,hashlib
from pathlib import Path
from residual_compiler import produce,certify
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
MOD=RESULTS/'session277_absolute_recovery_20261112/engine/modular_counter.cpp'
EXACT=RESULTS/'session270_chanpak_enumeration_20261111/enumeration/counter.cpp'
CACHE=RESULTS/'session270_chanpak_enumeration_20261111/enumeration/prefix.json'
def run(cases):
    rows=[]
    with tempfile.TemporaryDirectory() as tmp:
        mod=Path(tmp)/'mod';full=Path(tmp)/'full'
        for source,exe in [(MOD,mod),(EXACT,full)]:subprocess.run(['g++','-O3','-std=c++17',str(source),'-lgmpxx','-lgmp','-o',str(exe)],check=True)
        for n,top in cases:
            t=time.perf_counter();packet,stats=produce(n,tuple(range(3,top+1,2)));production=time.perf_counter()-t
            t=time.perf_counter();cert=certify(packet)['ground_with_R'];consumption=time.perf_counter()-t
            lo,hi=int(cert['low']),int(cert['high']);b=(hi-lo).bit_length()
            q=subprocess.run([str(mod),str(n),str(b)],capture_output=True,text=True,check=True,timeout=180)
            r=int(q.stdout);answer=lo+(r-lo)%(1<<b)
            if answer>hi:raise AssertionError('inconsistent interval/residue')
            # Only now read and compute independent exact answers.
            baseline=subprocess.run([str(full),str(n),'9'],capture_output=True,text=True,check=True,timeout=180)
            ref=int(baseline.stdout)
            cache=int(json.loads(CACHE.read_text())['values'][n])
            if answer!=ref or answer!=cache:raise AssertionError('exact comparison')
            raw=json.dumps(packet,sort_keys=True,separators=(',',':')).encode()
            rows.append({'n':n,'backbone_top':top,**stats,'production_seconds':production,'consumer_seconds':consumption,
              'residue_bits':b,'exact_answer_bits':answer.bit_length(),'residue':str(r),'interval_low':str(lo),'interval_high':str(hi),
              'recovered_integer':str(answer),'packet_sha256':hashlib.sha256(raw).hexdigest(),
              'modular_stats':json.loads(q.stderr) if q.stderr else {'seconds':0},'full_counter_stats':json.loads(baseline.stderr),
              'equals_independent_full_counter_and_cache':True})
            print(f'completed n={n}, B={top}, bits={b}',file=sys.stderr,flush=True)
    return {'status':'PASS','sources':{str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [MOD,EXACT,CACHE]},'rows':rows}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--cases',default='64:11,100:19,150:27,150:35,200:49');args=ap.parse_args()
    print(json.dumps(run([tuple(map(int,p.split(':'))) for p in args.cases.split(',')]),indent=2))
