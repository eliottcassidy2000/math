"""Normal/-O replay with pinned inherited inputs and identical result receipts."""
from pathlib import Path
import hashlib,json,subprocess,sys,time
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
    pins=json.loads((HERE/'source_pins.json').read_text())
    for path,expected in pins['sha256_lf'].items():
        if digest(ROOT/path)!=expected:raise RuntimeError(('changed inherited input',path))
    pairs=(('verify.py','verification.json'),('branch_probe.py','branches.json'),
           ('charges.py','charges.json'),('ear_bridge.py','ear_bridge.json'),
           ('zero_cycle_prepare.py','zero_cycles.json'))
    records=[]
    for script,result in pairs:
        hashes=[];runs=[]
        for optimized in (False,True):
            command=[sys.executable,'-X','utf8','-B']+(['-O'] if optimized else [])+[str(HERE/script)]
            start=time.monotonic();run=subprocess.run(command,cwd=HERE,capture_output=True,text=True,encoding='utf8')
            if run.returncode:raise RuntimeError((command,run.stdout,run.stderr))
            hashes.append(digest(HERE/result))
            runs.append({'optimized':optimized,'command':command,'seconds':round(time.monotonic()-start,3),
                         'stdout':run.stdout.strip(),'stderr':run.stderr.strip()})
        if hashes[0]!=hashes[1]:raise RuntimeError(('normal/-O mismatch',script))
        records.append({'script':script,'result':result,'sha256_lf':hashes[0],'runs':runs})
        print(script+': PASS normal/-O, identical result',flush=True)
    receipt={'status':'PASS','pinned_inputs':len(pins['sha256_lf']),'checks_survive_optimization':True,'runs':records}
    (HERE/'reproduction.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
if __name__=='__main__':main()
