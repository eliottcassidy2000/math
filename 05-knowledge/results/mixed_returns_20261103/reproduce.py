"""Pinned replay with normal/optimized exact result agreement."""
from pathlib import Path
import hashlib,json,subprocess,sys,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def main():
    pins=json.loads((HERE/'source_pins.json').read_text())
    for path,h in pins['sha256_lf'].items():
        if digest(ROOT/path)!=h:raise RuntimeError(('changed dependency',path))
    records=[]
    for script,out in [('return_census.py','return_census.json'),('verify.py','verification.json'),
                       ('carry_complete.py','carry_complete.json'),('controllers.py','controllers.json')]:
        runs=[];hashes=[]
        for optimized in (False,True):
            cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if optimized else [])+[str(HERE/script)]
            start=time.monotonic();p=subprocess.run(cmd,cwd=HERE,capture_output=True,text=True,encoding='utf8')
            if p.returncode:raise RuntimeError((script,optimized,p.stdout,p.stderr))
            hashes.append(digest(HERE/out));runs.append({'command':cmd,'seconds':round(time.monotonic()-start,3),'stdout':p.stdout.strip(),'stderr':p.stderr.strip()})
        if hashes[0]!=hashes[1]:raise RuntimeError(('optimized mismatch',script))
        records.append({'script':script,'output':out,'sha256_lf':hashes[0],'runs':runs})
        print(script+': PASS normal/-O',flush=True)
    (HERE/'reproduction.json').write_text(json.dumps({'status':'PASS','pinned_inputs':len(pins['sha256_lf']),'runs':records},indent=2)+'\n',encoding='utf8')
if __name__=='__main__':main()
