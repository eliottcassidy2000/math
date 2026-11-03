"""Pinned-input replay, two complete C++ transforms, and normal/-O agreement."""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]

def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def run(command):
    start=time.monotonic();r=subprocess.run(command,cwd=HERE,capture_output=True,text=True,encoding='utf8')
    if r.returncode:raise RuntimeError((command,r.returncode,r.stdout[-2000:],r.stderr))
    return r,round(time.monotonic()-start,3)

def main():
    pins=json.loads((HERE/'source_pins.json').read_text())
    for path,h in pins['sha256_lf'].items():
        if digest(ROOT/path)!=h:raise RuntimeError(('changed source input',path))
    compiler=shutil.which('g++') or 'C:/Users/Eliott/scoop/apps/gcc/current/bin/g++.exe'
    exe=HERE/('.reproduce-birth.exe' if sys.platform=='win32' else '.reproduce-birth')
    commands=[]
    try:
        cmd=[compiler,'-O2','-std=c++17',str(HERE/'two_birth.cpp'),'-o',str(exe)]
        r,seconds=run(cmd);commands.append({'command':cmd,'seconds':seconds,'exit_code':0})
        for source,born,out in (('parents55.txt',2,'two_birth.json'),('representatives55.txt',3,'three_birth.json')):
            cmd=[str(exe),str(HERE/source),str(born)];r,seconds=run(cmd)
            data=json.loads(r.stdout);(HERE/out).write_text(r.stdout,encoding='utf8')
            commands.append({'command':cmd,'seconds':seconds,'exit_code':0,
                             'summary':[(d['id'],d['H111_count']) for d in data]})
    finally:
        if exe.exists():exe.unlink()
    probes=(('channels.py','channels.json'),('arithmetic_probe.py','arithmetic.json'),
            ('port_guard.py','port_guard.json'),('birth_verify.py','birth_verification.json'),
            ('context_entry.py','context_entry.json'))
    records=[]
    for script,out in probes:
        hashes=[];runs=[]
        for optimized in (False,True):
            cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if optimized else [])+[str(HERE/script)]
            r,seconds=run(cmd);hashes.append(digest(HERE/out))
            runs.append({'command':cmd,'seconds':seconds,'stdout':r.stdout.strip(),'stderr':r.stderr.strip()})
        if hashes[0]!=hashes[1]:raise RuntimeError(('normal/optimized mismatch',script))
        records.append({'script':script,'result':out,'sha256_lf':hashes[0],'runs':runs})
        print(script+': PASS normal/-O with identical result',flush=True)
    result={'status':'PASS','pinned_inputs':len(pins['sha256_lf']),'cpp_runs':commands,'python_runs':records}
    (HERE/'reproduction.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
if __name__=='__main__':main()
