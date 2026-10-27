#!/usr/bin/env python3
"""Exercise explicit C++ universe rejection paths after normal/-O/UBSan replay."""
import argparse,hashlib,json,subprocess
from pathlib import Path

def need(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--normal',type=Path,required=True)
    p.add_argument('--optimized',type=Path,required=True)
    p.add_argument('--ubsan',type=Path,required=True)
    a=p.parse_args();result=[]
    for mode in ('normal','optimized','ubsan'):
        root=getattr(a,mode)
        receipt=json.loads((root/'REPLAY.json').read_text())
        need(receipt['status']=='PASS' and receipt['mode']==mode,'base replay not accepted')
        for component in ('boundary','boundary-audit'):
            exe=root/'bin'/(component+'-universe')
            data=(root/'logs'/(component+'-input.out')).read_bytes()
            lines=data.splitlines(keepends=True)
            need(lines[0]==b'4 8 1\n' and lines[-1].startswith(b'9 '),'unexpected universe input')
            tests=[
                ('outside_order',b'3 0 1\n',('order outside gate','bad n')),
                ('wrong_automorphism',b'4 8 2\n'+b''.join(lines[1:]),('automorphism mismatch',)),
                ('duplicated_ordinary_class',lines[0]+data,('isomorphic duplicate',)),
                ('deleted_ordinary_class',b''.join(lines[:-1]),('incomplete orbit mass',)),
                ('deleted_whole_order',b''.join(s for s in lines if not s.startswith(b'8 ')),('missing order',)),
            ]
            outcomes=[]
            for name,changed,expected in tests:
                run=subprocess.run([str(exe.resolve())],input=changed,capture_output=True,timeout=120)
                error=run.stderr.decode()
                need(run.returncode==1 and any(s in error for s in expected),'failed hostile rejection '+mode+'/'+component+'/'+name+': '+error)
                need('runtime error:' not in error and 'Sanitizer' not in error,'sanitizer error rather than controlled rejection')
                outcomes.append({'case':name,'returncode':run.returncode,'error':error.strip(),'changed_input_sha256':hashlib.sha256(changed).hexdigest()})
            result.append({'mode':mode,'component':component,'input_sha256':hashlib.sha256(data).hexdigest(),'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'rejected':outcomes})
    need(len(result)==6 and sum(len(r['rejected']) for r in result)==30,'control matrix incomplete')
    print(json.dumps({'status':'PASS','scope':'five fixed malformed/incomplete ordinary-universe inputs, two independent C++ authenticators, normal/O3-DNDEBUG/UBSan; no generalized parser or source-origin claim','rejections':30,'checks':result},sort_keys=True,indent=2))

if __name__=='__main__':main()
