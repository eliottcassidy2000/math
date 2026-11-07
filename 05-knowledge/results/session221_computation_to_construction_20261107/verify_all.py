#!/usr/bin/env python3
"""Portable full replay; generated audit receipts live in a temporary directory."""
import hashlib, json, pathlib, shutil, subprocess, sys, tempfile
ROOT=pathlib.Path(__file__).resolve().parent

def require(ok,msg):
    if not ok:raise RuntimeError(msg)

def verify_hashes():
    for manifest in ROOT.rglob('*SHA256SUMS'):
        for line in manifest.read_text().splitlines():
            digest,name=line.split('  ',1)
            require(hashlib.sha256((manifest.parent/name).read_bytes()).hexdigest()==digest,'hash '+str(manifest.parent/name))

def run_json(command,cwd=None):
    p=subprocess.run(command,cwd=cwd,check=True,capture_output=True,text=True)
    return json.loads(p.stdout)

def main():
    verify_hashes()
    with tempfile.TemporaryDirectory(prefix='session221-') as tmp:
        temp=pathlib.Path(tmp)
        for name,w,b in [('span2_B99',2,99),('span2_B199',2,199),('span3_B25',3,25)]:
            out=temp/(name+'.json')
            subprocess.run([sys.executable,'-B',str(ROOT/'compiler/verify.py'),'--span',str(w),'--bound',str(b),'--out',str(out)],check=True,stdout=subprocess.DEVNULL)
            require(json.loads(out.read_text())==json.loads((ROOT/'compiler'/(name+'.json')).read_text()),'compiler replay '+name)
        sys.path.insert(0,str(ROOT/'audit'))
        import audit_cap
        pair=([20,1,2,7,14],[20,9,2,5,14])
        require([audit_cap.independent_h(g) for g in pair]==[9,13],'global optimality hostile counts')
        require(sum((a^b).bit_count() for a,b in zip(*pair))==2,'global hostile one arc')
        for g in pair:
            for start in range(5):
                seen=1<<start
                while True:
                    nxt=seen
                    for v in range(5):
                        if seen>>v&1:nxt|=g[v]
                    if nxt==seen:break
                    seen=nxt
                require(seen==31,'global hostile strongness')
        for expected in json.loads((ROOT/'audit/cap_receipt.json').read_text()):
            actual=audit_cap.bfs(expected['span'],expected['cap'])
            want={k:v for k,v in expected.items() if k!='seconds'}
            require(json.loads(json.dumps(actual))==want,'independent cap BFS')
        for folder,script in [('physical','verify.py'),('automata','verify.py'),('gates','verify_gate.py')]:
            expected=json.loads((ROOT/folder/'VERIFY.json').read_text())
            for flags in ([],['-O']):
                actual=run_json([sys.executable,'-B']+flags+[str(ROOT/folder/script)])
                require(actual==expected,'candidate '+folder+' '+str(flags))
            require(expected==json.loads((ROOT/folder/'VERIFY_OPT.json').read_text()),'stored optimized '+folder)
        audit=temp/'audit';audit.mkdir()
        for source in (ROOT/'audit').glob('*.py'):shutil.copyfile(source,audit/source.name)
        for script,receipt,args in [
            ('audit_physical.py','physical_receipt.json',[]),
            ('audit_strong.py','strong_receipt.json',[]),
            ('audit_injection.py','injection_receipt.json',[]),
            ('audit_fixed_core.py','fixed_core_receipt.json',[str(ROOT/'automata')]),
            ('audit_gates.py','gates_receipt.json',[])]:
            subprocess.run([sys.executable,'-B',str(audit/script)]+args,cwd=audit,check=True,stdout=subprocess.DEVNULL)
            require(json.loads((audit/receipt).read_text())==json.loads((ROOT/'audit'/receipt).read_text()),'independent '+script)
    print('PASS: all hashes; three compiler fixed points; three independent cap BFS runs; three physical/inverse/macro verifiers in normal and optimized mode; five independent audit scripts.')
if __name__=='__main__':main()
