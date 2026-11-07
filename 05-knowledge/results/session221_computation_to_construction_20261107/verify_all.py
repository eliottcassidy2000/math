#!/usr/bin/env python3
"""Portable full checkpoint replay; all generated receipts live in a tempdir."""
import hashlib, importlib.util, json, pathlib, subprocess, sys, tempfile
ROOT=pathlib.Path(__file__).resolve().parent

def require(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    with tempfile.TemporaryDirectory(prefix='session221-') as tmp:
        for name,w,b in [('span2_B99',2,99),('span2_B199',2,199),('span3_B25',3,25)]:
            out=pathlib.Path(tmp)/(name+'.json')
            subprocess.run([sys.executable,'-B',str(ROOT/'compiler/verify.py'),'--span',str(w),'--bound',str(b),'--out',str(out)],check=True,stdout=subprocess.DEVNULL)
            require(json.loads(out.read_text())==json.loads((ROOT/'compiler'/(name+'.json')).read_text()),'compiler replay '+name)
        sys.path.insert(0,str(ROOT/'audit'))
        import audit_cap
        old=json.loads((ROOT/'audit/cap_receipt.json').read_text())
        for expected in old:
            actual=audit_cap.bfs(expected['span'],expected['cap'])
            want={k:v for k,v in expected.items() if k!='seconds'}
            require(json.loads(json.dumps(actual))==want,'independent cap BFS')
        physical=subprocess.run([sys.executable,'-B',str(ROOT/'physical/verify.py')],check=True,capture_output=True,text=True)
        require(json.loads(physical.stdout)==json.loads((ROOT/'physical/VERIFY.json').read_text()),'physical construction replay')
        for script,receipt in [('audit_physical.py','physical_receipt.json'),('audit_strong.py','strong_receipt.json')]:
            subprocess.run([sys.executable,'-B',str(ROOT/'audit'/script)],cwd=tmp,check=True,stdout=subprocess.DEVNULL)
            require(json.loads((pathlib.Path(tmp)/receipt).read_text())==json.loads((ROOT/'audit'/receipt).read_text()),'independent '+script)
    print('PASS: three compiler fixed points, three independent cap audits, physical construction and two independent physical audits.')
if __name__=='__main__':main()
