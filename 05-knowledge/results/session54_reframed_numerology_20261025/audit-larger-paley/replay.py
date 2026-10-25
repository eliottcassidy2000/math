#!/usr/bin/env python3
import pathlib,json,hashlib,subprocess,sys,tempfile,shutil,copy
ROOT=pathlib.Path(__file__).resolve().parent

def require(ok,reason):
    if not ok: raise RuntimeError(reason)

def integrity(root):
    manifest=json.loads((root/'MANIFEST.json').read_text())['files']
    require({str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}==set(manifest)|{'MANIFEST.json'},'closed file set')
    for f,h in manifest.items(): require(hashlib.sha256((root/f).read_bytes()).hexdigest()==h,'digest '+f)

def main():
    integrity(ROOT)
    for flags,name in ((['-B'],'result.json'),(['-B','-O'],'result-O.json')):
        p=subprocess.run([sys.executable,*flags,str(ROOT/'audit.py')],cwd='/',capture_output=True,check=True)
        require(p.stdout==(ROOT/name).read_bytes()==(ROOT/'result.json').read_bytes(),'normal optimized identity')
    if '--self-test' in sys.argv:
        import audit
        cert=json.loads((ROOT/'inputs/session54-larger-paley-certificates/certificate-19.json').read_text())
        for kind in ('cycle','endpoint','rank','minor','seed','modulus'):
            c=copy.deepcopy(cert)
            if kind=='cycle': c['base_cycle'][1]=c['base_cycle'][0]
            if kind=='endpoint': c['endpoint_path']=list(reversed(c['endpoint_path']))
            if kind=='rank': c['target_rank']-=1
            if kind=='minor': c['pivot_columns'][1]=c['pivot_columns'][0]
            if kind=='seed': c['seed_cycles'].append(c['base_cycle'])
            if kind=='modulus': c['modulus']=100
            rejected=False
            try: audit.verify(c)
            except ValueError: rejected=True
            require(rejected,'accepted semantic hostile '+kind)
        # Standalone singular/nonzero modular matrix controls for independent routine.
        require(audit.column_determinant([[1,2],[2,4]],101)==0,'singular hostile')
        require(audit.column_determinant([[0,1],[1,0]],101)==100,'column swap control')
        with tempfile.TemporaryDirectory() as tmp:
            dst=pathlib.Path(tmp)/'copy';shutil.copytree(ROOT,dst);integrity(dst)
            for kind in ('add','change','remove'):
                f=dst/('unexpected.txt' if kind=='add' else 'result.json');old=f.read_bytes() if f.exists() else None
                if kind=='remove':f.unlink()
                else:f.write_bytes((old or b'')+b'!')
                rejected=False
                try:integrity(dst)
                except RuntimeError:rejected=True
                require(rejected,'accepted file hostile '+kind)
                if old is None:f.unlink()
                else:f.write_bytes(old)
            subprocess.run([sys.executable,'-B',str(dst/'replay.py')],cwd='/',check=True,capture_output=True)
    print('PASS: closed text-only hashes, independent normal/-O, relocated replay, semantic and file hostiles')
if __name__=='__main__':main()
