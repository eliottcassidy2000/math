"""Relocate payloads; test closure and hostile manifest mutations in both modes."""
import json,pathlib,shutil,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent
m=json.loads((ROOT/'MANIFEST.json').read_text())
results={}
for mode in ([],['-O']):
    label='optimized' if mode else 'normal'
    for case in ('clean','missing','changed','extra','duplicate','unsafe','duplicate_json_key'):
        with tempfile.TemporaryDirectory() as td:
            out=pathlib.Path(td)/'relocated'; out.mkdir()
            for f in ['MANIFEST.json']+[r['path'] for r in m['files']]:
                dest=out/f;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/f,dest)
            manifest=json.loads((out/'MANIFEST.json').read_text())
            payload=out/'sources/inherited-proof.md'
            if case=='missing': payload.unlink()
            elif case=='changed': payload.write_bytes(payload.read_bytes()+b'\nchanged\n')
            elif case=='extra': (out/'extra.bin').write_bytes(b'extra')
            elif case=='duplicate': manifest['files'].append(manifest['files'][0])
            elif case=='unsafe': manifest['files'][0]['path']='../escape'
            if case in ('duplicate','unsafe'): (out/'MANIFEST.json').write_text(json.dumps(manifest))
            if case=='duplicate_json_key':
                text=(out/'MANIFEST.json').read_text(); text=text.replace('"schema":','"schema":"bad", "schema":',1); (out/'MANIFEST.json').write_text(text)
            proc=subprocess.run([sys.executable,*mode,str(out/'verify.py')],capture_output=True,text=True,cwd=out)
            ok=proc.returncode==0
            if ok!=(case=='clean'): raise ValueError(label+' '+case+' unexpected acceptance: '+proc.stdout+proc.stderr)
            results[label+'_'+case]='PASS' if ok else 'REJECTED_AS_REQUIRED'
print(json.dumps({'status':'PASS','cases':results},sort_keys=True,indent=2))
