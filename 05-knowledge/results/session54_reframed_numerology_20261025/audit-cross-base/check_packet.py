from pathlib import Path
import tempfile,hashlib,subprocess,json,sys
P=Path(__file__).resolve().parent.parent/'session54-cross-base-fractal'
expected='92db3b170771188fac5c9d55faec8c4c5d18e4fc8f34d93c74b78abfd0794e95'
def need(ok,label):
    if not ok:raise RuntimeError(label)
hash=lambda b:hashlib.sha256(b).hexdigest()
need(hash((P/'manifest.json').read_bytes())==expected,'manifest root')
m=json.loads((P/'manifest.json').read_text())
need({str(f.relative_to(P)) for f in P.rglob('*') if f.is_file()}==set(m)|{'manifest.json'},'source closure')
for f,h in m.items():need(hash((P/f).read_bytes())==h,'digest '+f)
root=P.parent
lineage=[('sources/mod3-digit-law.md','session53-cooper-prime-dichotomy/REPORT.md'),('sources/mod9-digit-law.md','session53-cooper-mod9-extension/REPORT.md'),('sources/mod4-dyadic-law.md','session53-dyadic-mod4-supplement/REPORT.md')]
for local,prior in lineage:
    need((P/local).read_bytes()==(root/prior).read_bytes(),'Session53 exact source '+local)
    oldm=json.loads((root/prior).parent.joinpath('MANIFEST.json').read_text())
    need(oldm['files']['REPORT.md']==hash((P/local).read_bytes()),'Session53 manifest source '+local)
repo=root/'session50-publication/clean-repo'
commit='bbe2ae278aeabdf707fdfca3c1bf0d3208564b03'
for local,path in [('sources/clock-fractal.md','05-knowledge/results/mod18_mod19_seven_sixtythree_fractal_20261006.md'),('sources/prior-transversality-published.md','05-knowledge/results/collatz_procgen_20260922_transversality_foundry.md')]:
    need(subprocess.check_output(['git','-C',str(repo),'show',commit+':'+path])==(P/local).read_bytes(),'published provenance '+local)
    if 'transversality' in local:
        need(subprocess.check_output(['git','-C',str(root/'math'),'show','1b91742ce55a202b73c80d56561ba56d00b503a6:'+path])==(P/local).read_bytes(),'original HEAD provenance')
        need((root/'math'/path).read_bytes()==(P/local).read_bytes(),'original worktree provenance')
with tempfile.TemporaryDirectory() as t:
    dst=Path(t)/'closed';dst.mkdir()
    # Only the manifest and explicitly named members are copied.
    for f in list(m)+['manifest.json']:
        (dst/f).parent.mkdir(parents=True,exist_ok=True);(dst/f).write_bytes((P/f).read_bytes())
    for opt in ([],['-O']):
        for program in ('replay.py','hostile.py'):
            result=subprocess.check_output([sys.executable,*opt,str(dst/program)],cwd=t)
            print((' '.join(opt) or 'normal')+' '+program+': '+result.decode().strip())
print('PASS assigned manifest, exact source lineage, isolated manifest-only normal/-O replay and hostiles')
