from pathlib import Path
import hashlib,json,subprocess,sys
r=Path(__file__).resolve().parent;m=json.loads((r/'MANIFEST.json').read_text());fs=m['files']
def need(c,m):
 if not c:raise RuntimeError(m)
need({str(p.relative_to(r)) for p in r.rglob('*') if p.is_file()}==set(fs)|{'MANIFEST.json'},'closed file set')
for n,h in fs.items():
 p=r/n;need(not p.is_symlink() and p.resolve().is_relative_to(r),'invalid path');need(hashlib.sha256(p.read_bytes()).hexdigest()==h,'digest '+n)
for opt,fixture in [([], 'arithmetic.normal.json'),(['-O'],'arithmetic.optimized.json')]:need(subprocess.check_output([sys.executable,'-B',*opt,str(r/'check_arithmetic.py')],cwd='/tmp')==(r/fixture).read_bytes(),'arithmetic replay')
# Frozen pre-audit map is evidence, never a hash of the eventual map containing this audit.
mp=json.loads((r/'inputs/pre-audit-map.json').read_text());d=mp['destination'];rows=mp['files'];need(len({x['destination'] for x in rows})==len(rows)==70,'snapshot map uniqueness');readme=next(x for x in rows if x['destination']==d+'/README.md');need(hashlib.sha256((r/'inputs/README.md').read_bytes()).hexdigest()==readme['sha256'],'snapshot summary')
print('PASS closed integration audit and normal/-O arithmetic; source/link integration receipts are historical workspace checks, not external runtime dependencies.')
