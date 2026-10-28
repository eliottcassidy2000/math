"""Validate the sealed research packet, optionally replaying all computations."""
from pathlib import Path
from hashlib import sha256
import argparse,json,re,subprocess,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def need(b,s):
    if not b:raise RuntimeError(s)
def digest(b):return sha256(b.replace(b'\r\n',b'\n')).hexdigest()
def read(n):return json.loads((HERE/n).read_text(encoding='utf-8'))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replay',action='store_true');args=ap.parse_args()
    sealed=[]
    for line in (HERE/'SHA256SUMS').read_text().splitlines():
        h,n=line.split('  ',1);need(digest((HERE/n).read_bytes())==h,'hash '+n);sealed.append(n)
    need(set(sealed)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='SHA256SUMS'},'sealed inventory')
    for src in read('SOURCES.json')['repository_sources']:
        b=subprocess.check_output(['git','cat-file','blob',src['git_blob']],cwd=ROOT)
        need(digest(b)==src['sha256_lf'],'source '+src['path'])
    for name in ('README.md','THEORY.md','RECOVERY.md','AUDIT.md'):
        for link in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text(encoding='utf-8')):
            if '://' not in link and not link.startswith('#'):
                need((HERE/link.split('#')[0]).resolve().exists(),'link '+link)
    for prefix in ('','audit-'):need(read(prefix+'normal.json')==read(prefix+'optimized.json'),'normal/O '+prefix)
    if args.replay:
        for script,name in [('check.py','normal.json'),('audit.py','audit-normal.json')]:
            for flags in ([],['-O']):
                r=subprocess.run([sys.executable,'-B',*flags,str(HERE/script)],cwd=ROOT,capture_output=True,check=True,timeout=60)
                need(json.loads(r.stdout)==read(name),'replay '+script)
    print(json.dumps({'sealed_files':len(sealed),'source_pins':'PASS','links':'PASS','receipts':'PASS','replay':args.replay},sort_keys=True))

if __name__=='__main__':main()
