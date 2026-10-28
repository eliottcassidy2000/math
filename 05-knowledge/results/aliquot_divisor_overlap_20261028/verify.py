"""Verify the sealed packet; optional complete normal/optimized replay."""
from pathlib import Path
from hashlib import sha256
import argparse,json,re,subprocess,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def need(b,s):
    if not b:raise RuntimeError(s)
def digest(b):return sha256(b.replace(b'\r\n',b'\n')).hexdigest()
def data(name):return json.loads((HERE/name).read_text(encoding='utf-8-sig'))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replay',action='store_true');args=ap.parse_args()
    names=[]
    for row in (HERE/'SHA256SUMS').read_text().splitlines():
        h,n=row.split('  ',1);need(digest((HERE/n).read_bytes())==h,'hash '+n);names.append(n)
    need(set(names)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='SHA256SUMS'},'sealed file set')
    for src in data('SOURCES.json')['repository_sources']:
        r=subprocess.run(['git','cat-file','blob',src['git_blob']],cwd=ROOT,capture_output=True,check=True)
        need(digest(r.stdout)==src['sha256_lf'],'source pin '+src['path'])
    for name in ['README.md','THEORY.md','AUDIT.md']:
        for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text(encoding='utf-8')):
            if '://' not in target and not target.startswith('#'):
                need((HERE/target.split('#')[0]).resolve().exists(),'local link '+target)
    need(data('normal.json')==data('optimized.json'),'producer normal/O')
    need(data('audit-normal.json')==data('audit-optimized.json'),'audit normal/O')
    if args.replay:
        for script,receipt in [('check.py','normal.json'),('audit.py','audit-normal.json')]:
            for flags in [[],['-O']]:
                r=subprocess.run([sys.executable,'-B',*flags,str(HERE/script)],cwd=ROOT,capture_output=True,check=True,text=True,encoding='utf-8',timeout=60)
                need(json.loads(r.stdout)==data(receipt),'replay '+script+' '+str(flags))
    print(json.dumps({'sealed_files':len(names),'source_pins':'PASS','links':'PASS','receipts':'PASS','replay':args.replay},sort_keys=True,indent=2))

if __name__=='__main__':main()
