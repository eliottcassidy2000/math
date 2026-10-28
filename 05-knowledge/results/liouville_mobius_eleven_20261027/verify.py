"""Verify sealed evidence; --replay reruns all producers and independent audits."""
from pathlib import Path
from hashlib import sha256
import argparse,json,os,re,subprocess,sys,tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def need(b,s):
    if not b:raise RuntimeError(s)

def digest(b):return sha256(b.replace(b'\r\n',b'\n')).hexdigest()
def read(name):return json.loads((HERE/name).read_text(encoding='utf-8-sig'))
def run(args):return subprocess.run(args,cwd=ROOT,capture_output=True,check=True,encoding='utf-8',timeout=120).stdout

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replay',action='store_true');args=ap.parse_args()
    sealed=[]
    for line in (HERE/'SHA256SUMS').read_text().splitlines():
        h,name=line.split('  ',1);need(digest((HERE/name).read_bytes())==h,'packet hash: '+name);sealed.append(name)
    need(set(sealed)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='SHA256SUMS'},'exact sealed file set')
    sources=read('SOURCES.json')
    for s in sources['repository_sources']:
        blob=subprocess.run(['git','cat-file','blob',s['git_blob']],cwd=ROOT,capture_output=True,check=True).stdout
        need(digest(blob)==s['sha256_lf'],'pinned source: '+s['path'])
    for name in ['README.md','THEORY.md','ELEVEN.md','INCOMING.md']:
        for target in re.findall(r'\]\(([^)]+)\)',(HERE/name).read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            need((HERE/target.split('#')[0]).exists(),'local link '+target)
    need(read('normal.json')==read('optimized.json'),'normal/optimized arithmetic agreement')
    need(read('hostile-audit-normal.json')==read('hostile-audit-optimized.json'),'normal/optimized audit agreement')
    if args.replay:
        for flags in [[],['-O']]:
            got=json.loads(run([sys.executable,'-B',*flags,str(HERE/'check.py')]))
            need(got==read('normal.json'),'replayed arithmetic '+str(flags))
            got=json.loads(run([sys.executable,'-B',*flags,str(HERE/'audit_hostile.py')]))
            need(got==read('hostile-audit-normal.json'),'replayed independent audit '+str(flags))
        with tempfile.TemporaryDirectory(prefix='liouville-mobius-hostile-') as d:
            exe=Path(d)/('hostile.exe' if os.name=='nt' else 'hostile')
            run(['g++','-std=c++17','-O3',str(HERE/'hostile.cpp'),'-o',str(exe)])
            need(json.loads(run([str(exe)]))==read('hostile.json'),'complete C++ replay')
    print(json.dumps({'sealed_files':len(sealed),'pinned_sources':len(sources['repository_sources']),
                      'local_links':'PASS','receipts':'PASS','replay':args.replay},sort_keys=True,indent=2))

if __name__=='__main__':main()
