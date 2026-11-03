"""Rebuild independent certificates, check pins, compare modes, seal files."""
from pathlib import Path
import hashlib,json,re,subprocess,sys,tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def need(ok,message):
    if not ok:raise RuntimeError(message)
def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def run(cmd):
    p=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
    need(p.returncode==0,(cmd,p.returncode,p.stdout,p.stderr))
    return p.stdout
def main():
    pins=json.loads((HERE/'sources.json').read_text(encoding='utf8'))
    for path,h in pins['sha256_lf'].items():need(digest(ROOT/path)==h,('changed source',path))
    old_census=json.loads((HERE/'census.json').read_text(encoding='utf-8-sig'))
    old_audit=json.loads((HERE/'audit.json').read_text(encoding='utf8'))
    with tempfile.TemporaryDirectory(prefix='math-exclusion-7-21-') as directory:
        # Only this freshly allocated, verified temp directory is cleaned up.
        d=Path(directory).resolve()
        need(d.parent==Path(tempfile.gettempdir()).resolve() and d.name.startswith('math-exclusion-7-21-'),'temporary path guard')
        exe=d/'census.exe'
        parsed=[]
        for flags in (['-O2','-DNDEBUG'],['-O0']):
            run(['g++','-std=c++17',*flags,str(HERE/'census.cpp'),'-o',str(exe)])
            parsed.append(json.loads(run([str(exe)])))
        need(parsed[0]==parsed[1]==old_census,'C++ modes or frozen census differ')
    (HERE/'census.json').write_text(json.dumps(old_census,indent=2)+'\n',encoding='utf8')
    modes=[]
    for optimized in (False,True):
        cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if optimized else [])+[str(HERE/'audit.py')]
        stdout=run(cmd);new=json.loads((HERE/'audit.json').read_text(encoding='utf8'))
        need(new==old_audit,'Python modes or frozen audit differ')
        modes.append(stdout.strip())
    need(modes[0]==modes[1],'mode stdout differs')
    receipt={'status':'PASS','compiler_modes':['O2 DNDEBUG','O0'],'python_modes':['normal','optimized'],
             'pinned_inputs':len(pins['sha256_lf']),'audit_stdout':modes[0],
             'census_sha256_lf':digest(HERE/'census.json'),'audit_sha256_lf':digest(HERE/'audit.json')}
    (HERE/'reproduction.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    for f in HERE.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',f.read_text(encoding='utf8')):
            if '://' not in link:need((f.parent/link.split('#')[0]).exists(),('broken link',f.name,link))
    files={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.suffix in {'.py','.cpp','.md','.json'} and p.name!='manifest.json'}
    (HERE/'manifest.json').write_text(json.dumps({'status':'VERIFIED','normalization':'CRLF to LF','sha256_lf':files},indent=2)+'\n',encoding='utf8')
    print('PASS: two C++ modes, normal/optimized Python, frozen certificates, source pins and links')
if __name__=='__main__':main()
