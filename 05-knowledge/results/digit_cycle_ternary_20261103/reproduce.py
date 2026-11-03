"""Pinned exact replay; normal and optimized runs use explicit checks."""
from pathlib import Path
import hashlib,json,sys,subprocess,platform
import sympy,numpy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'return701_algebra_20261103'
SCRIPTS=('ready2107','target2107','return2107','digits','clock','even_towers',
         'hostile_embedding','ternary_family','attachments','halving_normalization')
EXTRA=(
 'return701_algebra_20261103/native_return.json',
 'return701_algebra_20261103/radical_clocks.json',
 'return701_algebra_20261103/algebra_exact.json',
 'return701_algebra_20261103/ALGEBRA95.md',
 'return701_algebra_20261103/RETURN701.md',
 'return701_algebra_20261103/SPACE_FILLING.md',
 'return701_algebra_20261103/HANDOFF.md',
 'session143_magic_square_cyclic_bridges_20261103/magic_cycle/THEOREM.md',
 'two_edge_transport_20261102/ARITHMETIC_AND_SOURCE_AUDIT.md',
 'decoder_halving_20260925.md',
)

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf8',newline='\n')

def main():
    if '--pin' in sys.argv:
        prior=json.loads((PREV/'sources.json').read_text())
        paths={d['path'] for d in prior['dependencies']}
        paths.update(p.relative_to(ROOT).as_posix() for p in PREV.glob('*.py'))
        paths.update('05-knowledge/results/'+p for p in EXTRA)
        save('sources.json',{'status':'PINNED_INHERITANCE_WITH_SELECTED_FRESH_REPLAY',
          'startup_commit':'7fd9bc5df7e8','hash_scheme':'SHA256 with LF newline normalization',
          'dependencies':[{'path':p,'sha256':digest(ROOT/p)} for p in sorted(paths)],
          'hostile81_27':{'status':'FRESH_REPLAY','receipts':['hostile_embedding.json','attachments.json']},
          'one_birth55_orientations':prior['one_birth55_orientations'],
          'literature':[{'claim':'standard modular Jordan-block background only; application proofs supplied locally',
            'source':'Dave Benson, Seattle Summer School 2012, Modules of constant Jordan type and vector bundles on projective space, section 1',
            'url':'https://sites.math.washington.edu/~pischool/Notes/Benson-Notes.pdf'}]})
        save('manifest.json',{'scripts':[{'file':p.name,'sha256':digest(p)} for p in sorted(HERE.glob('*.py'))],
          'receipts':[{'file':k+'.json','sha256':digest(HERE/(k+'.json'))} for k in SCRIPTS]})
        print('Pinned',len(paths),'inherited files and',len(SCRIPTS),'script/receipt pairs');return
    sources=json.loads((HERE/'sources.json').read_text());manifest=json.loads((HERE/'manifest.json').read_text())
    for row in sources['dependencies']:need(digest(ROOT/row['path'])==row['sha256'],'changed inherited input: '+row['path'])
    for row in manifest['scripts']:need(digest(HERE/row['file'])==row['sha256'],'changed script: '+row['file'])
    expected={r['file']:r['sha256'] for r in manifest['receipts']};out=[]
    for k in SCRIPTS:
        cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if sys.flags.optimize else [])+[str(HERE/(k+'.py'))]
        r=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
        need(r.returncode==0,k+' failed:\n'+r.stdout+r.stderr)
        need(digest(HERE/(k+'.json'))==expected[k+'.json'],'receipt drift: '+k)
        out.append({'script':k,'stdout':r.stdout.strip(),'sha256':expected[k+'.json']})
        print('PASS',k,flush=True)
    mode='optimized' if sys.flags.optimize else 'normal'
    save('verification.'+mode+'.json',{'status':'PASS','mode':mode,'python':platform.python_version(),
      'sympy':sympy.__version__,'numpy':numpy.__version__,'inherited_files':len(sources['dependencies']),'receipts':out})
    print('PASS',mode,'all source pins and deterministic exact receipts',flush=True)
if __name__=='__main__':main()
