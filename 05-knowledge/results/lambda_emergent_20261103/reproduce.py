"""Pinned normal/optimized replay of the emergent effect and ready-grid packet."""
from pathlib import Path
import hashlib,json,sys,subprocess,platform
import sympy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'rooted_carry_closure_20261103'
SCRIPTS=('effects','algebra','lattice','controls','ready','one_hole','legalizer','clock_cube')
EXTRA=(
 'rooted_carry_closure_20261103/ROOTED_RETURN.md',
 'rooted_carry_closure_20261103/PHASE_49_50_51.md',
 'rooted_carry_closure_20261103/SINGULAR_RETURN.md',
 'rooted_carry_closure_20261103/HANDOFF.md',
 'two_selector_boolean_20261102/CARRIER.md',
 'two_selector_boolean_20261102/INSPIRATION.md',
 'session143_magic_square_cyclic_bridges_20261103/lambda_octonions/PROOF.md',
 'session143_frobenius_magic_followup_20261103/producer/PROOF.md',
 'session164_emergent_20261103/triple_gluing/NOTE.md',
 'ready_return_pentagon_20261103/README.md',
 'ready_return_pentagon_20261103/OBSTRUCTION.md',
 'session164_emergent_20261103/profile_algebras/CANONICAL_RETURNS.md',
 'session165_geometry_representation_20261103/guard_adjunctions/NOTE.md',
 'session165_geometry_representation_20261103/slot_representation/UNIFORM_PLANCHEREL.md',
)
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf8',newline='\n')
def main():
    if '--pin' in sys.argv:
        prior=json.loads((PREV/'sources.json').read_text())
        paths={d['path'] for d in prior['dependencies']}
        paths.update(str(p.relative_to(ROOT)).replace('\\','/') for p in PREV.glob('*.py'))
        paths.update('05-knowledge/results/'+p for p in EXTRA)
        # Inherited hostiles remain explicit controls; they are not re-enumerated here.
        save('sources.json',{'status':'INHERITED_PINNED_NOT_ALL_RETESTED','startup_commit':'1700971fb579',
             'hash_scheme':'SHA256 with LF newline normalization',
             'dependencies':[{'path':p,'sha256':digest(ROOT/p)} for p in sorted(paths)],
             'hostile81_27':prior['hostile81_27'],'one_birth55_orientations':prior['one_birth55_orientations']})
        save('manifest.json',{'receipts':[{'file':k+'.json','sha256':digest(HERE/(k+'.json'))} for k in SCRIPTS]})
        print('Pinned',len(paths),'inherited files and',len(SCRIPTS),'result receipts');return
    sources=json.loads((HERE/'sources.json').read_text())
    for row in sources['dependencies']:need(digest(ROOT/row['path'])==row['sha256'],'changed inherited input: '+row['path'])
    expected={r['file']:r['sha256'] for r in json.loads((HERE/'manifest.json').read_text())['receipts']}
    out=[]
    for k in SCRIPTS:
        cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if sys.flags.optimize else [])+[str(HERE/(k+'.py'))]
        r=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
        need(r.returncode==0,k+' failed:\n'+r.stdout+r.stderr)
        need(digest(HERE/(k+'.json'))==expected[k+'.json'],'receipt drift: '+k)
        out.append({'script':k,'stdout':r.stdout.strip(),'sha256':expected[k+'.json']})
        print('PASS',k,flush=True)
    mode='optimized' if sys.flags.optimize else 'normal'
    save('verification.'+mode+'.json',{'status':'PASS','mode':mode,'python':platform.python_version(),'sympy':sympy.__version__,
         'inherited_files':len(sources['dependencies']),'receipts':out})
    print('PASS',mode,'all source pins and',len(SCRIPTS),'deterministic exact receipts')
if __name__=='__main__':main()
