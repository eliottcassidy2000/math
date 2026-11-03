"""Pinned exact replay of the return, observer and contraction certificates."""
from pathlib import Path
import hashlib,json,sys,subprocess,platform
import sympy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PREV=HERE.parent/'digit_cycle_ternary_20261103'
SCRIPTS=('ready6325','target6325','return6325','matching_carry','completion_hostile','observer_filtration','halving_response')
EXTRA=('digit_cycle_ternary_20261103/return2107.json','digit_cycle_ternary_20261103/hostile_embedding.json',
 'digit_cycle_ternary_20261103/halving_normalization.json','digit_cycle_ternary_20261103/SYMMETRY_AND_PORTS.md',
 'digit_cycle_ternary_20261103/HANDOFF.md','joint_port_geometry_20261103/persistent_port.json',
 'joint_port_geometry_20261103/JOINT_EFFECTS.md','return701_algebra_20261103/algebra_exact.json',
 'boundary_carry_returns_20261103/PROOFS.md','six_output_return_20261103/COMMON_GRAMMAR.md',
 'session170_holonomy_homotopy_20261103/guards/NOTE.md','session171_knots_information_20261103/information/NOTE.md')
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def save(name,d):(HERE/name).write_text(json.dumps(d,indent=2)+'\n',encoding='utf8',newline='\n')
def main():
    if '--pin' in sys.argv:
        old=json.loads((PREV/'sources.json').read_text());paths={d['path'] for d in old['dependencies']}
        paths.update(p.relative_to(ROOT).as_posix() for p in PREV.glob('*.py'))
        paths.update('05-knowledge/results/'+p for p in EXTRA)
        save('sources.json',{'status':'PINNED_INHERITANCE_WITH_SELECTED_FRESH_REPLAY','startup_commit':'31b13dd9a2c8',
          'hash_scheme':'SHA256 with LF newline normalization','dependencies':[{'path':p,'sha256':digest(ROOT/p)} for p in sorted(paths)],
          'hostile81_27':{'status':'FRESH_REPLAY','receipt':'completion_hostile.json'},'one_birth55_orientations':old['one_birth55_orientations']})
        files=sorted(list(HERE.glob('*.py'))+list(HERE.glob('*.cpp')))
        save('manifest.json',{'code':[{'file':p.name,'sha256':digest(p)} for p in files],
          'receipts':[{'file':k+'.json','sha256':digest(HERE/(k+'.json'))} for k in SCRIPTS]})
        print('Pinned',len(paths),'inherited files and',len(SCRIPTS),'receipts');return
    sources=json.loads((HERE/'sources.json').read_text());manifest=json.loads((HERE/'manifest.json').read_text())
    for row in sources['dependencies']:need(digest(ROOT/row['path'])==row['sha256'],'changed inherited input: '+row['path'])
    for row in manifest['code']:need(digest(HERE/row['file'])==row['sha256'],'changed code: '+row['file'])
    expected={r['file']:r['sha256'] for r in manifest['receipts']};out=[]
    for k in SCRIPTS:
        cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if sys.flags.optimize else [])+[str(HERE/(k+'.py'))]
        r=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
        need(r.returncode==0,k+' failed:\n'+r.stdout+r.stderr)
        need(digest(HERE/(k+'.json'))==expected[k+'.json'],'receipt drift: '+k)
        out.append({'script':k,'stdout':r.stdout.strip(),'sha256':expected[k+'.json']});print('PASS',k,flush=True)
    mode='optimized' if sys.flags.optimize else 'normal'
    compiler=subprocess.run(['g++','--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0]
    save('verification.'+mode+'.json',{'status':'PASS','mode':mode,'python':platform.python_version(),'sympy':sympy.__version__,
      'compiler':compiler,'inherited_files':len(sources['dependencies']),'receipts':out})
    print('PASS',mode,'all exact replay certificates',flush=True)
if __name__=='__main__':main()
