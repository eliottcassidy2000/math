"""Pinned exact replay. Run normally and with -O; checks do not use assert."""
from pathlib import Path
import hashlib,json,sys,subprocess,platform
import sympy,numpy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'joint_port_geometry_20261103'
SCRIPTS=('algebra_exact','radical','semisimple','multiplication','radical_clocks','clock_split',
         'ready701','two_ready701','escape_probe','escape_birth','ternary_escapes',
         'target_match','native_return','hilbert_collatz','hostiles')
DIRS=('joint_port_geometry_20261103','ready_return_pentagon_20261103',
      'two_entry_funnels_20261103','polarized_path_calculus_20261103','defect_offset_dynamics_20261103')
EXTRA=(
 'joint_port_geometry_20261103/persistent_port.json',
 'joint_port_geometry_20261103/persistent_certificate.json',
 'joint_port_geometry_20261103/return_native.json',
 'joint_port_geometry_20261103/RETURN_233.md',
 'joint_port_geometry_20261103/JOINT_EFFECTS.md',
 'joint_port_geometry_20261103/JACOBIAN_PROBES.md',
 'two_entry_funnels_20261103/READY.md',
 'two_entry_funnels_20261103/TWO_ENTRY.md',
 'polarized_path_calculus_20261103/FUNNELS.md',
 'ready_return_pentagon_20261103/OBSTRUCTION.md',
 'regular_parent_symmetry_20261029/verify.py',
 'cayley_dickson_integral_20261030/PROOFS.md',
 'clocks_carriers_genus_20261102/COLLATZ_GEOMETRY.md',
 'joint_incidence_clocked_frames_20261102/COLLATZ_LATTICE.md',
 'collatz_hexagonal_growth_20261102/GEOMETRY.md',
 'session167_recurring_structures_20261103/honeycomb/NOTE.md',
 'lenses_tm_ca_fractran_tournaments_20261101.md',
 'coverage_machine_induction_20261103/README.md',
 'cursor_cm_phase_20261028/TRACE_CARRIER.md',
 'cursor_cm_phase_20261028/RAMANUJAN.md',
 'overnight_research_20261029/PERIODIC_FROBENIUS_LIFT.md',
)

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf8',newline='\n')

def main():
    if '--pin' in sys.argv:
        prior=json.loads((PREV/'sources.json').read_text());paths={d['path'] for d in prior['dependencies']}
        for dirname in DIRS:
            paths.update(p.relative_to(ROOT).as_posix() for p in (HERE.parent/dirname).glob('*.py'))
        paths.update('05-knowledge/results/'+p for p in EXTRA)
        save('sources.json',{'status':'PINNED_INHERITANCE_WITH_SELECTED_FRESH_REPLAY','startup_commit':'031b85ab30df',
          'hash_scheme':'SHA256 with LF newline normalization',
          'dependencies':[{'path':p,'sha256':digest(ROOT/p)} for p in sorted(paths)],
          'hostile81_27':{'status':'FRESH_REPLAY','receipt':'hostiles.json'},
          'one_birth55_orientations':prior['one_birth55_orientations'],
          'literature':[{'claim':'sharp strong-order Hamilton-path lower bound','source':'Arthur H. Busch, 2006, Corollary 1',
            'url':'https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf'},
            {'claim':'classical recursive space-filling curve','source':'David Hilbert, 1891, Mathematische Annalen 38, 459-460',
            'url':'https://www.digizeitschriften.de/download/pdf/235181684_0038/log31.pdf'}]})
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
