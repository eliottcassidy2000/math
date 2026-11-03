"""Pinned deterministic normal/optimized replay of the joint-port packet."""
from pathlib import Path
import hashlib,json,sys,subprocess,platform
import sympy,numpy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'lambda_emergent_20261103'
SCRIPTS=('return_probe','return_native','word_return','joint_probe','rank_transfer',
         'algebra_probe','filtration','split_algebra','persistent_port',
         'persistent_certificate','literal_interfaces','recurrences',
         'geometric_phases','jacobian_observer')
EXTRA=(
 '05-knowledge/results/lambda_emergent_20261103/JOINT_CARRIER.md',
 '05-knowledge/results/lambda_emergent_20261103/PHASE_LEGALIZER.md',
 '05-knowledge/results/lambda_emergent_20261103/READY_REPAIR.md',
 '05-knowledge/results/mixed_returns_20261103/controllers.json',
 '05-knowledge/results/ready_return_pentagon_20261103/OBSTRUCTION.md',
 '05-knowledge/results/hadwiger_six_seven_eight_20261025/SYNTHESIS.md',
 '05-knowledge/results/hadwiger_six_seven_eight_20261025/metric_spindle/RESULT.md',
 '05-knowledge/results/hadwiger_six_seven_eight_20261025/fano_eight/RESULT.md',
 '05-knowledge/results/spectrum_descent_jacobian_20261102/JACOBIAN_FIBRES.md',
 '05-knowledge/results/planar_jc48_sep06_torsion_audit.md',
 '01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md',
 '05-knowledge/results/session166_planar_metagraph_20261103/reflection_voltages/REFLECTION_VOLTAGES.md',
 '05-knowledge/results/session166_planar_metagraph_20261103/unit_diamonds/NOTE.md',
 '05-knowledge/results/defect_offset_dynamics_20261103/ESCAPES.md',
 '05-knowledge/results/six_roots_hidden_bits_20261103/PORT_MEMORY.md',
)

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf8',newline='\n')

def main():
    if '--pin' in sys.argv:
        prior=json.loads((PREV/'sources.json').read_text());paths={d['path'] for d in prior['dependencies']}
        paths.update(str(p.relative_to(ROOT)).replace('\\','/') for p in PREV.glob('*.py'))
        paths.update(EXTRA)
        save('sources.json',{'status':'INHERITED_PINNED_NOT_ALL_RETESTED','startup_commit':'2af5e7549da1',
          'note':'PHASE_LEGALIZER.md includes the live-port correction; incoming integration through596ae1ea47 changes only the older ready README route among prior pins',
          'hash_scheme':'SHA256 with LF newline normalization',
          'dependencies':[{'path':p,'sha256':digest(ROOT/p)} for p in sorted(paths)],
          'hostile81_27':prior['hostile81_27'],'one_birth55_orientations':prior['one_birth55_orientations']})
        save('manifest.json',{'scripts':[{'file':k+'.py','sha256':digest(HERE/(k+'.py'))} for k in SCRIPTS],
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
