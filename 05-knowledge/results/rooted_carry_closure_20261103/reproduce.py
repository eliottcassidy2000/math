"""Source-pinned normal/optimized replay of the rooted-carry research packet."""
from pathlib import Path
import hashlib,json,sys,subprocess,platform
import sympy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SCRIPTS=('observer','carrier6','returns','field_clock','ternary_probe','golden_transport','phase_recurrence')
RELATIVE=(
 'ordered_response_ten_20261103/NEUTRAL_CARRIES.md',
 'ordered_response_ten_20261103/neutral_carries.json',
 'apex_wedge_constructor_20261103/mechanisms.py',
 'apex_wedge_constructor_20261103/PROFILE.md',
 'coverage_machine_induction_20261103/gluing/check.py',
 'coverage_machine_induction_20261103/gluing/NOTE.md',
 'coverage_machine_induction_20261103/junction.py',
 'coverage_machine_induction_20261103/RESET_JUNCTION.md',
 'source_selection_closure_20261103/native_port_selector.py',
 'source_selection_closure_20261103/NATIVE_SELECTION.md',
 'source_selection_closure_20261103/graded_carrier.py',
 'source_selection_closure_20261103/GRADED_OBSERVER.md',
 'source_selection_closure_20261103/BALANCED_BOUNDARY.md',
 'adaptive_state_induction_20261103/live_port.py',
 'adaptive_state_induction_20261103/live_port.json',
 'adaptive_metagraph_reformulation_20261103/ADAPTIVE.md',
 'two_port_return_induction_20261103/ROOTED_RESETS.md',
 'golden_inverse_defects_20261103/inverse.py',
 'golden_inverse_defects_20261103/TRANSPORT.md',
 'golden_seam_recurrence_20261103/atom_atlas.json',
 'mixed_returns_20261103/controllers.json',
 'session156_spectrum_20261103/returns/NOTE.md',
 'two_entry_funnels_20261103/READY.md',
 'two_entry_funnels_20261103/TRANSPORT.md',
 'two_entry_funnels_20261103/carrier.py',
 'two_entry_funnels_20261103/AUDIT.md',
 'twin_pythagorean_124_20261026/finite_field/README.md',
 'session55_parker_fields_integer_lifts_20261025/finite-fields/README.md',
 'golden_prime_clocks_20261003.md',
 'exception_pathcover_transfers_20261021/session24-pathcovers/PATHCOVERS.md',
 'aliquot_table_untouchable_primes_chains_seven_twentyone_20261026.md',
 'session163_ternary_triples_20261103/strong_returns/NOTE.md',
 'session163_ternary_triples_20261103/nonuniform/MARKED_ESCAPE.md',
 'session162_mrb_square_20261103/square_transfer/NOTE.md',
)
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def save(name,data):(HERE/name).write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n')

def main():
    if '--pin' in sys.argv:
        deps=[{'path':'05-knowledge/results/'+p,'sha256':digest(ROOT/'05-knowledge/results'/p)} for p in RELATIVE]
        controls=json.loads((ROOT/'05-knowledge/results/mixed_returns_20261103/controllers.json').read_text())
        save('sources.json',{'status':'INHERITED_PINNED_NOT_ALL_RETESTED','startup_commit':'d79221bc31b8',
            'incoming_through':'5891f002e2','ready_source_commit':'4ed40bd587','recovery':'RECOVERY.md',
            'hash_scheme':'SHA256 UTF8 text with LF newlines','dependencies':deps,
            'hostile81_27':controls['hostile81_27'],'one_birth55_orientations':controls['one_birth55_orientations']})
        save('manifest.json',{'receipts':[{'file':name+'.json','sha256':digest(HERE/(name+'.json'))} for name in SCRIPTS]})
        print('Pinned',len(deps),'dependencies and',len(SCRIPTS),'receipts');return
    for entry in json.loads((HERE/'sources.json').read_text())['dependencies']:
        need(digest(ROOT/entry['path'])==entry['sha256'],'changed inherited input: '+entry['path'])
    expected={r['file']:r['sha256'] for r in json.loads((HERE/'manifest.json').read_text())['receipts']}
    outputs=[]
    for name in SCRIPTS:
        cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if sys.flags.optimize else [])+[str(HERE/(name+'.py'))]
        run=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
        need(run.returncode==0,name+' failed:\n'+run.stdout+run.stderr)
        need(digest(HERE/(name+'.json'))==expected[name+'.json'],'receipt drift: '+name)
        outputs.append({'script':name,'stdout':run.stdout.strip(),'receipt_sha256':expected[name+'.json']})
        print(run.stdout.strip(),flush=True)
    mode='optimized' if sys.flags.optimize else 'normal'
    save('verification.'+mode+'.json',{'status':'PASS','mode':mode,'python':platform.python_version(),'sympy':sympy.__version__,
        'dependency_files':len(RELATIVE),'receipts':outputs})
    print('PASS:',mode,'source pins and all seven reproducible receipts')
if __name__=='__main__':main()
