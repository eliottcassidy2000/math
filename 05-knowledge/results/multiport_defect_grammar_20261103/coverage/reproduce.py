"""Reproduce owned controls only; pin separately reviewed proof sources."""
from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def main():
    output=[];receipts=[]
    for optimized in (False,True):
        pre=[sys.executable]+(['-O'] if optimized else []);text=''
        for name in ('check_coverage.py','audit_c2.py'):
            r=subprocess.run(pre+[str(HERE/name)],cwd=ROOT,text=True,capture_output=True,check=True);need(not r.stderr,'unexpected stderr');text+=r.stdout
        output.append(text.replace('\r\n','\n'));receipts.append([(HERE/name).read_text() for name in ('coverage.json','c2_audit.json')])
    need(output[0]==output[1] and receipts[0]==receipts[1],'normal/-O receipts differ')
    for name,txt in zip(('normal.out','optimized.out'),output):(HERE/name).write_text(txt,encoding='utf-8',newline='\n')
    files=['COVERAGE.md','REVIEW.md','check_coverage.py','audit_c2.py','reproduce.py','coverage.json','c2_audit.json','normal.out','optimized.out']
    refs=['01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md','05-knowledge/results/balanced_ternary_returns_20261103/conditional/NOTE.md','05-knowledge/results/two_entry_funnels_20261103/READY.md','05-knowledge/results/source_selection_closure_20261103/GRAMMAR.md','05-knowledge/results/source_selection_closure_20261103/BRANCH.md','05-knowledge/results/multiport_defect_grammar_20261103/construction/algebra.py','05-knowledge/results/multiport_defect_grammar_20261103/construction/RESULT.md','05-knowledge/results/multiport_defect_grammar_20261103/transport/PROOF.md','05-knowledge/results/multiport_defect_grammar_20261103/transport/check_transport.py','05-knowledge/results/zero_ready_strong_switch_20261103/literal/realized_generation2_seed4_depth3_all_prefixes.json']
    report=dict(status='NORMAL_OPTIMIZED_FINITE_CONTROLS_IDENTICAL_AND_SCOPED_MATHEMATICAL_REVIEW',normal_optimized_identical=True,files={name:digest(HERE/name) for name in files},reviewed_or_inherited_source_pins={name:digest(ROOT/name) for name in refs},reproduction='python 05-knowledge/results/multiport_defect_grammar_20261103/coverage/reproduce.py',scope='Arithmetic induction and explicit source-type obligations; independent small C2/C4 and literal start-defect controls. Does not claim arbitrary-depth digit closure or global tournament spectrum completeness.')
    (HERE/'manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: owned arithmetic/C2/entry controls match normal and -O; reviewed transport/construction sources pinned')
if __name__=='__main__':main()
