"""Normal/-O reproducibility without writing to any inherited packet."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def main():
    transcripts=[];receipts=[];probes=[];collars=[]
    with tempfile.TemporaryDirectory(prefix='history-phase-receipts-') as td:
        for optimized in (False,True):
            output=Path(td)/('optimized.json' if optimized else 'normal.json')
            pre=[sys.executable]+(['-O'] if optimized else [])
            commands=[pre+[str(HERE/'validate.py'),'--input',str(HERE/'new_states.json'),'--native-plan',str(HERE/'native_plan.json'),'--output',str(output)],pre+[str(HERE/'probe_marked_substitution.py')],pre+[str(HERE/'check_collar.py')]]
            text=''
            for cmd in commands:
                r=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,check=True);text+=r.stdout
                need(not r.stderr,'no unexpected stderr')
            transcripts.append(text.replace('\r\n','\n'));receipts.append(output.read_text());probes.append((HERE/'substitution_probe.json').read_text());collars.append((HERE/'collar.json').read_text())
        need(transcripts[0]==transcripts[1],'normal/-O stdout')
        need(receipts[0]==receipts[1],'normal/-O endpoint/native/phase receipt')
        need(probes[0]==probes[1],'normal/-O independent substitution witness')
        need(collars[0]==collars[1],'normal/-O inner-collar and all-root control')
    (HERE/'validation_new.json').write_text(receipts[0],encoding='utf-8',newline='\n')
    for name,text in zip(('normal.out','optimized.out'),transcripts):(HERE/name).write_text(text,encoding='utf-8',newline='\n')
    files=['HISTORY.md','REVIEW.md','endpoint_backward.cpp','validate.py','probe_marked_substitution.py','check_collar.py','reproduce.py','native_plan.json','new_states.json','validation_new.json','substitution_probe.json','collar.json','normal.out','optimized.out']
    manifest=dict(status='FINITE_EXACT_NORMAL_OPTIMIZED_MATCH',normal_optimized_identical=True,files={name:sha(HERE/name) for name in files},reproduction='python 05-knowledge/results/zero_ready_strong_switch_20261103/history_phase/reproduce.py',scope='Independent native guards, backward endpoints, named current codecs, old81/27 controls and small marked substitution counterexample; no all-value coverage theorem.')
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: independent new interfaces and marked-clone/C3 obstruction; normal and -O receipts identical')
if __name__=='__main__':main()
