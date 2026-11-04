"""Default: independent controls. --full also replays all discovery menus."""
from pathlib import Path
from hashlib import sha256
import argparse,json,subprocess,sys,tempfile,shutil,os
HERE=Path(__file__).resolve().parent
def need(p,msg):
    if not p:raise ValueError(msg)
def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def discovery_receipt():
    names=['one_edge_search.json','two_edge_search.json','edge_birth.json','edge_two_births.json','two_edges_birth.json','three_births.json','two_births_degree2.json','two_births_degree3.json','two_births_all.json']
    return dict(status='PASS_FULL_DISCOVERY_REPLAY',cpp_sha256={p.name:digest(p) for p in sorted(HERE.glob('*.cpp'))},receipts_sha256={name:digest(HERE/name) for name in names},input_sha256=digest(HERE.parent/'return18979_palindrome_20261103/returns.json'),uniform_input_sha256=digest(HERE.parent/'return56941_two_feedback_20261103/feedback_ports.json'))
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true');args=parser.parse_args()
    compiler=shutil.which('g++');need(compiler is not None,'C++20 compiler required')
    with tempfile.TemporaryDirectory(prefix='triangle-return-') as tmp:
        bin=Path(tmp);built={}
        def compile(name,src=None):
            if name not in built:
                target=bin/(name+'.exe');subprocess.run([compiler,'-O3','-std=c++20',str(src or HERE/(name+'.cpp')),'-o',str(target)],check=True);built[name]=target
            return built[name]
        for name in ('hp','inspect'):compile(name)
        compile('endpoint',HERE.parent/'return56941_two_feedback_20261103/endpoint_checked.cpp')
        def inp(T,target=None):return ' '.join(map(str,[len(T)]+([] if target is None else [target])+list(T)))+'\n'
        def run(name,T,target,extra=()):
            return json.loads(subprocess.run([str(compile(name)),*extra],input=inp(T,target),text=True,capture_output=True,check=True).stdout)
        def compare(file,value):need(json.loads((HERE/file).read_text())==value,'replayed '+file)
        inherited=json.loads((HERE.parent/'return18979_palindrome_20261103/returns.json').read_text());T=next(z['next_rows'] for z in inherited if z['next_H']==56935)
        if args.full:
            cases=[('source56935',T,170807)]
            for z in json.loads((HERE.parent/'return56941_two_feedback_20261103/feedback_ports.json').read_text())['states']:
                if z['stage']=='second_port' and 2<=z['m']<=4:
                    cases.extend(('U'+str(z['m'])+'_'+str(d),z['rows'],3*z['H']+d) for d in (2,4))
            out=[]
            for tag,U,target in cases:
                row=run('search',U,target);row['tag']=tag;out.append(row)
            compare('one_edge_search.json',out)
            for name,file in [('two_edges','two_edge_search.json'),('edge_birth','edge_birth.json'),('edge_two_births','edge_two_births.json'),('two_edges_birth','two_edges_birth.json'),('three_births','three_births.json')]:
                compare(file,run(name,T,170807));print('PASS discovery',name,flush=True)
            for d in (2,3,-1):compare('two_births_all.json' if d==-1 else 'two_births_degree'+str(d)+'.json',run('two_births',T,170807,[str(d)]))
        # Nearby actual plus4 branch is a positive control for pruning/search.
        positive=run('two_births',T,170809,['-1']);need(any(z['strong'] for z in positive['hits']),'positive two-birth control')
        (HERE/'positive_birth_control.json').write_text(json.dumps(positive,indent=2)+'\n')
        for name,argv,outputs in [('realize.py',[],['interface.json']),('audit.py',[tmp],['audit.json','hostile_pair_deletions.json','zero_ready.json']),('clone_probe.py',[str(bin/'hp.exe')],['clone_probe.json']),('perfect_completion.py',[],['perfect_completion.json'])]:
            for optimized in (False,True):
                subprocess.run([sys.executable,*(['-O'] if optimized else []),'-B','-X','utf8',str(HERE/name),*argv],check=True,stdout=subprocess.DEVNULL)
                hashes={f:digest(HERE/f) for f in outputs}
                if not optimized:normal=hashes
                else:need(hashes==normal,'normal/-O receipts '+name)
            print('PASS independent controls',name,flush=True)
    discovery=discovery_receipt()
    if args.full:(HERE/'discovery_replay.json').write_text(json.dumps(discovery,indent=2)+'\n')
    else:need(json.loads((HERE/'discovery_replay.json').read_text())==discovery,'previous complete discovery still pinned to current C++ and input sources')
    pins=['return18979_palindrome_20261103/returns.json','return56941_two_feedback_20261103/feedback_ports.json','return56941_two_feedback_20261103/feedback_transport.json','return6325_completion_20261103/halving_response.json','return6325_completion_20261103/transport.json','return6325_completion_20261103/audit.py','return6325_completion_20261103/transport.py','return56941_two_feedback_20261103/audit_moving.py','return18979_palindrome_20261103/phases_bridge.py']
    report=dict(status='PASS',full_discovery_replayed=args.full,compiler=Path(compiler).name,python=sys.version.split()[0],source_pins={p:digest(HERE.parent/p) for p in pins},owned_sources={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.cpp')},scope='Same-agent independent algorithms; mathematical proofs and finite universes are in RETURN,TRIANGLE,PERFECT and AUDIT.')
    (HERE/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS reproduction')
if __name__=='__main__':main()
