"""Compare the graded compiler to independent literal endpoint subset DP."""
from pathlib import Path
import sys,json,subprocess
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'polarized_path_calculus_20261103'))
from recognize import unrank
def need(b,s):
    if not b:raise ValueError(s)
def main():
    p=subprocess.run([sys.argv[1],'20','audit'],capture_output=True,text=True,check=True)
    rows=[json.loads(s) for s in p.stdout.splitlines()];old=json.loads((HERE/'one_start_bank.json').read_text())['rows']
    need([x for x in rows if x['kind']!='query']==old,'entire compiler search reproduces')
    Q=[x for x in rows if x['kind']=='query'];text=''
    for z in Q:
        T=unrank(z['n'],z['rank']);text+=' '.join(map(str,[len(T),0,1,*T]))+'\n'
    p=subprocess.run([sys.argv[2]],input=text,capture_output=True,text=True,check=True)
    for z,line in zip(Q,p.stdout.splitlines(),strict=True):
        v=list(map(int,line.split()));need(z['values']==[v[0],v[2],v[3],v[1],v[5],v[6]],('literal path-bin compiler',z))
    targets=json.loads((HERE/'targets.json').read_text());text=''.join(' '.join(map(str,[len(z['rows']),z['root'],z['inport'],*z['rows']]))+'\n' for z in targets)
    p=subprocess.run([sys.argv[2]],input=text,text=True,capture_output=True,check=True)
    for z,line in zip(targets,p.stdout.splitlines(),strict=True):
        v=list(map(int,line.split()));need([v[0],v[3],v[4],v[5]]==z['profile'],'all positive target profiles independently')
    result=dict(status='FINITE_EXACT_INDEPENDENT_COMPILER_CHECK',full_parent_classes=sum(z['classes'] for z in old if z['kind']=='census'),
        subset_queries=len(Q),query_scope='every class through9; first/middle/last rank at each order10..20',positive_targets=len(targets),
        one_ready_target_orders={str(h):sorted(z['n']+1 for z in old if z['kind']=='hit' and z['target']==h) for h in (56933,56935,56939,56941)},
        all_order_bound='Inherited complete READY classification; Hparent <= (H+1)/2 <=28471; Busch f(21)=46875')
    (HERE/'audit_bank.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS complete compiler replay and',len(Q),'independent literal queries',flush=True)
if __name__=='__main__':main()
