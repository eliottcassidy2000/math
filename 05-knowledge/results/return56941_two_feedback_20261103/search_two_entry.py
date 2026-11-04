"""All tails through parent18 over a bounded, sourced two-entry core bank."""
from search_funnels import *
from itertools import combinations
sys.path.insert(0,str(HERE.parent/'ready_return_pentagon_20261103'))
import ready_inverse as inv
def main():
    cores=[];seen=set()
    for n in range(3,11):
        for T in inv.two_start(n):
            z=(T,(0,))
            if z not in seen:seen.add(z);cores.append(dict(rows=T,ports=[0]))
        for T in inv.unique_cycle_classes(n):
            roots=[r for r in range(n) if inv.unique_start_path(T,r) is not None]
            for J in combinations(roots,2):
                z=(T,J)
                if z not in seen:seen.add(z);cores.append(dict(rows=T,ports=J))
        print('core order',n,'bank',len(cores),flush=True)
    text=''.join(' '.join(map(str,[len(c['rows']),sum(1<<j for j in c['ports']),i,*c['rows']]))+'\n' for i,c in enumerate(cores))
    exe=Path(tempfile.gettempdir())/'two-entry56941.exe';subprocess.run(['g++','-std=c++20','-O3',str(HERE/'two_entry.cpp'),'-o',str(exe)],check=True)
    p=subprocess.run([str(exe)],input=text,text=True,capture_output=True);need(p.returncode==0,p.stderr)
    hits=[json.loads(s) for s in p.stdout.splitlines()]
    (HERE/'two_entry_bank.json').write_text(json.dumps(dict(status='FINITE_EXACT_BOUNDED_CORES_ALL_TAILS',core_order_bound=10,parent_order_bound=20,cores=cores,hits=hits),indent=2)+'\n')
    print('two-entry hits',hits,flush=True)
if __name__=='__main__':main()
