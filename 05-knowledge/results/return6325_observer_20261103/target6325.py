"""Exact arc-distance alignment in the single bounded ready-target bank."""
import json,subprocess
from support import HERE,source,need,save,junction
def main():
    A,r,b,names=source();bank=json.loads((HERE/'ready6325.json').read_text())['hits'];need(len(bank)==1,'declared target bank')
    d=bank[0];T=tuple(d['rows']);tr=d['root'];n=len(A);need(len(T)==n,'same source order')
    si=[v for v in range(n) if A[v]>>r&1];ti=[v for v in range(n) if T[v]>>tr&1]
    need(len(si)==len(ti)==1,'fixed root and unique in-port')
    exe=HERE/'match6325.exe';subprocess.run(['g++','-std=c++20','-O2',str(HERE/'match6325.cpp'),'-o',str(exe)],check=True)
    args=(n,r,si[0],tr,ti[0])+A+T
    run=subprocess.run([str(exe)],input=' '.join(map(str,args)),capture_output=True,text=True,check=True)
    lines=run.stdout.splitlines();cost,nodes=map(int,lines[0].split());perm=tuple(map(int,lines[1].split()))
    check=subprocess.run([str(exe),'--without-degree-bound'],input=' '.join(map(str,args)),capture_output=True,text=True,check=True)
    other,other_nodes=map(int,check.stdout.splitlines()[0].split());need(other==cost,'independent pruning control')
    need(sorted(perm)==list(range(n)) and perm[r]==tr and perm[si[0]]==ti[0],'named alignment')
    U=junction.sub(T,perm);edits=[(i,j) for i in range(n) for j in range(i+1,n) if ((A[i]>>j)^(U[i]>>j))&1]
    need(len(edits)==cost and junction.profile(U,r)==tuple(d['profile']),'exact aligned ready source')
    save('target6325.json',{'status':'FINITE_EXACT_MINIMUM_IN_ONE_SPECIFIED_ROOTED_TARGET_CLASS','nodes':nodes,'unbounded_degree_control_nodes':other_nodes,'cost':cost,
      'target_to_source_order':perm,'final_rows':U,'flips':edits,'source_names':names,'root':r,'spare':b,'fixed_in_port':si[0],
      'scope':'all root/in-port-fixed identifications of this bounded-search target; not all ready targets or all native programs'})
    print('PASS target6325 minimum arcs',cost,'search nodes',nodes,flush=True)
if __name__=='__main__':main()
