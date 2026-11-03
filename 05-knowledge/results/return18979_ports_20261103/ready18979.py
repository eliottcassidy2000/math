"""Positive witness search at matching source order; no census claim."""
import sys,subprocess,random,json
from support import HERE,PREV,ROOT,need,save
sys.path.insert(0,str(ROOT/'05-knowledge/results/ready_return_pentagon_20261103'))
import ready_inverse as inv
def main():
    exe=HERE/'ready_filter.exe';subprocess.run(['g++','-std=c++20','-O2',str(PREV/'ready_filter.cpp'),'-o',str(exe)],check=True)
    n=int(sys.argv[1]) if len(sys.argv)>1 else 14
    budget=int(sys.argv[2]) if len(sys.argv)>2 else 1
    bank=inv.two_start(n) if budget==2 else tuple(inv.unrank(n,j) for j in range(inv.class_counts(n)[n]))
    order=list(range(len(bank)));random.Random(18979).shuffle(order)
    print('parent bank',n,'start budget',budget,'classes',len(bank),flush=True);examined=0;hit=None
    for pos in range(0,len(order),256):
        indices=order[pos:pos+256];lines=[' '.join(map(str,(n,budget,j,18979)+bank[j])) for j in indices]
        r=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',capture_output=True,text=True,check=True);examined+=len(indices)
        if r.stdout.strip():
            j,h,q,F=map(int,r.stdout.splitlines()[0].split());T=bank[j]
            need(inv.rooted_data(T,0)==(h,budget,q,F),'independent positive parent check')
            U,root=inv.birth(T,0);z=inv.rooted_data(U,root);need(z[0]==18979 and z[3]-z[0]==budget,'actual ready18979 witness')
            hit={'rows':U,'root':root,'profile':z,'parent_order':n,'budget':budget,'rank':j};break
        if pos%4096==0:print('examined',examined,flush=True)
    name='ready18979.json' if hit else f'search_n{n}_b{budget}.json'
    save(name,{'status':'FINITE_EXACT_POSITIVE_WITNESS' if hit else 'FINITE_EXACT_ABSENT_FROM_DECLARED_BANK','parent_order':n,'budget':budget,'bank_size':len(bank),'examined':examined,'seed':18979,'hit':hit,
      'scope':'classified parent bank at specified order/budget; stop at first positive batch; no all-order absence or native-cost minimality claim'})
    print('PASS ready18979',None if hit is None else hit['profile'],'after',examined,flush=True)
if __name__=='__main__':main()
