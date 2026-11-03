"""Bounded ready-target search using exact C++ endpoint counts and Python controls."""
import sys,subprocess,json
from support import HERE,ROOT,need,save
sys.path.insert(0,str(ROOT/'05-knowledge/results/ready_return_pentagon_20261103'))
import ready_inverse as inv

def main():
    exe=HERE/'ready_filter.exe'
    subprocess.run(['g++','-std=c++20','-O2',str(HERE/'ready_filter.cpp'),'-o',str(exe)],check=True)
    hits=[];tested=[]
    for n in range(3,13):
        for b,bank in ((1,[inv.unrank(n,j) for j in range(inv.class_counts(n)[n])]),(2,inv.two_start(n))):
            lines=[' '.join(map(str,(n,b,j,6325)+tuple(T))) for j,T in enumerate(bank)]
            r=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',capture_output=True,text=True,check=True)
            for line in r.stdout.splitlines():
                j,h,q,F=map(int,line.split());T=bank[j]
                need(inv.rooted_data(T,0)==(h,b,q,F),'independent Python profile')
                U,root=inv.birth(T,0);z=inv.rooted_data(U,root)
                need(z[0]==6325 and z[3]-z[0]==b,'literal ready target')
                hits.append({'rows':U,'root':root,'profile':z,'parent_order':n,'budget':b,'rank':j})
                inv.tables.cache_clear()
            tested.append([n,b,len(bank)]);print('parent',n,'budget',b,'count',len(bank),'hits',len(hits),flush=True)
        if hits:break
    save('ready6325.json',{'status':'FINITE_EXACT_BOUNDED_GRAMMAR_SEARCH','tested':tested,'hits':hits,
      'scope':'one/two-start parent grammars only; no all-size absence or minimum native cost claim'})
if __name__=='__main__':main()
