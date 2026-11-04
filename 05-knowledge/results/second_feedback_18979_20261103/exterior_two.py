"""Actual two-ready targets and matching-one exterior clone checks."""
import json,sys,subprocess
from support import HERE,PREV,INCOMING,save,need,junction,clone,matching
sys.path.insert(0,str(HERE.parent/'ready_return_pentagon_20261103'))
import ready_inverse as inv

def main():
    cores=json.loads((INCOMING/'two_entry_bank.json').read_text())['cores']
    exe=HERE/'two_entry4277.exe';subprocess.run(['g++','-std=c++20','-O3',str(HERE/'two_entry4277.cpp'),'-o',str(exe)],check=True)
    inp=''.join(' '.join(map(str,(len(c['rows']),sum(1<<j for j in c['ports']),i,*c['rows'])))+'\n' for i,c in enumerate(cores))
    p=subprocess.run([str(exe)],input=inp,capture_output=True,text=True,check=True)
    bases=[];hits=[]
    for line in p.stdout.splitlines():
        i,t,bits,h,F=map(int,line.split());c=cores[i];word=tuple((bits>>j)&1 for j in range(t))
        P,w=inv.funnel(tuple(c['rows']),set(c['ports']),word);T,r=inv.birth(P,w);z=junction.profile(T,r)
        need(z[0]==4277 and z[3]-z[0]==2,'literal two-ready target4277')
        query=[]
        for q in range(len(T)):
            if q in (r,w):continue
            M,C=matching(T,q)
            if len(M)!=1:continue
            U,ids=clone(T,q,2);rr=ids.index((r,0));zz=junction.profile(U,rr)
            query.append({'q':q,'H2':zz[0],'defect2':zz[3]-zz[0]})
            if zz[0]==11735 and zz[3]-zz[0] in (1,2):hits.append({'rows':T,'root':r,'clone_vertex':q,'profile1':z,'profile2':zz,'matching':M,'cover':C})
            junction.endpoints.cache_clear()
        bases.append({'core_id':i,'tail':word,'rows':T,'root':r,'profile':z,'queries':query})
    save('exterior_two.json',{'status':'FINITE_EXACT_BOUNDED_TWO_ENTRY_CORE_BANK','core_count':len(cores),'core_bound':9,'parent_bound':15,'targets':bases,'uniform_hits':hits,
      'scope':'all legal binary tails on saved two-entry core presentations through order9; no complete all-order two-ready exclusion'})
    print('PASS exterior two-ready search',len(bases),'bases',len(hits),'uniform hits',flush=True)
if __name__=='__main__':main()
