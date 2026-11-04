"""Search an exterior ready mark at the exact (3729,-3181) dyadic target."""
import sys
from support import HERE,PREV,save,need,junction,clone,matching
sys.path.insert(0,str(PREV))
import moment_search as mom

def main():
    hits=[];scanned=[];basehits=[]
    for n in range(3,16):
        for j in range(mom.COUNTS[n]):
            h=mom.moment(n,j,1,0);F=mom.moment(n,j,2,1)
            if h+F-1!=4277:continue
            P=mom.unrank(n,j);T=junction.ear(P,1);r=n;z=junction.profile(T,r)
            need(z[0]==4277 and z[3]-z[0]==1,'literal one-ready base')
            queries=[]
            for q in range(1,n):
                M,C=matching(T,q)
                if len(M)!=1:continue
                U,ids=clone(T,q,2);rr=ids.index((r,0));zz=junction.profile(U,rr)
                queries.append({'clone_vertex':q,'H2':zz[0],'defect2':zz[3]-zz[0]})
                if zz[0]==11735 and zz[3]-zz[0]==1:
                    hits.append({'parent_n':n,'parent_rank':j,'rows':T,'root':r,'clone_vertex':q,'profile1':z,'profile2':zz,'matching':M,'cover':C})
                junction.endpoints.cache_clear()
            basehits.append({'parent_n':n,'rank':j,'rows':T,'profile':z,'matching_one_queries':queries})
            print('4277 parent',n,j,'eligible clones',len(queries),'hits',len(hits),flush=True)
        scanned.append({'parent_n':n,'classes':mom.COUNTS[n]})
    save('exterior_search.json',{'status':'FINITE_EXACT_DECLARED_ONE_READY_BANK','target_formula':'3729*2^m-3181','samples':[4277,11735],
      'parent_banks':scanned,'base_targets':basehits,'uniform_hits':hits,'scope':'all one-start parent classes orders3..15; inspect every exterior matching-one clone except unique in-port; no two-ready or split-incidence absence claim'})
    print('PASS exterior search',len(basehits),'bases',len(hits),'uniform hits',flush=True)
if __name__=='__main__':main()
