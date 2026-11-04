"""Test the next two ternary lanes from all five exterior-ready 4277 laws."""
import sys,json,subprocess
from support import HERE,PREV,INCOMING,save,need,junction,clone,matching
sys.path.insert(0,str(PREV))
import moment_search as mom

def main():
    code=(INCOMING/'funnel_bank.cpp').read_text()
    needle='vector<I>targets={18977,18979};';need(code.count(needle)==1,'pinned compiler target site')
    code=code.replace(needle,'vector<I>targets={12833,12835};')
    file=HERE/'next_bank.generated.cpp';file.write_text(code,encoding='utf8',newline='\n')
    exe=HERE/'next_bank.exe';subprocess.run(['g++','-std=c++20','-O3',str(file),'-o',str(exe)],check=True)
    p=subprocess.run([str(exe),'17'],text=True,capture_output=True,check=True)
    bank=[json.loads(line) for line in p.stdout.splitlines()];families=[]
    source=json.loads((HERE/'exterior_search.json').read_text())['base_targets'][1]
    wanted={d:{(3*((q['H2']-4277)//2),3*(8554-q['H2'])+d) for q in source['matching_one_queries']} for d in (2,4)}
    for z in bank:
        if z['kind']!='hit':continue
        n,j,H=z['n'],z['rank'],z['target'];d=H-12831;P=mom.unrank(n,j);T=junction.ear(P,1);r=n
        need(junction.profile(T,r)[0::3]==(H,H+1),'independent literal ready target')
        queries=[]
        for q in range(1,n):
            M,C=matching(T,q)
            if len(M)!=1:continue
            U,ids=clone(T,q,2);rr=ids.index((r,0));h2=junction.H(U);a=(h2-H)//2;b=2*H-h2
            # Only evaluate the costlier root profile for a coefficient hit.
            hit=(a,b) in wanted[d];defect=None
            if hit:
                zz=junction.profile(U,rr);defect=zz[3]-zz[0]
            queries.append({'q':q,'coefficients':[a,b],'coefficient_hit':hit,'defect2_if_hit':defect})
            junction.endpoints.cache_clear()
        families.append({'offset':d,'parent_order':n,'parent_rank':j,'rows':T,'root':r,'H':H,'queries':queries})
        print('next fibre',H,n,j,'clones',len(queries),'hits',sum(q['coefficient_hit'] for q in queries),flush=True)
    save('next_fibres.json',{'status':'FINITE_EXACT_COMPLETE_ONE_READY_FIBRES_USING_CITED_ORDER_BOUND','bank':bank,'targets':families,
        'wanted_coefficients':{d:sorted(v) for d,v in wanted.items()},'scope':'all one-start parents through17; strong parent bound (12835+1)/2=6418<Busch f18=9375; exterior clone queries only, no native source substitution'})
    print('PASS next fibres',len(families),'ready targets',flush=True)
if __name__=='__main__':main()
