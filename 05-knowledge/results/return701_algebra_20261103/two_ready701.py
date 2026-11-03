"""Complete second ready lane at701 under the inherited order reduction."""
from ready701 import inv
from common import save,need

def main():
    hits=[];counts=[]
    for n in range(3,12):
        bank=inv.two_start(n);count=len(bank);counts.append([n,count]);found=0
        for j,T in enumerate(bank):
            h,ss,q,F=inv.rooted_data(T,0)
            if h<=351 and h+F-2==701:
                U,r=inv.birth(T,0);need(inv.rooted_data(U,r)[3]-701==2,'literal2-ready701')
                hits.append({'rows':U,'root':r,'parent':T,'parent_order':n,'rank':j,'parent_profile':[h,ss,q,F]});found+=1
            inv.tables.cache_clear()
        print('two-start parent',n,'classes',count,'ready701',found,flush=True)
    exceptional=[]
    for k in range(1,11):
        exceptional.extend(h for h in (3*2**k-1,6*2**k+1,5*2**k-1,13*2**k+3,14*2**k+3) if h==701)
    need(not exceptional,'five monotone exceptional families omit701')
    save('two_ready701.json',{'status':'FINITE_EXACT_COMPLETE_TWO_READY_INVERSE_AT701','parent_order_bound':11,
       'parent_H_bound':351,'grammar_counts':counts,'hits':hits,'exceptional_hits':exceptional})
if __name__=='__main__':main()
