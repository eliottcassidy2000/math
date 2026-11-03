"""Recover actual701 ready sources from the inherited complete grammar."""
import sys
from common import HERE,ROOT,save,need
sys.path.insert(0,str(ROOT/'05-knowledge/results/ready_return_pentagon_20261103'))
import ready_inverse as inv

def main():
    hits=[];counts=[]
    for n in range(3,12):
        count=inv.class_counts(n)[n];counts.append([n,count]);found=0
        for j in range(count):
            T=inv.unrank(n,j);h,ss,q,F=inv.rooted_data(T,0)
            if h<=351 and h+F-1==701:
                U,r=inv.birth(T,0);need(inv.rooted_data(U,r)[3]-701==1,'literal1-ready701')
                hits.append({'rows':U,'root':r,'parent':T,'parent_order':n,'rank':j,'parent_profile':[h,ss,q,F]});found+=1
            inv.tables.cache_clear()
        print('parent',n,'classes',count,'ready701',found,flush=True)
    save('ready701.json',{'status':'FINITE_EXACT_COMPLETE_ONE_READY_INVERSE_AT701','parent_order_bound':11,
       'parent_H_bound':351,'grammar_counts':counts,'hits':hits})
if __name__=='__main__':main()
