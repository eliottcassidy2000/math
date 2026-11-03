"""Find literal ready targets; bounded search, not an all-size absence claim."""
from support import HERE,ROOT,need,save
import sys
sys.path.insert(0,str(ROOT/'05-knowledge/results/ready_return_pentagon_20261103'))
import ready_inverse as inv

def main():
    hits=[];tested=[]
    for n in range(3,11):
        for budget,bank in ((1,[inv.unrank(n,j) for j in range(inv.class_counts(n)[n])]),(2,inv.two_start(n))):
            count=0
            for j,T in enumerate(bank):
                h,ss,q,F=inv.rooted_data(T,0);count+=1
                if h+F-budget==2107:
                    U,r=inv.birth(T,0);profile=inv.rooted_data(U,r)
                    need(profile[0]==2107 and profile[3]-2107==budget,'literal ready2107')
                    hits.append({'rows':U,'root':r,'parent_order':n,'budget':budget,'rank':j,'profile':profile})
                inv.tables.cache_clear()
            tested.append([n,budget,count]);print('parent',n,'budget',budget,'count',count,'hits',len(hits),flush=True)
        if hits:break
    save('ready2107.json',{'status':'FINITE_EXACT_SHORTEST_WITHIN_SCANNED_GRAMMARS','tested':tested,'hits':hits,
        'scope':'all smaller orders in one/two-start grammars scanned; exceptional families and larger parent orders not needed for witness'})
if __name__=='__main__':main()
