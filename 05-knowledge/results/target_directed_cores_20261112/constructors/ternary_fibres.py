"""Focused same-base source-change gate; every reported universe is bounded."""
from pathlib import Path
from itertools import product
import json
from paired_scale_probe import evaluate, inflate
from selective_split import compiler
HERE=Path(__file__).resolve().parent

def main():
 ns=compiler();bank={1:[(ns['SINGLE'],())]};wanted={431,551,893};fibres=[];cases=[];meta=[]
 for n in range(3,13):
  bank[n]=[]
  for m in (1,*range(3,n)):
   t=n-m-1;words=[''] if t==0 else [''.join(w)+'1' for w in product('01',repeat=t-1)]
   for S,address in bank[m]:
    for word in words:
     N=ns['advance'](S,word);a=address+(word,);bank[n].append((N,a));h=N[0]+N[4]-1
     if h not in wanted:continue
     T=ns['make'](a);fibres.append({'h':h,'address':a,'rows':T})
     for v in range(n-1):
      sizes=[1]*len(T);sizes[v]=2;R,_=inflate(T,sizes);cases.append((R,-1));meta.append((h,a,T,v))
 hits=[]
 for (h,a,T,v),(h2,_,_) in zip(meta,evaluate(cases)):
  if h2!=3*h+2:continue
  sizes=[1]*len(T);sizes[v]=3;R,g=inflate(T,sizes)
  h3,q,s=evaluate([(R,g[-1][0])])[0]
  hits.append({'h':h,'address':a,'rows':T,'module':v,'H2':h2,'profile3':[h3,q,s],'pure3_at3':h3==9*(h+1)-1})
 out={'status':'FINITE-EXACT','scope':'Strong unique-start parents through12, target base431/551/893, every module avoiding readyroot/predecessor','fibres':fibres,'tests':len(cases),'hits':hits}
 (HERE/'pure_ternary_fibres_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'fibres':len(fibres),'tests':len(cases),'first_step_hits':len(hits),'third_step_hits':sum(x['pure3_at3'] for x in hits)}))

if __name__=='__main__':main()
