import itertools,json
out={}
for p in (19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101):
 sq={r*r%p for r in range(1,p)}
 found=None
 for C,x,y in itertools.product(range(p),range(1,p),range(1,p)):
  vals={(C+i*x+j*y)%p for i,j in itertools.product((-1,0,1),repeat=2)}
  if len(vals)==9 and vals<=sq:
   A=[(-x)%p,0,x];B=[(C-y)%p,C,(C+y)%p]
   bad={((b-a)*pow(2,-1,p))%p for a,b in itertools.product(A,B)}
   good=set(range(p))-bad
   assert len(bad)==9
   t=min(good);Ag=[(a+t)%p for a in A];Bg=[(b-t)%p for b in B]
   assert len(set(Ag+Bg))==6
   found={'C':C,'x':x,'y':y,'nine_square_values':sorted(vals),'roots_by_value':{v:min(r for r in range(1,p) if r*r%p==v) for v in sorted(vals)},'bad_gauges':sorted(bad),'valid_gauge_count':len(good),'chosen_gauge':t,'A':Ag,'B':Bg}
   break
 out[p]=found
 if found:break
print(json.dumps(out,indent=2));open(__import__('pathlib').Path(__file__).with_name('finite-field-results.json'),'w').write(json.dumps(out,indent=2))
