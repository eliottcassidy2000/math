"""Complete fixed-address ternary response probe through core order seven."""
from search_central import *
def main():
 old=[(0,)];bank={};records=[];class_counts={}
 for n in range(1,8):
  if n>1:old=expand(old)
  class_counts[str(n)]=len(old)
  for rows in old:
   for O,s in cuts(rows):
    key=s[:7];q=s[7]+s[8];w={'rows':rows,'O':O,'signature':s}
    records.append((key,q,w));bank.setdefault(key,{}).setdefault(q,w)
 found=[];scaled=[]
 for key,q,w in records:
  if all(v%3==0 for v in key):
   small=tuple(v//3 for v in key)
   for qs,ws in bank.get(small,{}).items():
    delta=q-3*qs
    scaled.append(delta)
    if delta in (2,4,-2,-4):found.append({'delta':delta,'source':ws,'target':w})
 result={'scope':'All labelled cuts of one representative of each tournament isomorphism class through order7, both source and target; no source-order restriction. F_target(A,B)=3F_source(A,B)+delta for the same binary addresses. Exact seven-coefficient comparison.', 'class_counts':class_counts,'core_cut_records':len(records),'scaled_response_matches':len(scaled),'delta_histogram':dict(sorted(collections.Counter(scaled).items())),'nonzero_unit_ternary_lifts':found}
 (HERE/'scaled_search.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
