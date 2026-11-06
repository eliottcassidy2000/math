"""Check every recipe by exact transfer, and emit literal rows on request."""
import json,sys
missing={7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647}
seen=set(); largest_order=0
for line in open('WITNESSES.tsv'):
 target,word=line.split();target=int(target)
 if target in seen:raise ValueError('Duplicate target')
 seen.add(target);h,r,e,d=1,2,0,1;out=[2,0];root=0
 for c in word.replace('-',''):
  n=len(out)
  if c=='O':
   for i in range(n):out[i]|=1<<n
   out.append(0);h,r,e,d=h,2*r,0,d
  elif c in 'ZJ':
   for i in range(n):
    if i!=root:out[i]|=1<<n
   out.append(1<<root)
   if c=='Z':h,r,e,d=h+r-e,2*r-d,d,d
   else:h,r,e,d=h+r-e,h+r,h-e,h;root=n
  else:raise ValueError('Invalid letter')
 for i in range(len(out)):
  if out[i]>>i&1:raise ValueError('Loop')
  for j in range(i):
   if (out[i]>>j&1)+(out[j]>>i&1)!=1:raise ValueError('Invalid orientation')
 if h!=target:raise ValueError('Count mismatch')
 largest_order=max(largest_order,len(out))
 if len(sys.argv)>1 and target==int(sys.argv[1]):print(json.dumps({'target':target,'root':root,'rows':out,'word':word}))
expected=set(range(1,50000,2))-missing
if seen!=expected:raise ValueError('Incorrect scalar coverage')
print(json.dumps({'status':'PASS','recipes':len(seen),'max_order':largest_order,'literal_rows_validated':len(seen),'counts_verified_by':'inherited exact physical transfer; not literal subset DP for every large witness'}))
