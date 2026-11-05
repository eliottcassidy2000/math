import json, collections, sys
cap=int(sys.argv[1]) if len(sys.argv)>1 else 2500
t=(1,0,2,1)
def children(v):
 D,E,U,d=v
 return ((D,0,2*U,d),(D-E+U,d,2*U-d,d),(D-E+U,D-E,D+U,D))
def F(v):D,E,U,d=v;return D-E+U
bank={t}; todo=[t]
for v in todo:
 for w in children(v):
  if w[2]<=cap and w not in bank:bank.add(w);todo.append(w)
fibers=collections.defaultdict(list)
for v in bank:
 n=F(v)
 if n<=cap+1:fibers[n].append(v)
# F=U+(D-E)>=U+1, so this threshold is a complete fiber census.
for v in bank:
 o,z,j=children(v)
 if F(o)!=F(j) or o==j:raise RuntimeError(('two-profile failure',v))
result={'U_cap':cap,'profiles':len(bank),'complete_F_target_cap':cap+1,'distinct_F_values':len(fibers),'singleton_F_fibers':{str(k):v for k,v in sorted(fibers.items()) if len(v)==1},'two_profile_identity_checked':len(bank)}
print(json.dumps(result,sort_keys=True,indent=2))
