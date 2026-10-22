import numpy as np,json
from decorated import evaluate
out=[]
for name,n in [('transitive4',4),('transitive6',6),('regular7',7),('deleted_regular7',6)]:
 if 'transitive' in name:S=np.array([[0 if i==j else 1 if j>i else -1 for j in range(n)] for i in range(n)])
 else:
  S=np.array([[0 if i==j else 1 if 0<(j-i)%7<=3 else -1 for j in range(7)] for i in range(7)])[:n,:n]
 out.append({'name':name,'S':S.tolist(),'terms':{str(w):[dict(d,value=float(v)) for d,v in evaluate(S/n,w,detail=True)] for w in [1,2]}})
print(json.dumps(out,indent=2))
