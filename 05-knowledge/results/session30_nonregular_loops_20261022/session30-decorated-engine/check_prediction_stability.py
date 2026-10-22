import json,numpy as np
from predict import predict
out=[]
for q in [7,11]:
 r={i*i%q for i in range(1,q)};n=q-1
 S=np.array([[0 if i==j else 1 if (j-i)%q in r else -1 for j in range(n)] for i in range(n)])
 a=predict(S);b=predict(S,derivative_radius=.05,derivative_points=32)
 dif={k:abs(a[k]-b[k]) for k in ['D_A1','D2_A1','y2','kappa2']}
 if max(dif.values())>2e-9:raise RuntimeError(dif)
 out.append({'q':q,'differences':dif})
print(json.dumps({'status':'PASS','checks':out},indent=2))
