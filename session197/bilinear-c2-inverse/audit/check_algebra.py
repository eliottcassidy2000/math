import json

def need(b,s):
 if not b:raise RuntimeError(s)
def raw(k,v,w):
 D,E,U,d=v;g,s,L,e=w
 if k=='N':h,A,R,K=D-E+U,D-E+d,D+3*U-d,U+E
 else:h,A,R,K=D+2*U,U+d,2*D-E+6*U-d,U+E
 return h*g-A*s+R*L-K*e
basis=((1,0,1,0),(0,0,1,0),(1,1,1,0),(1,0,1,1))
matrices={k:[[raw(k,v,w) for w in basis] for v in basis] for k in ('N','M')}
need(matrices['N']==[[6,4,5,5],[4,3,4,3],[5,4,5,3],[5,3,3,4]],'exact N coefficients')
need(matrices['M']==[[11,8,10,10],[8,6,7,7],[10,7,9,8],[10,7,8,9]],'exact M coefficients')
for i in range(4):
 for j in range(4):
  need(matrices['N'][i][j]-3-(i==j==3)>=0,'N positive remainder')
  need(matrices['M'][i][j]-6-(i==3)-(j==3)-(i==j==3)>=0,'M positive remainder')
K=120;b=[0]*(K+1);b[1]=1
for n in range(3,K+1):b[n]=b[n-1]+sum((1<<(n-m-2))*b[m] for m in range(1,n-1))
W=[sum(b[:n+1]) for n in range(K+1)];a=[0]*(K+1);c=[0]*(K+1)
for n in range(K+1):
 a[n]=sum(W[i]*W[j] for add in (1,2) for i in range(2,n) for j in [n-add-i] if j>=2)
 c[n]=(6*c[n-1] if n>=1 else 0)-(11*c[n-2] if n>=2 else 0)+(6*c[n-3] if n>=3 else 0)-(c[n-4] if n>=4 else 0)+({5:1,6:-1,7:-1,8:1}.get(n,0))
 need(a[n]==c[n],'ordered presentations match independent nonregular C2 GF')
print(json.dumps({'status':'PASS','matrix_entries':32,'positive_remainder_entries':32,'graph_class_coefficients':K+1,'nonregular_C2_orders_5_through_11':a[5:12],'weak_rooted_source_counts_1_through_8':W[1:9]},sort_keys=True,indent=2))
