# Residual analysis of P_d(1) = total strong classes with H = d (all orders), odd d < 375, vs a local trend.
import json, math
tab={int(k):{int(a):b for a,b in v.items()} for k,v in json.load(open('strong_table_lt375.json')).items()}
tot={d:sum(v.values()) for d,v in tab.items()}
ds=sorted(tot)
def trend(d,frac=0.15):
    w=[tot[e] for e in ds if abs(e-d)<=frac*d and e!=d]
    return sum(w)/len(w) if w else None
rows=[]
for d in ds:
    if d<25: continue
    t=trend(d)
    r=math.log((tot[d]+0.5)/(t+0.5))
    rows.append((r,d,tot[d],round(t,1)))
rows.sort()
print('most negative log-residuals (d, P_d(1), local mean):')
for r,d,v,t in rows[:30]: print(f'  d={d:4d} total={v:4d} mean={t:6.1f} logres={r:+.2f}  Q={(d-1)//2:b}  d={d:b}')
print('ranks of special values (rank 1 = rarest relative to trend) of', len(rows))
rank={d:i+1 for i,(r,d,v,t) in enumerate(rows)}
for d in [31,63,127,255,85,341,33,65,129,257,43,171,21*4+1]:
    if d in rank: print(f'  d={d} rank {rank[d]}  total {tot[d]}')
