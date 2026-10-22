# Parse S50's strong polynomials P_d(z) (strong classes of order j with H = d), odd d < 375, from strong_ears_375.out.
import re, json
src='../h_spectrum_multiplicities_20261021/strong_ears_375.out'
P={}
for line in open(src):
    m=re.match(r'\s+P_(\d+) = (.*)$',line)
    if not m: continue
    d=int(m.group(1)); terms={}
    for c,e in re.findall(r'(\d+)z\^(\d+)',m.group(2)): terms[int(e)]=int(c)
    P[d]=terms
# P_1: K1 only (order 1), plus the transitive? Strong classes with H=1: only K1.  P_3: C3.
json.dump({str(k):v for k,v in P.items()},open('strong_polys_lt375.json','w'))
odd=[d for d in range(1,375,2)]
missing=[d for d in odd if d not in P]
print('parsed',len(P),'polynomials; odd d<375 without a line:',missing[:20])
tot={d:sum(P.get(d,{}).values()) for d in odd}
for d in [3,5,7,9,15,17,21,31,33,63,65,85,127,129,255,257,341]: print(d, P.get(d), 'total', tot.get(d))
