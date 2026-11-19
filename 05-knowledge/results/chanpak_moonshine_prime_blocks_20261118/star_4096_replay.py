try:
    from star_factor_replay import cached_profile, bank, need, H, star
except ImportError:
    from _replay_base import cached_profile, bank, need, H, star
from itertools import product
from pathlib import Path
bank=dict(bank)
for row in Path(__file__).with_name('star_seam_bank_4096.tsv').read_text().splitlines():
    if row.startswith('#') or not row.strip():continue
    q,c,a,A,b,B=row.split()
    q,c,a,b=map(int,(q,c,a,b))
    A=tuple(map(int,A.split(','))) if A!='-' else ()
    B=tuple(map(int,B.split(','))) if B!='-' else ()
    ha,ka,ea,oa=cached_profile(A,a)
    hb,kb,eb,ob=cached_profile(B,b)
    need(ha*(kb-2*eb)+hb*(ka-ha)==1,(q,A,B))
    need(c==(ha+ea+oa)*(eb+ob) and c//(c&-c)==q,(q,c))
    need(q not in bank,q)
    bank[q]=c
seen=set()
for n in range(1,(2*4096+1).bit_length()):
    for p in product('01',repeat=n-1):
        h,e,t=1,0,2
        for b in (*p,'0','1'):
            h,e,t=(h+t-e,1,2*t-1) if b=='1' else (h,0,2*t)
        if t-h<=4096:
            c=t-1
            q=c//(c&-c)
            need(q in bank,(t-h,p,q))
            need(bank[q]+1<2**32,(q,bank[q]))
            seen.add(q)
need(seen==set(bank), (len(seen),len(bank)))
print('PASS: all',len(seen),'slope classes through defect 4096; largest inherited early term',max(bank.values())+1)

# Recompute the one-arm endpoint equation without using apex counts
# in the weighted even, double-tagged path enumeration.
B=(24,1,3,6,14)
from itertools import permutations
for tags,p in [(6,29),(12,37),(24,47)]:
    h,k,e,o=cached_profile(B,tags)
    J=0
    for r in (2,4):
        for path in permutations(range(5),r):
            if not(tags>>path[0]&1 and tags>>path[-1]&1):continue
            if not all(B[u]>>v&1 for u,v in zip(path,path[1:])):continue
            outside=[v for v in range(5) if v not in path]
            J+=H([sum(1<<i for i,v in enumerate(outside) if B[u]>>v&1)
                  for u in outside])
    need((h,e+o,J)==(11,p,5) and k==h+2*(e-J),(tags,h,k,e,o,J))
print('PASS: double-tagged path identity on slopes 29, 37, 47')
