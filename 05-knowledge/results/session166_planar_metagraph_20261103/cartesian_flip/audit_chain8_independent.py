"""Independent bit-action, endpoint-DP, and separator audit of the G5 witness."""
import json
from itertools import combinations,permutations
from pathlib import Path
P=Path(__file__).resolve().parent
CERT=json.loads((P/'chain8_certificate.json').read_text())
MASKS=[17928585,9314500,10846530,9320707,1712324,3244354,1718531,54952202,25475209,3342737,54542928,67668801,83954712,54804816,67684801,84478232]
def req(x,msg):
    if not x:raise RuntimeError(msg)
req(CERT['masks']==MASKS,'mask list')
pairs=list(combinations(range(8),2));positions={p:i for i,p in enumerate(pairs)}
actions=[]
for p in permutations(range(8)):
    toggle=0;images=[]
    for a,b in pairs:
        x,y=p[a],p[b];bit=1<<positions[tuple(sorted((x,y)))];images.append(bit)
        if x>y:toggle|=bit
    actions.append((toggle,images))
def image(t,act):
    u,imgs=act
    while t:
        z=t&-t;u^=imgs[z.bit_length()-1];t-=z
    return u
expected=set()
for i in range(5):
    t,a,b=1+3*i,2+3*i,3+3*i
    expected|={tuple(sorted(e))for e in [(0,a),(0,b),(a,b),(a,t),(b,t)]}
    if i:expected.add((t-3,t))
req(len(expected)==29,'edge count')
dist=[[0]*16 for _ in MASKS]
for i,mask in enumerate(MASKS):
    orbit=[image(mask,act)for act in actions]
    for j in range(i+1,16):
        d=min((x^MASKS[j]).bit_count()for x in orbit);dist[i][j]=dist[j][i]=d
        req(d>0,'nonisomorphism')
        req((d==1)==((i,j)in expected),'pair adjacency')
req(dist==CERT['minimum_quotient_hamming_distances'],'certificate distance matrix')
H=[]
for ti,mask in enumerate(MASKS):
    out=[0]*8
    for k,(a,b)in enumerate(pairs):
        if mask>>k&1:out[a]|=1<<b
        else:out[b]|=1<<a
    reach=[out[i]|1<<i for i in range(8)]
    for k in range(8):
        for i in range(8):
            if reach[i]>>k&1:reach[i]|=reach[k]
    req(all(x==255 for x in reach),'strong connectivity')
    dp=[[0]*8 for _ in range(256)]
    for i in range(8):dp[1<<i][i]=1
    for subset in range(1,256):
        for last in range(8):
            available=out[last]&~subset
            while available:
                bit=available&-available;available-=bit
                dp[subset|bit][bit.bit_length()-1]+=dp[subset][last]
    hs=[sum(row)for row in dp];H.append(hs[-1])
    covers=[[0]*9 for _ in range(256)];covers[0][0]=1
    for subset in range(1,256):
        first=subset&-subset;part=subset
        while part:
            if part&first:
                for k in range(1,9):covers[subset][k]+=hs[part]*covers[subset^part][k-1]
            part=(part-1)&subset
    req(covers[-1]==CERT['profiles'][ti],'independent full profile')
    req(hs[-1]>0 and hs[-1]%2==1,'tournament Hamilton path count')
# The undirected sixteen-class graph has no Hamilton path: a three-vertex separator leaves five components.
adj=[set()for _ in range(16)]
for a,b in expected:adj[a].add(b);adj[b].add(a)
removed={0,4,10};remaining=set(range(16))-removed;components=[]
while remaining:
    todo=[min(remaining)];comp=set(todo);remaining-=comp
    while todo:
        v=todo.pop();new=adj[v]&remaining;remaining-=new;comp|=new;todo+=list(new)
    components.append(sorted(comp))
req(len(components)==5,'nontraceability separator')
coloring=[0]*16
for i in range(5):coloring[1+3*i]=3*(i%2);coloring[2+3*i]=1;coloring[3+3*i]=2
req(all(coloring[a]!=coloring[b]for a,b in expected),'four coloring')
# Any 3-coloring of a diamond makes its opposite tips equal. Adjacent chain tips contradict this.
report={'status':'PASS','pairs':120,'permutations_per_source':40320,'strong_tournaments':16,'induced_edges':29,'tournament_H':H,'separator':[0,4,10],'components_after_deletion':components,'metagraph_has_Hamilton_path':False,'proper_four_coloring':coloring,'independent_profile_algorithms':'endpoint subset DP plus first-block cover recurrence'}
(P/'audit_chain8_independent.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
