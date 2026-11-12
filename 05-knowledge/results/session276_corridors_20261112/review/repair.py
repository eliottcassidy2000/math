import contextlib,io,json
with contextlib.redirect_stdout(io.StringIO()):
 from independent import subset_count,endpoint,cover
 from sandwich import one,fuse

def mat(rows):return [[bool(r>>j&1) for j in range(len(rows))] for r in rows]
def flip(G,edges):
 G=[r[:] for r in G]
 for a,b in edges:G[a][b],G[b][a]=G[b][a],G[a][b]
 return G
T43=mat((52,13,24,49,34,6));T45=mat((52,13,24,33,42,6))
P=fuse(one(2),T43);T=fuse(one(2),T45);R=flip(T,[(2,3),(3,7)])
rows=[sum(1<<j for j,b in enumerate(r) if b) for r in R]
assert rows==[210,252,241,181,96,129,168,16]
assert [subset_count(G) for G in [P,T,R]]==[157,159,159]
for X,Y,p in [(T,R,(1,2,5,7,4,6,3,0)),(P,R,(0,1,2,4,6,3,5,7))]:
 assert all(G[a][b] for G in [X,Y] for a,b in zip(p,p[1:]))
U=endpoint(R);V=endpoint(R,True);core=[r[1:] for r in R[1:]];h=subset_count(core);cs=cover(core)
profile=dict(H=subset_count(R),h=h,Q=2*h+2*cs[2],U=sum(U),V=sum(V),d=sum((-1)**(j+1)*x for j,x in enumerate(U)),e=sum((-1)**(j+1)*x for j,x in enumerate(V)))
assert profile==dict(H=159,h=17,Q=448,U=312,V=213,d=8,e=17)
changed=[(i,j) for i in range(8) for j in range(i+1,8) if P[i][j]!=R[i][j]]
assert changed==[(2,3),(3,7),(5,6)]
print(json.dumps({'status':'PASS','repaired_rows':rows,'profile':profile,'same_V_direct_edge_changes':changed,'scope':'finite repair; three-free-edge enlarged move only for157toR159'},indent=2))
