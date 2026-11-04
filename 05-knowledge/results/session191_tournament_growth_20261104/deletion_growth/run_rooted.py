import itertools,json,pathlib,subprocess,collections,sys,time
P=pathlib.Path(__file__).resolve().parent
cap=int(sys.argv[1])if len(sys.argv)>1 else 44
end=int(sys.argv[2])if len(sys.argv)>2 else 10
D=P/('rooted_cap'+str(cap));D.mkdir(exist_ok=True)
def canonical(r):
 g=collections.defaultdict(list)
 for i in range(1,len(r)):g[r[i].bit_count()].append(i)
 blocks=[list(itertools.permutations(g[k]))for k in sorted(g)]
 return min(tuple(sum(((r[x]>>y)&1)<<j for j,y in enumerate(p))for x in p)for pp in itertools.product(*blocks)for p in[(0,)+sum(pp,())])
def strong(r,S):
 n=len(r)
 for v in range(n):
  if not(S>>v&1):continue
  reach=1<<v
  while True:
   nxt=reach
   for i in range(n):
    if reach>>i&1:nxt|=r[i]&S
   if nxt==reach:break
   reach=nxt
  if reach!=S:return False
 return True
# Exhaust the 64 labelled tournaments of order four, fixing the marked root at 0.
base=set()
pairs=list(itertools.combinations(range(4),2))
for bits in range(1<<len(pairs)):
 r=[0]*4
 for k,(a,b) in enumerate(pairs):
  if bits>>k&1:r[a]|=1<<b
  else:r[b]|=1<<a
 if strong(r,15) and strong(r,14):base.add(canonical(r))
json.dump(sorted(base),open(D/'base4.json','w'),indent=2)
current={r:(5,3)for r in base}
summary={'cap':cap,'basis_order4_count':len(base),'levels':[]}
for n in range(4,end):
 inp=D/('bank'+str(n)+'.txt');out=D/('raw'+str(n+1)+'.txt');log=D/('step'+str(n+1)+'.log')
 with open(inp,'w')as f:
  f.write(str(len(current))+'\n')
  for r in sorted(current):f.write(' '.join(map(str,r))+'\n')
 with open(inp)as fi,open(out,'w')as fo,open(log,'w')as fe:subprocess.run([str(P/'rooted_ears'),str(n),str(cap)],stdin=fi,stdout=fo,stderr=fe,check=True)
 nxt={}
 for line in open(out):
  h,h0,*r=map(int,line.split());c=canonical(r)
  if c in nxt and nxt[c]!=(h,h0):raise RuntimeError('canonical count mismatch')
  nxt[c]=(h,h0)
 current=nxt
 row={'order':n+1,'classes':len(nxt),'minimum':min((h-h0 for h,h0 in nxt.values()),default=None),'maxH':max((h for h,h0 in nxt.values()),default=None),'counts':dict(sorted(collections.Counter(h-h0 for h,h0 in nxt.values()).items())),'extension_log':log.read_text()}
 summary['levels'].append(row);print(row,flush=True)
 json.dump([{'rows':r,'H':hh[0],'Hminusroot':hh[1]}for r,hh in sorted(nxt.items())],open(D/('bank'+str(n+1)+'.json'),'w'),indent=2)
 json.dump(summary,open(D/'summary.json','w'),indent=2)
 if not nxt:break
