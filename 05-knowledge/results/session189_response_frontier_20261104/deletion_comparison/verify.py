import collections,hashlib,json,math,pathlib
P=pathlib.Path(__file__).resolve().parent
def need(b,msg):
    if not b:raise RuntimeError(msg)
bank=json.loads((P/'bank.json').read_text());cards=json.loads((P/'literal_cards.json').read_text());reps=json.loads((P/'representatives.json').read_text())
need(len(bank)==len(cards)==len(reps['response_cases'])==12,'case count')
raw=list(map(int,(P/'independent_values.txt').read_text().split()))
need(len(raw)==25392,'independent query count')
source=(P/'source.out').read_text().split('\nreal\t')[0].strip()
need((P/'reproduced.out').read_text().strip()==source,'exact census output mismatch')
meta=json.loads((P/'source_metadata.json').read_text())
need(meta['hypothesis']['sha']==meta['hypothesis']['latest_sha'],'hypothesis changed at source check')
for name,key in [('source.c','code'),('source.out','output'),('HYP-9455.md','hypothesis')]:
    data=(P/name).read_bytes();blob=b'blob '+str(len(data)).encode()+b'\0'+data
    need(hashlib.sha1(blob).hexdigest()==meta[key]['sha'],'pinned Git blob mismatch '+name)
hist=[];degree=collections.Counter();pos=0
for g,lit,rep in zip(bank,cards,reps['response_cases']):
    need(g['rows']==rep['rows'],'representative rows')
    n=len(g['rows']);need(n==9,'old order')
    for v,row in enumerate(g['rows']):
        need(not(row>>v&1),'loop')
        for w in range(v):need(((row>>w)&1)+((g['rows'][w]>>v)&1)==1,'not tournament')
    need(g['coefficients'][0]==lit['H']==int(rep['source_line'].split()[3]),'independent Hamilton count')
    need(lit['aut']==int(rep['source_line'].split()[2]),'independent automorphism count')
    need(lit['self_converse']==int(rep['source_line'].split()[1]),'independent self-converse status')
    masks=g['masks'];L=len(masks);need(L==46,'complete support basis')
    vals=raw[pos:pos+L*L];pos+=L*L;need(vals==g['values'],'independent predecessor counts')
    index={m:i for i,m in enumerate(masks)};coefs=[];h=collections.defaultdict(collections.Counter)
    for i,a in enumerate(masks):
        for j,c in enumerate(masks):
            z=sum((-1)**((a^s).bit_count()+(c^t).bit_count())*vals[index[s]*L+index[t]]for s in masks if s&a==s for t in masks if t&c==t)
            coefs.append(z);degree[a.bit_count()+c.bit_count()]+=1
            h[a.bit_count(),c.bit_count(),(a&c).bit_count()][z]+=1
    need(coefs==g['coefficients'],'full coefficient inversion')
    hist.append(h)
def deck(g):return collections.Counter(x['merged']for x in g['cards'])
def weighted(g):
    z=collections.Counter()
    for x in g['cards']:z[x['merged']]+=x['starts']+x['ends']
    return z
def primitive(z):
    d=math.gcd(*z.values());return {k:v//d for k,v in z.items()}
comparisons=[];invariants=[]
for j,pair in enumerate(reps['pairs']):
    a,b=cards[4*j],cards[4*j+2]
    need(deck(a)==deck(b),'independent merged deck collision')
    need(a['merged']!=b['merged'],'not distinct merged nodes')
    need(primitive(weighted(a))!=primitive(weighted(b)),'primitive weighted decks collide')
    for u in range(2):
        for v in range(2):
            i,k=4*j+u,4*j+2+v
            need(hist[i][0,1,0]!=hist[k][0,1,0],'linear c histogram must separate')
            dif=[list(key)for key in sorted(hist[i])if hist[i][key]!=hist[k][key]]
            comparisons.append(dict(pair=pair,converse_choices=[u,v],different_histogram_types=dif,all_old_name_bijections_excluded=True))
    for k in range(4*j,4*j+4):
        g=bank[k];L=len(g['masks']);invariants.append(dict(node=reps['response_cases'][k]['node'],converse=reps['response_cases'][k]['converse'],H=g['coefficients'][0],linear_c_coefficients=sorted(g['coefficients'][t]for t,m in enumerate(g['masks'])if m.bit_count()==1)))
    for i,k in [(4*j,4*j+1),(4*j+2,4*j+3)]:
        need(cards[i]['merged']==cards[k]['merged']and deck(cards[i])==deck(cards[k])and weighted(cards[i])==weighted(cards[k]),'converse invariance')
result=dict(status='PASS',census_reproduction='All non-timing lines through order 9 match exactly; this reuses incoming census code.',graphs=12,independent_old_card_canonicalizations=108,full_named_coefficients=len(raw),degree_counts=dict(sorted(degree.items())),orientation_comparisons=12,comparisons=comparisons,invariants=invariants,scope='The three reported order-9 merged unweighted-deck collision pairs are excluded as full named two-port collisions under all relabelings and both converse choices. No global reconstruction theorem is asserted.')
print(json.dumps(result,indent=2,sort_keys=True))
