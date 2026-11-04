#!/usr/bin/env python3
"""Literal named-port polynomial and independently encoded deck audit."""
from pathlib import Path
from collections import Counter
from itertools import product
from math import gcd,isqrt
import json,hashlib,subprocess,tempfile,shutil,sys,re
P=Path(__file__).resolve().parent

def need(x,m):
    if not x:raise RuntimeError(m)
def opposite(r):return [((1<<len(r))-1)^s^(1<<i)for i,s in enumerate(r)]
def decode(n,z):
    r=[0]*n;b=n*(n-1)//2-1
    for j in range(1,n):
        for i in range(j):
            if z>>b&1:r[i]|=1<<j
            else:r[j]|=1<<i
            b-=1
    return r
def primitive(d):
    d={k:v for k,v in d.items()if v};q=gcd(*d.values());return {k:v//q for k,v in d.items()}
def main():
    src=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else P.parent/'session189-deletion-incoming';manifest=json.loads((src/'RELEASE_FILES.json').read_text());pins={f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in manifest};need(all(pins[f]==v['sha256']for f,v in manifest.items()),'frozen producer pins');meta=json.loads((src/'source_metadata.json').read_text())
    blobs={}
    for f,key in [('source.c','code'),('source.out','output'),('HYP-9455.md','hypothesis')]:
        data=(src/f).read_bytes();h=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest();need(h==meta[key]['sha'],'exact original Git blob');blobs[f]=h
    need(2147483647%2 and all(2147483647%d for d in range(3,isqrt(2147483647)+1,2)),'rank modulus prime')
    for source,exe in [('literal_polynomial.cpp','literal_polynomial'),('canonical.cpp','canonical')]:subprocess.run(['g++','-O3','-std=c++17',str(P/source),'-o',str(P/exe)],check=True)
    input_data=(src/'graphs.txt').read_text();literal_out=subprocess.check_output([str(P/'literal_polynomial')],input=input_data,text=True);canon_out=subprocess.check_output([str(P/'canonical')],input=input_data,text=True);(P/'literal_bank.json').write_text(literal_out);(P/'row_major_cards.json').write_text(canon_out);bank=json.loads(literal_out);cards=json.loads(canon_out);expected=json.loads((src/'bank.json').read_text());oldcards=json.loads((src/'literal_cards.json').read_text());reps=json.loads((src/'representatives.json').read_text());need(len(bank)==len(cards)==len(expected)==len(reps['response_cases'])==12,'all twelve graphs');degrees=[0]*5;codemap={};reversemap={};index={};weighted=[];decks=[];invariants=[]
    for k,(B,C,E,old,R)in enumerate(zip(bank,cards,expected,oldcards,reps['response_cases'])):
        need(B['rows']==E['rows']==R['rows'],'physical rows')
        want=[0]+[1<<i for i in range(9)]+[m for m in range(512)if m.bit_count()==2];need(B['masks']==E['masks']==want,'complete canonical support basis');need(B['coefficients']==E['coefficients']and B['values']==E['values'],'complete literal polynomial');w=len(want)
        for a in want:
            for c in want:degrees[a.bit_count()+c.bit_count()]+=1
        h=sum(B['ends']);need(sum(B['starts'])==h==E['coefficients'][0]==old['H'],'literal old Hamilton endpoints');source_line=R['source_line'].split('|')[0].split();need((C['aut'],C['self_converse'],h)==(int(source_line[2]),int(source_line[1]),int(source_line[3])),'source automorphisms/converse/H');need(C['aut']==old['aut']and C['self_converse']==old['self_converse'],'independent whole canon')
        weights=Counter();deck=Counter()
        for v,(card,oc)in enumerate(zip(C['cards'],old['cards'])):
            need(card['deleted']==oc['deleted']==v and(B['starts'][v],B['ends'][v])==(oc['starts'],oc['ends']),'same deleted physical name/endpoints');x,y=oc['merged'],card['merged'];need(x not in codemap or codemap[x]==y,'old canonical class maps consistently');need(y not in reversemap or reversemap[y]==x,'new canonical map injective');codemap[x]=y;reversemap[y]=x;deck[y]+=1;weights[y]+=B['starts'][v]+B['ends'][v]
        need(sum(weights.values())==2*h,'weighted deck total');weighted.append(primitive(weights));decks.append(deck);lin=sorted(B['coefficients'][i]for i,m in enumerate(want)if m.bit_count()==1);need(len(lin)==9,'linear multiset multiplicity');index[R['node'],R['converse']]=k;invariants.append({'node':R['node'],'converse':R['converse'],'H':h,'linear_c_multiset':lin,'primitive_weighted_deck':sorted(weighted[-1].items())})
    need(len(index)==12,'all node/orientation labels')
    for C in reps['cases']:
        need(decode(9,C['code'])==C['rows'],'decode original census code');need(C['code']==int(C['source_line'].split('|')[-1]),'source line code');i=index[C['node'],False];j=index[C['node'],True];need(bank[i]['rows']==C['rows']and bank[j]['rows']==opposite(C['rows']),'both converse choices exact');need(cards[i]['merged']==cards[j]['merged']and decks[i]==decks[j]and weighted[i]==weighted[j],'merged observer converse invariance')
    comparisons=[]
    for pair in reps['pairs']:
        a=index[pair[0],False];b=index[pair[1],False];need(cards[a]['merged']!=cards[b]['merged'],'distinct whole merged classes');need(decks[a]==decks[b],'independent equal unweighted deck');need(weighted[a]!=weighted[b],'distinct primitive weighted decks')
        for u,v in product([False,True],repeat=2):
            i=index[pair[0],u];j=index[pair[1],v];left=invariants[i]['linear_c_multiset'];right=invariants[j]['linear_c_multiset'];need(Counter(left)!=Counter(right),'all-name-independent exclusion');comparisons.append({'pair':pair,'converse_choices':[u,v],'left_linear_multiset':left,'right_linear_multiset':right,'all_old_name_bijections_excluded':True})
    need(len(comparisons)==12,'all twelve orientation comparisons');source_lines=(src/'source.out').read_text();part=source_lines.split('n=9 ordinary_classes=')[1];reported=[list(map(int,re.findall(r'node (\d+)',line)))for line in part.splitlines()if'merged-deck collision class'in line];need(reported==reps['pairs'],'exact incoming reported pairs')
    (P/'exclusions.json').write_text(json.dumps({'invariants':invariants,'orientation_comparisons':comparisons,'canonical_card_class_map':sorted(codemap.items())},sort_keys=True,indent=2)+'\n')
    # Source reproduction is intentionally labeled as reproduction, not an independent full census.
    with tempfile.TemporaryDirectory(prefix='deletion-source-audit-')as td:
        D=Path(td)
        for f in manifest:shutil.copyfile(src/f,D/f)
        (D/'data').mkdir();subprocess.run(['gcc','-O3','-std=c11',str(D/'source.c'),'-o',str(D/'census')],check=True);replayed=subprocess.check_output([str(D/'census'),'9',str(D/'data')]);need(replayed.decode().strip()==source_lines.split('\nreal\t')[0].strip(),'exact incoming source replay through nine');(P/'census_replay.out').write_bytes(replayed);(D/'reproduced.out').write_bytes(replayed)
        subprocess.check_output([sys.executable,str(D/'prepare.py')]);need((D/'representatives.json').read_bytes()==(src/'representatives.json').read_bytes()and(D/'graphs.txt').read_bytes()==(src/'graphs.txt').read_bytes(),'exact original source representatives')
        for opt,label in [([], 'normal'),(['-O'],'optimized')]:
            out=subprocess.check_output([sys.executable,*opt,str(D/'verify.py')]);need(out==(src/'receipt.json').read_bytes(),'producer normal optimized replay');(P/('producer.'+label+'.json')).write_bytes(out)
    need(pins=={f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in manifest},'frozen producer unchanged')
    receipt={'status':'PASS','literal_response_graphs':12,'literal_named_coefficients':sum(degrees),'degree_counts':degrees,'independent_row_major_card_canonicalizations':108,'orientation_comparisons':12,'all_bijections_excluded_by_linear_multisets':True,'source_census_status':'Exact original-code reproduction through9, not independent census/rank implementation','source_git_blobs':blobs,'producer_pins':pins,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['audit.py','literal_polynomial.cpp','canonical.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['literal_bank.json','row_major_cards.json','exclusions.json','census_replay.out','producer.normal.json','producer.optimized.json']}}
    (P/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
