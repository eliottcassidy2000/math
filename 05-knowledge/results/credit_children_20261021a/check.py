"""Ground every declared scale-credit source and paid child by supplied words."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import spec_from_file_location,module_from_spec
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('grounded_native',HERE.parent/'session18_native_annulus/verify_annulus.py')
native=module_from_spec(spec);sys.modules[spec.name]=native;spec.loader.exec_module(native)
CHECKS=0

def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)

def v2(x):return (x&-x).bit_length()-1

def replay(n,w):
    need(type(n) is int and n>0 and n%2,'literal positive odd source')
    x=n;D=F(0);xs=[n];ds=[D]
    for a in w:
        need(x!=1,'strict first ROOT stop')
        need(type(a) is int and a>=1 and a==v2(3*x+1),'actual valuation letter')
        x=(3*x+1)>>a
        D=F(3,2**a)*D+max(F(0),1-F(4,2**a))
        xs.append(x);ds.append(D)
    need(x==1,'supplied word reaches ROOT')
    z=1
    for a in reversed(w):
        q=(z<<a)-1
        need(q>0 and q%3==0,'independent inverse integrality')
        z=q//3
    need(z==n,'independent source reconstruction')
    return xs,ds

def inverse(u,a,j):
    need(type(u) is int and native.source(u),'primitive parent')
    need(type(a) is int and a in (1,2),'literal primitive branch')
    need(type(j) is int and j>=0,'literal strip count')
    q=(1<<a)*(3*u+1)*4**j-((1<<a)+3)
    need(q%9==0,'inverse branch integrality')
    s=q//9
    need(native.source(s) and native.edge(s)==(u,a+2*j),'actual primitive edge and clock')
    return s

def compile_child(u,a,j,w):
    replay(u,w)  # authenticate the parent before replacing its first letter
    s=inverse(u,a,j)
    need(s>u and s>1,'original source exceeds parent')
    b=v2(3*u+1)
    tail=list(w[1:]) if w else []
    if j:
        z=(4**j*(3*u+1)-1)//3
        word=[a,b+2*j]+tail
    else:
        z=u;word=[a]+list(w)
    xs,ds=replay(s,word)
    need(xs[1]==z and word[0]==a,'literal first paid child')
    stop=1 if a==2 else 2
    need(stop<len(xs) and ds[stop]<s-xs[stop],'strict original-source payment')
    child=xs[stop]
    replay(child,word[stop:])
    return {'source':s,'mod15':s%15,'mod35':s%35,'parent':u,'branch':a,'strips':j,
            'native_edge_clock':a+2*j,'payment_steps':stop,'paid_child':child,
            'paid_child_mod15':child%15,'paid_child_mod35':child%35,
            'payment_margin':str(s-child-ds[stop]),'payment_clock':sum(word[:stop]),
            'ROOT_word':word,'ROOT_clock':sum(word),'child_ROOT_word':word[stop:]}

def rejects(fn):
    try:fn()
    except ValueError:return True
    return False

def run():
    raw=(HERE.parent/'scale_credit_supply_20261021a/receipt.json').read_bytes().replace(b'\r\n',b'\n')
    need(sha256(raw).hexdigest()=='5b9b010d19c5b1d38f380d1e31567d44750f5277e786a5d6b6c446a7785ff64a','frozen source-owned credit bank')
    data=json.loads(raw)['case8192'];owners=data['linear_rows']+data['frontier_supplement']
    inp=json.loads((HERE/'inputs.json').read_text());words={}
    for row in inp['parents']:
        u=row['parent']
        need(type(row['mod15']) is int and row['mod15']==u%15,'retained parent mod15')
        need(type(row['mod35']) is int and row['mod35']==u%35,'retained parent mod35')
        need(u not in words,'one supplied word per parent')
        replay(u,row['ROOT_word']);words[u]=row['ROOT_word']
    need(set(words)=={r['parent'] for r in owners},'exact declared parent universe')
    need(len(owners)==380 and len({r['source'] for r in owners})==380,'380 distinct source owners')
    rows=[];counts=Counter();largest=0
    for owner in owners:
        u,a,j=owner['parent'],owner['branch'],owner['strips']
        result=compile_child(u,a,j,words[u])
        need(result['source']==owner['source'],'original credit source retained')
        need(result['native_edge_clock']==owner['child_clock']-owner['parent_clock'],'native receipt clock retained')
        need(owner['source_row'][1:3]==[result['mod15'],result['mod35']],'both credit source registers retained')
        need(result['mod15']%5==result['mod35']%5,'compatible source registers')
        result['family']=owner['family'];result['terminal']=owner['terminal']
        result['arrival_clock']=owner['child_clock'];result['arrival_depth']=owner['child_depth']
        counts[str(result['payment_steps'])]+=1;rows.append(result)
        for k in (1,2,10):
            lifted=compile_child(u,a,j+9*k,words[u]);F64=64**(3*k)
            need(lifted['source']==F64*result['source']+((1<<a)+3)*(F64-1)//9,'exact three-scale affine lift')
            need(lifted['native_edge_clock']==result['native_edge_clock']+18*k,'native integer lift clock')
            need(4*8192*F64<=lifted['source']<256*8192*F64,'source annulus preserved')
            need(lifted['source']%3==result['source']%3,'charge phase restored')
            largest=max(largest,lifted['source'].bit_length())
    need(rejects(lambda:compile_child(3,True,8,words[3])),'Boolean branch rejected')
    need(rejects(lambda:compile_child(3,1,0,words[3])),'ungrounded prefix outside source-greater guard rejected')
    bad=list(words[3]);bad[0]+=1
    need(rejects(lambda:replay(3,bad)),'false supplied valuation rejected')
    need(rejects(lambda:compile_child(3,1,8,[999,4])),'discarded parent letter must authenticate')
    need(rejects(lambda:compile_child(1,1,1,[2])),'ROOT-padded parent word rejected')
    return {'status':'FINITE-EXACT 380 base owners; PROVED inherited inverse-word compiler extends each guarded ray',
            'parent_count':len(words),'source_count':len(rows),'payment_step_counts':dict(counts),
            'largest_probe_bits':largest,'scope':'ROOT grounding does not prove timely arrival counts or uniform153N',
            'input_sha256_lf':sha256((HERE/'inputs.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
            'rows':rows,'checks':CHECKS}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=run();raw=(json.dumps(result,indent=2)+'\n').encode();path=HERE/'receipt.json'
    if a.write:path.write_bytes(raw)
    else:need(path.read_bytes().replace(b'\r\n',b'\n')==raw,'immutable receipt comparison')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
