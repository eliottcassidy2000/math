#!/usr/bin/env python3
"""Complete gap2-ready classification at orders14,15,16.

Generate the inherited unique-cycle/two-entry class grammar, evaluate EVERY
class with an exact C++ subset table, add all five exceptional scalar classes.
The C++ is merely an independent integer counter, not a generation oracle.
"""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import surgery_ready as ready

GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok: raise RuntimeError(label)


def induced(rows,order):
    return tuple(sum(((rows[i]>>j)&1)<<b for b,j in enumerate(order)) for i in order)


def rational_series(numerator,denominator,n):
    out=[]
    for k in range(n+1):
        out.append((numerator[k] if k<len(numerator) else 0)
                   -sum(denominator[j]*out[k-j] for j in range(1,min(k,len(denominator)-1)+1)))
    return out


U=rational_series([0,0,0,1],[1,-3,1],15)
P=rational_series([0,0,0,1,0,-1],[1,-3,1],15)
V=rational_series([0,0,0,0,1,-1,-1,1],[1,-6,11,-6,1],15)


@lru_cache(None)
def unrank_one(n,index):
    """Inherited complete strong unique-start rooted grammar, root0."""
    need(3<=n<=14 and 0<=index<U[n], 'one-start rank domain')
    for m in [1]+list(range(3,n)):
        t=n-m-1
        words=1 if t==0 else 1<<(t-1)
        core_count=1 if m==1 else U[m]
        if index>=words*core_count:
            index-=words*core_count
            continue
        wi,ci=divmod(index,core_count)
        core=(0,) if m==1 else unrank_one(m,ci)
        word=() if not t else tuple((wi>>i)&1 for i in range(t-1))+(1,)
        rows,r=ready.funnel(core,(0,),word)
        return induced(rows,[r]+list(range(r)))
    raise RuntimeError('complete one-start grammar exhausted')


def cyclic_code(rows):
    n=len(rows)
    return min(induced(rows,list(range(i,n))+list(range(i))) for i in range(n))


def generate_bank():
    cycles={};pairs={};two={};records=[];bank=[]
    A,B=0,1
    for n in range(3,15):
        cycles[n]=sorted({cyclic_code(unrank_one(n,i)) for i in range(U[n])})
        if n==3:
            need(len(cycles[n])==1, 'one cyclic triple class')
        else:
            if n>4:A,B=A+B,A+2*B
            need(len(cycles[n])==A+B, 'independent classical two-type ordinary-class count')
        pairs[n]=[]
        for T in cycles[n]:
            low={i for i,row in enumerate(T) if row.bit_count()==1}
            entries=sorted(low|{(i+1)%n for i in low})
            choices=[(0,1)] if n==3 else combinations(entries,2)
            for J in choices:pairs[n].append((T,tuple(J)))
        need(len(pairs[n])==P[n], 'complete two-unique-entry class count')
        d={}
        if n>=4:
            for m in range(3,n):
                t=n-m-1
                suppliers=[(T,J,[tuple(range(j,m))+tuple(range(j)) for j in J])
                           for T,J in pairs[m]]
                suppliers += [(T,(0,),paths) for T,paths in two.get(m,{}).items()]
                for T,J,paths in suppliers:
                    for pre in product((0,1),repeat=max(t-1,0)):
                        word=pre+(1,) if t else ()
                        Q,r=ready.funnel(T,J,word)
                        hp=[(r,)+tuple(path)+tuple(range(m,m+t)) for path in paths]
                        code,path=min((induced(Q,path),path) for path in hp)
                        inverse={v:i for i,v in enumerate(path)}
                        canon_paths=tuple(sorted(tuple(inverse[v] for v in old) for old in hp))
                        need(len(set(canon_paths))==2, 'two distinct actual canonical starting paths')
                        if code in d:need(d[code]==canon_paths, 'duplicate construction retains both actual paths')
                        d[code]=canon_paths
            need(len(d)==V[n], 'complete rooted two-start class count')
            two[n]=d
        bank += [(n,T,J,'pair') for T,J in pairs[n]]
        bank += [(n,T,(0,),'single') for T in sorted(d)]
        records.append({'n':n,'one_start_rooted':U[n],'ordinary_unique_cycle':len(cycles[n]),
                        'two_unique_entry_pairs':len(pairs[n]),'two_start_rooted':len(d)})
    return bank,records


def exceptional_values(n):
    return [3*2**(n-3)-1,6*2**(n-4)+1,5*2**(n-4)-1,
            13*2**(n-5)+3,14*2**(n-5)+3]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    base=Path(__file__).parent
    need(sha256((base/'surgery_ready.py').read_bytes()).hexdigest()
         =='024ac0b483d832865c6f12eabf0fe0f8676ac39188631e44910a0a4c9e66ca2f',
         'frozen inherited word-transfer dependency')
    bank,records=generate_bank()
    need(len(bank)==258900, 'full finite core-entry universe')
    hist={n:Counter() for n in (14,15,16)}
    winners={};target_hits=[];cpp_rows=0;controls=0
    with tempfile.TemporaryDirectory(prefix='surgery-order-proof-') as td:
        td=Path(td);exe=td/'counter';inp=td/'input.txt';out=td/'data.txt'
        compiler=shutil.which('c++') or shutil.which('clang++') or shutil.which('g++')
        need(compiler is not None, 'C++17 compiler available')
        subprocess.run([compiler,'-O3','-std=c++17',str(base/'surgery_orders.cpp'),'-o',str(exe)],check=True)
        with inp.open('w') as f:
            for m,T,J,_ in bank:f.write(' '.join(map(str,(m,sum(1<<j for j in J),*T)))+'\n')
        with inp.open('rb') as f,out.open('wb') as g:
            proc=subprocess.run([str(exe)],stdin=f,stdout=g,stderr=subprocess.PIPE,check=True)
        cpp_summary=json.loads(proc.stderr)
        need(cpp_summary['rows']==len(bank), 'counter processes every declared core')
        for hostile in ('3 3 3 4 1\n','3 1 2 4 1\n'):
            h=subprocess.run([str(exe)],input=hostile.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            need(h.returncode==2 and b'REJECT' in h.stderr, 'invalid adjacency or wrong actual start budget rejected')
        sample_indices={0,len(bank)-1}
        for n in range(3,15):
            ids=[i for i,row in enumerate(bank) if row[0]==n]
            sample_indices.update((ids[0],ids[-1],ids[len(ids)//2]))
        with out.open() as f:
            for index,((m,T,J,kind),line) in enumerate(zip(bank,f)):
                values=list(map(int,line.split()))
                need(len(values)==7, 'complete exact core response row')
                data=dict(zip(('h','p2','e','F','G','b','C'),values));cpp_rows+=1
                if index in sample_indices:
                    need(data==ready.core_data(T,J), 'independent Python core-table reconstruction')
                    controls+=1
                for n in hist:
                    length=n-m-2
                    if length<0:continue
                    for pre in product((0,1),repeat=max(length-1,0)):
                        word=pre+(1,) if length else ()
                        H=ready.ready_count(data,word)
                        hist[n][H]+=1
                        if n not in winners or H<winners[n]['H']:
                            winners[n]={'H':H,'core':T,'ports':J,'word':word,'kind':kind}
                        if H==4289:target_hits.append({'order':n,'core':T,'ports':J,'word':word})
        need(cpp_rows==len(bank), 'no omitted/truncated integer responses')
        input_pin=sha256(inp.read_bytes()).hexdigest();data_pin=sha256(out.read_bytes()).hexdigest()
    result_orders=[]
    for n in (14,15,16):
        need(sum(hist[n].values())==V[n-1], 'complete degree-one ready class universe')
        win=winners[n]
        T,r,a=ready.ready_target(tuple(win['core']),tuple(win['ports']),tuple(win['word']))
        d=ready.counts(T,(r,));s=ready.counts(ready.converse(T),(a,))
        need(len(T)==n and ready.strong(T), 'minimum attained by an actual strong tournament')
        need(d['H']==win['H'] and s['F'][a]-d['H']==2, 'independent literal minimum count and actual gap2')
        need(T[a]==1<<r, 'minimum also satisfies exact joint sole-exit guard')
        exc=exceptional_values(n)
        overall=min(min(hist[n]),min(exc))
        need(overall>4289, 'complete scalar gap2-ready rank excludes target4289')
        ordered=sorted(hist[n].items());raw=(json.dumps(ordered,separators=(',',':'))+'\n').encode()
        result_orders.append({'order':n,'degree_one_classes':sum(hist[n].values()),
                              'exceptional_classes':5,'all_scalar_classes':sum(hist[n].values())+5,
                              'degree_one_minimum':win['H'],'exceptional_values':exc,
                              'scalar_minimum':overall,'distinct_counts':len(hist[n]),
                              'histogram_sha256':sha256(raw).hexdigest(),'lowest_ten':ordered[:10],
                              'minimum_constructor':win,'minimum_rows':T,'minimum_root':r,'minimum_aux':a,
                              'minimum_F':d['F'][r]})
    need(not target_hits, 'no target4289 in any of the three complete ranks')
    # Busch's cited exact minimum at n17 pays the remaining unbounded tail.
    need(9*5**4==5625 and 5625>4289, 'cited strong-order cutoff arithmetic')
    total=GATES+ready.GATES+cpp_summary['gates']
    result={'scope':'complete strong scalar gap2-ready minima at orders14,15,16; cited strong-order bound excludes higher4289 targets',
            'grammar_counts':records,'core_entries':len(bank),'input_sha256':input_pin,
            'core_responses_sha256':data_pin,'python_core_controls':controls,'orders':result_orders,
            'target':4289,'target_hits':target_hits,'gates':total,'python_gates':GATES,
            'shared_counter_gates':ready.GATES,'cpp':cpp_summary}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();args.output.write_bytes(raw)
    print(json.dumps({'gates':total,'core_entries':len(bank),'ready_classes':[r['all_scalar_classes'] for r in result_orders],
                      'minima':[r['scalar_minimum'] for r in result_orders],
                      'receipt_bytes':len(raw),'receipt_sha256':sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__':main()
