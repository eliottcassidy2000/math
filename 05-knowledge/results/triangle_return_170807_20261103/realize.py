"""Compile the discovered two-edge interface into actual named complements."""
from pathlib import Path
from functools import lru_cache
from itertools import permutations
import json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'return6325_completion_20261103'))
from audit import need,matrix,birth,native,strong
def guard(T,a,b,ports):
    full=(1<<len(T))-1
    @lru_cache(None)
    def walk(mask,last):
        if mask==full:return () if last not in (b,*ports) else None
        avail=(T[last]&~mask) if last>=0 else full
        if last==a:avail&=1<<b
        else:avail&=~(1<<b)
        while avail:
            bit=avail&-avail;avail-=bit;v=bit.bit_length()-1;tail=walk(mask|bit,v)
            if tail is not None:return (v,)+tail
        return None
    return walk(0,-1)
def main():
    old=next(z for z in json.loads((HERE.parent/'return18979_palindrome_20261103/returns.json').read_text()) if z['next_H']==56935)
    found=json.loads((HERE/'two_edges_birth.json').read_text());need(found['found'],'actual discovery required')
    source=tuple(old['next_rows']);names=old['next_names'];r=old['next_root'];sp=names.index(6);full=(1<<len(source))-1
    for order in permutations(found['arcs']):
        T=source;trace=[]
        for step,(a,b) in enumerate(order):
            if not (T[a]>>b)&1:a,b=b,a
            P=guard(T,a,b,(r,sp))
            if P is None:break
            helper=52+step;N=len(T);cut=full^(1<<a);X=birth(T,cut);i=P.index(a)
            P1=(*P,N);P2=(*P[:i+1],N,*P[i+1:]);Y=native(X,P1);V=native(Y,P2);U=tuple(x&full for x in V[:N])
            expected=list(T);expected[a]^=1<<b;expected[b]^=1<<a
            need(U==tuple(expected),'isolated named edge reversal')
            if not strong(U):break
            def frame(path,port):
                j=path.index(port);return (path[j-1] if j else None,path[j+1] if j+1<len(path) else None)
            need(all(frame(P1,p)==frame(P2,p) for p in (r,sp)),'common physical joint frame')
            local=names+[helper];trace.append(dict(arc=[names[a],names[b]],helper=helper,cut=cut,paths=[[local[v] for v in path] for path in (P1,P2)],rows=U));T=U
        if len(trace)==2:break
    need(len(trace)==2,'strong guarded order found')
    target=birth(T,found['cut']);need(target==tuple(found['rows']) and strong(target),'literal strong endpoint')
    record=dict(status='FINITE_EXACT_NATIVE_INTERFACE',source_H=56935,target_H=170807,source=source,source_names=names,source_root=r,source_spare=sp,
        trace=trace,native_complements=4,temporary_helpers=[52,53],deleted_old_names=[],birth_cut=found['cut'],birth_name=54,target=target,target_names=names+[54],target_root=r,target_spare=sp,
        scope='Count-changing repair then one birth; not a ready reset or closed ternary grammar.')
    (HERE/'interface.json').write_text(json.dumps(record,indent=2)+'\n');print('PASS four actual complements, common root/spare frames, all old names retained')
if __name__=='__main__':main()
