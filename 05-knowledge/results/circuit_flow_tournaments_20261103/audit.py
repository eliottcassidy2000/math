"""Separate graph replay and inherited hostile controls; same-agent implementation."""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import hashlib,importlib.util,json,sys

HERE=Path(__file__).resolve().parent


def need(ok,why):
    if not ok:raise RuntimeError(why)


def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod


def mat(a):return [[bool(row>>j&1) for j in range(len(a))] for row in a]
def hp(a):
    n=len(a)
    @lru_cache(None)
    def visit(left,v):
        if not left:return 1
        return sum(visit(left^(1<<u),u) for u in range(n) if left>>u&1 and a[v][u])
    return sum(visit(((1<<n)-1)^(1<<v),v) for v in range(n)) if n else 1


def birth(a,cut):
    n=len(a);b=[list(row)+[bool(cut>>i&1)] for i,row in enumerate(a)]
    b.append([not bool(cut>>i&1) for i in range(n)]+[False]);return b


def native(a,path):
    n=len(a);need(sorted(path)==list(range(n)),'independent spine')
    need(all(a[u][v] for u,v in zip(path,path[1:])),'independent path legality')
    fixed={tuple(sorted((u,v))) for u,v in zip(path,path[1:])}
    return [[False if i==j else a[i][j] if tuple(sorted((i,j))) in fixed else not a[i][j] for j in range(n)] for i in range(n)]


def hostile(reverse):
    a=[[False]*12 for _ in range(12)]
    for i,j in combinations(range(12),2):
        b,s=divmod(i,4);c,t=divmod(j,4)
        arc=(c-b)%3==1 if b!=c else s==0 or (t-s)%3==1
        if reverse and b==c==2:arc=not arc
        a[i][j]=arc;a[j][i]=not arc
    return tuple(sum(int(z)<<j for j,z in enumerate(row)) for row in a)


def main():
    d=json.loads((HERE/'certificate.json').read_text());counts=[]
    for x in d['neutral_square_and_readers']:
        a=[[False]]
        for cut in x['birth_cuts']:a=birth(a,cut)
        need(a==mat(x['rows']) and hp(a)==x['H'],'independent literal birth and H count')
        counts.append(x['H'])
    q=json.loads((HERE/'native_reset.json').read_text());a=mat(q['source']['rows']);n=len(a)
    for step in q['trace']:
        for cut in step['cuts']:a=birth(a,cut)
        for p in step['paths']:a=native(a,p)
        a=[row[:n] for row in a[:n]]
        need(a==mat(step['completed_rows']) and hp(a)==step['completed_H'],'independent native event replay')
    need(a==mat(q['target']['rows']),'actual target')
    sym=load_module('symmetry_control',HERE.parent/'regular_parent_symmetry_20261029/verify.py')
    junc=load_module('junction_control',HERE.parent/'coverage_machine_induction_20261103/junction.py')
    A,B=hostile(False),hostile(True)
    triples=lambda a:[t for t in combinations(range(12),3) if all(sum(a[u]>>v&1 for v in t)==1 for u in t)]
    need(triples(A)==triples(B) and len(triples(A))==67,'actual labelled triple equality')
    autos=[len(sym.iso_maps(a,a)) for a in (A,B)];need(autos==[81,27],'actual full automorphism orders')
    f=[junc.slot_values(a) for a in (A,B)];need(f[0]==f[1] and f[0][0]==462915,'full degree12 profile equality')
    read=[hp(birth(mat(a),291)) for a in (A,B)];need(read[0]!=read[1],'named nonuniform read distinguishes sources')
    profiles=[];suffix_reads=[]
    for x in d['neutral_square_and_readers'][:4]:
        a=tuple(x['rows']);profiles.append({'word':x['word_inner_to_outer'],
                  'root_profiles':[junc.profile(a,i) for i in range(len(a))]})
        if x['word_inner_to_outer'] in ('UL','LU'):
            b=birth(birth(mat(a),1),1<<len(a));z=hp(b)
            need(z==junc.H(a)+2*junc.profile(a,0)[3],'actual rooted two-birth read')
            suffix_reads.append(z)
        junc.endpoints.cache_clear()
    need(suffix_reads==[6467,6379],'named-root correction44, child correction88')
    inherited=HERE.parent/'second_feedback_18979_20261103/source_pins.json'
    out={'status':'FINITE_EXACT same-agent separate implementations','reader_H':counts,'native_reset_replayed':True,
         'all_nine_named_root_profiles':profiles,'root0_two_birth_H':suffix_reads,
         'hostile':{'rows':[A,B],'automorphism_orders':autos,'same_labelled_triples':67,
                    'common_slot_values_1_to_12':f[0],'nonuniform_cut':291,'read_H':read},
         'inherited_source_pins_sha256':hashlib.sha256(inherited.read_bytes()).hexdigest(),
         'scope':'independent counting/native implementations; exact automorphisms and full profiles; source pins retain failed55->111 and prior source returns'}
    (HERE/'audit.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
    print('PASS separate native/H audit; actual81/27 profiles equal; cut291 reads',read,flush=True)


if __name__=='__main__':main()
