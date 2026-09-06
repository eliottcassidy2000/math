#!/usr/bin/env python3
"""Independent consumer audit of the pinned full-word cost maxima.

Reuses only the pinned independent baseline AUDITOR's literal interval
geometry and bounded-knapsack helper engine; neither producer is imported.
The maximum-cost supplier is separately proved/audited.
"""
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path
import importlib.util
import json
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
BASE_AUDIT_SHA='ec841b2a3e9898dcf675d4c0212b78e3d16e21b8d13a05aa54b535a231226cc0'
BASE_SET_SHA='a25d83f0eeb630bb82e84cdfac4e3cf7312f892f6c426d6affd5239a064e4b58'
MAXIMA_SHA='ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca'
GATES=0
DIGEST=sha256()


def gate(ok,label,data=None):
    global GATES
    if not ok:raise RuntimeError(f'{label}: {data!r}')
    GATES+=1;DIGEST.update((label+':'+repr(data)+'\n').encode())


def canonical(value):return json.dumps(value,separators=(',',':'))

def digest(value):return sha256(canonical(value).encode()).hexdigest()


REFINED_MAIN=r'''
int main(){
    int n;cin>>n;vector<int> allowed(n);for(int&d:allowed)cin>>d;
    int nr;cin>>nr;vector<Row> rows(nr);
    for(Row&r:rows){int k;cin>>r.p>>r.q>>r.J>>r.D>>r.M>>k;
        for(int i=0;i<k;++i){int a,b;cin>>a>>b;r.lengths.emplace_back(a,b);}}
    map<int,pair<int,int>> slopes;
    for(const Row&r:rows){auto it=slopes.find(r.J);
        if(it==slopes.end() || ll(r.M)*it->second.second<ll(it->second.first)*r.D)
            slopes[r.J]={r.M,r.D};}
    const int caps[7]={90,30,9,4,2,1,1};
    int clocks;cin>>clocks;
    vector<int> aggregate,components,coupled;
    for(int index=0;index<clocks;++index){
        int t,E;cin>>t>>E;
        vector<int> values;array<int,91> capacity{},cost{};
        for(int d:allowed)if(t%d==0){
            values.push_back(d);for(int cap:caps)capacity[d]+=(d<=cap);
            cost[d]=7*d*up(t,7*d)-t;
            need(cost[d]==d*((7-(t/d)%7)%7),"literal marginal cost");
        }
        int score=knapsack(values,capacity,cost,7);
        need(score>=0 && score%7==0,"unconditioned seven-slot maximum");
        int old_E=score/7;
        need(E>=0 && E<=old_E,"full-word maximum strengthens the bag");
        vector<int> es;
        for(int e=1;e<=6;++e)if(t%e==0){
            ll least=1000000000;
            for(auto item:slopes){int J=item.first,M=item.second.first,D=item.second.second;
                least=min(least,ll(e)*(up(ll(t/e)*M,D)-J));}
            if(least<=E)es.push_back(e);
        }
        if(es.empty())continue;
        aggregate.push_back(t);bool has_component=false,has_coupled=false;
        map<pair<int,int>,int> conditioned;
        for(int e:es){
            int order=t/e;
            for(const Row&r:rows){
                ll merged=ll(e)*(up(ll(order)*r.M,r.D)-r.J);
                if(merged>E)continue;
                ll credit=-ll(e)*r.J;
                for(auto item:r.lengths)credit+=ll(e)*item.second*up(ll(order)*item.first,r.D);
                need(credit>=max(0LL,merged),"nonnegative full component credit");
                if(credit>E)continue;
                has_component=true;
                int dp=e*gcd(order,r.p),dq=e*gcd(order,r.q);
                need(gcd(dp,dq)==e,"forced states preserve the sheet gcd");
                if(dp>90 || dq>90 || !capacity[dp] || !capacity[dq])continue;
                if(dp>dq)swap(dp,dq);
                auto key=make_pair(dp,dq);auto it=conditioned.find(key);int EE;
                if(it==conditioned.end()){
                    int rest=knapsack(values,capacity,cost,5,dp,dq);
                    if(rest<0)EE=-1;
                    else{int s=rest+cost[dp]+cost[dq];need(s%7==0,"reserved-pair seven-slot integrality");EE=s/7;}
                    conditioned[key]=EE;
                }else EE=it->second;
                need(EE<=old_E,"conditioned bag ceiling");
                if(EE>=0 && credit<=min(E,EE)){has_coupled=true;break;}
            }
            if(has_coupled)break;
        }
        if(has_component)components.push_back(t);
        if(has_coupled)coupled.push_back(t);
        need(!has_coupled || has_component,"finite stage nesting");
    }
    emit("full_words_aggregate",aggregate);
    emit("full_words_components",components);
    emit("full_words_coupled",coupled);
    cout<<"CPP_GATES "<<gates<<'\n';
}
'''


def main():
    source=ROOT/'04-computation/third_20260906_grid_bootstrap_audit.py'
    gate(sha256(source.read_bytes()).hexdigest()==BASE_AUDIT_SHA,'pinned-independent-audit-module')
    spec=importlib.util.spec_from_file_location('independent_grid_audit_helpers',source)
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    baseline=None
    for line in (ROOT/'05-knowledge/results/third_20260906_grid_bootstrap_audit.out').read_text().splitlines():
        if line.startswith('AUDIT_SCALE_SET '):
            item=json.loads(line[len('AUDIT_SCALE_SET '):])
            if item['name']=='profile6_coupled':baseline=item['survivors']
    gate(baseline is not None and len(baseline)==8301 and digest(baseline)==BASE_SET_SHA,'complete-audited-baseline-domain')
    lines=(ROOT/'05-knowledge/results/third_20260906_grid_full_words.out').read_text().splitlines()
    maxima=[json.loads(line[len('MAXIMUM '):]) for line in lines if line.startswith('MAXIMUM ')]
    normalized=[[row[0],row[1]] for row in maxima]
    gate([row[0] for row in maxima]==baseline and digest(normalized)==MAXIMA_SHA,'pinned-complete-full-cost-table')
    gate('MAXIMA_SHA256 '+MAXIMA_SHA in lines and any(line.startswith('PASS gates=') for line in lines),'full-cost-supplier-completed')
    for t,E,owner,visits in maxima:
        gate(len(owner)==7 and gcd(*owner)==1 and all(d>0 and t%d==0 for d in owner),'table-owner-type',t)
        gate(sum(d*((t+7*d-1)//(7*d)) for d in owner)-t==E,'table-owner-literal-cost',t)
    raw=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    gate(sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','pinned-profile-supplier')
    levels=json.loads(raw)['levels']
    profiles={int(k):{(c,tuple(w)) for c,w in level['profiles']} for k,level in levels.items()}
    allowed=sorted({c for c,w in profiles[6]})
    rows=audit.make_atlas();lookup={r[:2]:r for r in rows}
    data=[str(len(allowed)),' '.join(map(str,allowed)),str(len(rows))]
    for p,q,J,D,lengths,M in rows:
        data.append(' '.join(map(str,(p,q,J,D,M,len(lengths)))))
        data.extend(f'{n} {mult}' for n,mult in lengths)
    data.append(str(len(normalized)));data.extend(f'{t} {E}' for t,E in normalized)
    gate(audit.CPP.count('int main(){')==1,'independent-engine-helper-boundary')
    engine=audit.CPP.split('int main(){')[0]+REFINED_MAIN
    with tempfile.TemporaryDirectory(prefix='third-grid-refined-audit-') as directory:
        directory=Path(directory);cpp=directory/'engine.cpp';binary=directory/'engine'
        cpp.write_text(engine)
        subprocess.run(['g++','-std=c++17','-O2',str(cpp),'-o',str(binary)],check=True,capture_output=True,text=True)
        output=subprocess.run([str(binary)],input='\n'.join(data)+'\n',check=True,capture_output=True,text=True).stdout
    recovered={};cpp_gates=0
    for line in output.splitlines():
        name,payload=line.split(' ',1)
        if name=='CPP_GATES':cpp_gates=int(payload)
        else:recovered[name]=json.loads(payload)
    producer={}
    for line in (ROOT/'05-knowledge/results/third_20260906_grid_refined.out').read_text().splitlines():
        if line.startswith('SCALE_SET '):
            item=json.loads(line[10:]);producer[item['name']]=item['survivors']
    expected={'full_words_aggregate':(8288,21600),'full_words_components':(8202,16704),'full_words_coupled':(8202,16704)}
    for name in sorted(expected):
        values=recovered[name]
        gate(values==sorted(set(values)) and set(values)<=set(baseline),'canonical-complete-refined-subset',name)
        gate((len(values),max(values))==expected[name],'complete-refined-census',name)
        gate(values==producer[name],'independent-refined-array-comparison',name)
        print('AUDIT_SCALE_SET',canonical({'name':name,'count':len(values),'maximum':max(values),'sha256':digest(values),'survivors':values}))
    survivors=recovered['full_words_coupled']
    gate(survivors==recovered['full_words_components'],'component-and-conditioned-arrays-equal')
    removed=sorted(set(baseline)-set(survivors))
    gate(len(removed)==99 and [t for t in survivors if t>14904]==[16704],'99-deletions-and-isolated-upper-clock')
    # Positive compatibility witness in the FULL profile quotient, not a
    # realized physical decoder. This independently checks all126 words.
    t,e,p,q=16704,4,3,308
    owner=(12,16,72,58,64,9,9)
    gate(gcd(*owner)==1 and all(t%d==0 for d in owner),'boundary-owner-arithmetic')
    for size in range(1,7):
        for I in combinations(range(7),size):
            c=gcd(*(owner[i] for i in I));word=tuple(sorted(gcd(c,owner[j]) for j in range(7) if j not in I))
            gate((c,word) in profiles[7-size],'boundary-full-complement-word',(I,c,word))
    gate((e*gcd(t//e,p),e*gcd(t//e,q))==owner[:2],'boundary-forced-pair-compatibility')
    cost=sum(d*((t+7*d-1)//(7*d)) for d in owner)-t
    credit=audit.component_credit(t,e,lookup[(p,q)])
    gate(dict(normalized)[t]==cost==188 and credit==172 and credit<cost,'boundary-sharp-stopping-control')
    print('BOUNDARY',canonical([t,e,p,q,list(owner),cost,credit]))
    print('REMOVED_FROM_BASELINE',canonical(removed))
    print('MAXIMA_SHA256',MAXIMA_SHA)
    print('CPP_EXACT_GATES',cpp_gates)
    print('PYTHON_EXACT_GATES',GATES)
    print('INHERITED_AUDIT_GEOMETRY_GATES',audit.GATES)
    print('SEMANTIC_SHA256',DIGEST.hexdigest())
    print('PASS')


if __name__=='__main__':main()
