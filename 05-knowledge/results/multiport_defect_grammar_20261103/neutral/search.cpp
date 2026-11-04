// Exact local test at the failed inner collar. Inherit only the tested
// subset-counting routines; this main has a different target and filter.
#include <algorithm>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main inherited_search_main
#include "../../zero_ready_strong_switch_20261103/literal/search.cpp"
#undef main
#pragma GCC diagnostic pop

vector<int> zero_by_components(const T& t) {
    unsigned full=(1u<<t.size())-1;
    vector<int> answer;
    for(int r=0;r<(int)t.size();++r) {
        unsigned incoming=full^t[r]^(1u<<r);
        if(popcount(incoming)!=1) continue;
        unsigned allowed=full^(1u<<r), reached=incoming, previous=0;
        while(reached!=previous) {
            previous=reached;
            for(unsigned q=reached;q;q&=q-1)reached|=t[countr_zero(q)]&allowed;
        }
        // In a strong tournament the sole inport fails to start a spanning
        // path after deleting r exactly when it cannot reach that deletion.
        if(reached!=allowed) answer.push_back(r);
    }
    return answer;
}

int main() {
    int n, oldmark, minimum, maximum_edits, forced_a, forced_b;
    I target;
    if (!(cin >> n >> oldmark >> minimum >> maximum_edits >> forced_a >> forced_b >> target)) return 2;
    if (n != 14 || minimum < 1 || maximum_edits < 0 || maximum_edits > 3) return 2;
    T source(n); for (auto &row : source) cin >> row;
    const unsigned full=(1u<<n)-1;
    for (int a=0;a<n;++a) {
        if(source[a]>full || (source[a]>>a&1)) return 3;
        for(int b=0;b<a;++b) if(((source[a]>>b&1)+(source[b]>>a&1))!=1) return 3;
    }
    vector<pair<int,int>> edges;
    int forced=-1;
    for(int a=0;a<n;++a)for(int b=a+1;b<n;++b) {
        if(a==forced_a && b==forced_b) forced=edges.size();
        edges.push_back({a,b});
    }
    if(forced_a>=0 && forced<0) return 4;
    int examined=0, strongs=0, repairable=0, hits=0;
    I nearest=numeric_limits<I>::max();
    cout<<"{\"target_H\":"<<target<<",\"minimum_prefix\":"<<minimum
        <<",\"maximum_edits\":"<<maximum_edits<<",\"forced_edge\":["
        <<forced_a<<','<<forced_b<<"],\"hits\":[";
    auto examine=[&](const vector<int>& edit) {
        if(forced>=0 && find(edit.begin(),edit.end(),forced)==edit.end()) return;
        ++examined; T t=source;
        for(int e:edit) {auto[a,b]=edges[e];t[a]^=1u<<b;t[b]^=1u<<a;}
        if(!strong(t)) return; ++strongs;
        auto component_zero=zero_by_components(t);
        if(component_zero.empty()) return;
        auto q=profile(t,oldmark);
        vector<int> zero;
        for(int r=0;r<n;++r) if(q.all_f[r]==q.h) zero.push_back(r);
        if(zero!=component_zero) throw runtime_error("component and endpoint zero-ready tests disagree");
        if(zero.empty()) return; ++repairable;
        for(int w=0;w<n;++w)for(int s=minimum;s<62;++s) {
            __int128 value=q.h+(__int128(1)<<s)*q.all_f[w];
            if(value>numeric_limits<I>::max()) break;
            I delta=I(value)-target;
            nearest=min(nearest,delta<0?-delta:delta);
            if(value==target) {
                if(hits++) cout<<',';
                cout<<"{\"prefix_size\":"<<s<<",\"mark\":"<<w
                    <<",\"h\":"<<q.h<<",\"F\":"<<q.all_f[w]<<",\"zero_roots\":[";
                for(int i=0;i<(int)zero.size();++i){if(i)cout<<',';cout<<zero[i];}
                cout<<"],\"edges\":[";
                for(int i=0;i<(int)edit.size();++i){if(i)cout<<',';auto[a,b]=edges[edit[i]];cout<<'['<<a<<','<<b<<']';}
                cout<<"],\"rows\":";rows(t);cout<<'}'<<flush;
            }
            if(value>target) break;
        }
    };
    examine({});
    if(maximum_edits>=1) for(int i=0;i<(int)edges.size();++i)examine({i});
    if(maximum_edits>=2) for(int i=0;i<(int)edges.size();++i)for(int j=0;j<i;++j)examine({j,i});
    if(maximum_edits>=3) for(int i=0;i<(int)edges.size();++i)for(int j=0;j<i;++j)for(int k=0;k<j;++k)examine({k,j,i});
    cout<<"],\"examined\":"<<examined<<",\"strong\":"<<strongs
        <<",\"with_zero_ready_mark\":"<<repairable<<",\"nearest_distance\":"<<nearest<<"}\n";
}
