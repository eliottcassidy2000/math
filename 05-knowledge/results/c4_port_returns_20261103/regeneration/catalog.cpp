#define main small_bank_original_main
#include "small_bank.cpp"
#undef main

// Every n-vertex isomorphism class has an (n-1)-vertex deletion.  All birth
// cuts from one representative of every deletion class therefore suffice.
map<unsigned,T> extend_bank(const map<unsigned,T>& bank) {
    map<unsigned,T> next;
    for (const auto& [code,t]:bank)
        for(unsigned cut=0;cut<(1u<<t.size());cut++) {
            auto child=birth(t,cut);
            next.emplace(canonical(child),child);
        }
    return next;
}

int main(int argc,char**argv) {
    const int max_n=argc>1?stoi(argv[1]):7;
    const I target_f=argc>2?stoll(argv[2]):63;
    const vector<int> targets=argc>3?vector<int>{stoi(argv[3])}:vector<int>{-35,-29};
    if(max_n<1||max_n>8||target_f<1)return 2;
    vector<map<unsigned,T>> banks(max_n+1);
    banks[1].emplace(0,T{0});
    for(int n=2;n<=max_n;n++)banks[n]=extend_bank(banks[n-1]);
    cout<<"{\"max_order\":"<<max_n<<",\"target_f\":"<<target_f<<",\"counts\":[";
    for(int n=1;n<=max_n;n++){if(n>1)cout<<',';cout<<banks[n].size();}
    cout<<"],\"target_intercepts\":[";
    for(int i=0;i<(int)targets.size();i++){if(i)cout<<',';cout<<targets[i];}
    cout<<"],\"representatives\":[";
    int printed=0;
    struct Hit {T t;Profile p;unsigned W;int c;};
    vector<Hit>hits;
    for(int n=1;n<=max_n;n++)for(const auto&[code,t]:banks[n]){
        auto p=profile(t);auto term=terminal(t);
        if(printed++)cout<<',';
        cout<<"{\"n\":"<<n<<",\"code\":"<<code<<",\"h\":"<<p.h
            <<",\"p2\":"<<p.p<<",\"strong\":"<<(term==(1u<<n)-1?"true":"false")
            <<",\"rows\":[";
        for(int i=0;i<n;i++){if(i)cout<<',';cout<<t[i];}
        cout<<"]}";
        for(unsigned W=1;W<(1u<<n);W++){
            if(!(W&term))continue;
            I f=0,j=0;
            for(int u=0;u<n;u++)if(W>>u&1){
                f+=p.f[u];
                for(int v=0;v<n;v++)if(W>>v&1)j+=p.j[u][v];
            }
            if(f==target_f&&find(targets.begin(),targets.end(),p.h-j)!=targets.end())
                hits.push_back({t,p,W,int(p.h-j)});
        }
    }
    cout<<"],\"hits\":[";
    printed=0;
    for(const auto&hit:hits){if(printed++)cout<<',';print_hit(hit.t,hit.p,hit.W,hit.c,0,-1);}
    cout<<"]}\n";
}
