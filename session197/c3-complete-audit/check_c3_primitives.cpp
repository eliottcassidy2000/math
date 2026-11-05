#define main census_main_not_called
#include "check_c3_complete.cpp"
#undef main
#include <random>
int main(){mt19937_64 rng(197312);long long packed=0,graphs=0,cuts=0,relabellings=0;for(int z=0;z<100000;z++){uint64_t a=0,b=0,want=0;for(int j=0;j<12;j++){int x=rng()%5,y=rng()%5;a|=uint64_t(x)<<(4*j);b|=uint64_t(y)<<(4*j);want|=uint64_t(min(4,x+y))<<(4*j);}need(capadd(a,b)==want,"independent lane arithmetic");packed++;}
 vector<map<string,string>>canonbanks(7);
 for(int n=1;n<=6;n++){vector<pair<int,int>>pairs;for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)pairs.push_back({u,v});for(int code=0;code<(1<<(int)pairs.size());code++){vector<int>a(n);for(int b=0;b<(int)pairs.size();b++){auto[u,v]=pairs[b];if(code>>b&1)a[u]|=1<<v;else a[v]|=1<<u;}auto M=endpoints(a);vector<vector<int>>literal(n,vector<int>(n));vector<int>p(n);iota(p.begin(),p.end(),0);string ordinary;do{bool path=true;for(int j=0;j+1<n;j++)path&=a[p[j]]>>p[j+1]&1;if(path)literal[p.front()][p.back()]++;if(n<=5){string s(2*n,'\0');for(int u=0;u<n;u++){int r=0;for(int v=0;v<n;v++)if(a[p[u]]>>p[v]&1)r|=1<<v;s[2*u]=r&255;s[2*u+1]=r>>8;}if(ordinary.empty()||s<ordinary)ordinary=s;}}while(next_permutation(p.begin(),p.end()));for(int u=0;u<n;u++)for(int v=0;v<n;v++)need(M[u][v]==min(4,literal[u][v]),"literal all-start path counts");graphs++;
 if(n<=5){string key=canon(a);auto it=canonbanks[n].find(ordinary);if(it==canonbanks[n].end())canonbanks[n][ordinary]=key;else need(it->second==key,"full-permutation isomorphism invariant");for(int mask=1;mask<(1<<n)-1;mask++){int exact=0,capped=0;for(int u=0;u<n;u++)if(mask>>u&1)for(int v=0;v<n;v++)if(!(mask>>v&1)){exact+=literal[u][v];capped+=M[u][v];}need(min(4,exact)==min(4,capped),"capped cut threshold");cuts++;}}
 }}
 for(int n=1;n<=5;n++){map<string,int>seen;for(auto&kv:canonbanks[n])need(seen.emplace(kv.second,1).second,"refinement canonical false isomorphism");}
 cout<<"PASS packed_lane_pairs "<<packed<<" all_labelled_endpoint_graphs_through6 "<<graphs<<" literal_cut_controls "<<cuts<<" full_permutation_canonical_classes";for(int n=1;n<=5;n++)cout<<" "<<canonbanks[n].size();cout<<"\n";
}
