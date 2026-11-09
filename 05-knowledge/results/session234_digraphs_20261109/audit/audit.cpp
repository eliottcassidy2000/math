#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>
using namespace std;
int n,base,m; vector<uint16_t> counts; vector<int> powers,required; int idx[8][8]; vector<pair<int,int>> pairs;
void scatter(int k,int code){if(k==m){counts[code]++;return;} if(required[k]<0){for(int d=0;d<base;d++)scatter(k+1,code+d*powers[k]);}else {scatter(k+1,code+required[k]*powers[k]);if(base==3)scatter(k+1,code+2*powers[k]);}}
bool strong(const array<unsigned,8>& a){for(int s=0;s<n;s++){unsigned seen=1u<<s,old;do{old=seen;for(int i=0;i<n;i++)if(seen>>i&1)seen|=a[i];}while(old!=seen);if(seen!=((1u<<n)-1))return false;}return true;}
uint64_t canon(const array<unsigned,8>& a){vector<int> p(n);iota(p.begin(),p.end(),0);uint64_t best=UINT64_MAX;do{uint64_t c=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)c=(c<<1)|((a[p[i]]>>p[j])&1);best=min(best,c);}while(next_permutation(p.begin(),p.end()));return best;}
int main(int argc,char**argv){base=atoi(argv[1]);n=atoi(argv[2]);for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){idx[i][j]=idx[j][i]=pairs.size();pairs.push_back({i,j});}m=pairs.size();int size=1;for(int i=0;i<m;i++){powers.push_back(size);size*=base;}counts.resize(size);required.resize(m);vector<int> p(n);iota(p.begin(),p.end(),0);do{fill(required.begin(),required.end(),-1);for(int i=1;i<n;i++)required[idx[p[i-1]][p[i]]]=(p[i-1]>p[i]);scatter(0,0);}while(next_permutation(p.begin(),p.end()));map<int,long long> hist;map<pair<int,uint64_t>,pair<int,int>> classes;int floor=100000,strongct=0;long long sum=0;for(int code=0;code<size;code++){int h=counts[code];hist[h]++;sum+=h;if(base==3 && h!=3 && h!=7 && h!=21)continue;array<unsigned,8>a{};int t=code,recip=0;for(auto [i,j]:pairs){int d=t%base;t/=base;if(d!=1)a[i]|=1u<<j;if(d!=0)a[j]|=1u<<i;recip+=d==2;}bool sc=strong(a);if(base==2){if(sc){floor=min(floor,h);strongct++;}}else {auto key=make_pair(h,canon(a));classes[key]={sc,recip};}}
cout<<"base="<<base<<" n="<<n<<" graphs="<<size<<" total_path_incidences="<<sum<<" H7="<<hist[7]<<" H21="<<hist[21]<<"\n";if(base==2)cout<<"strong_count="<<strongct<<" strong_floor="<<floor<<"\n";map<tuple<int,int,int>,int> summary;for(auto [key,v]:classes){summary[{key.first,v.first,v.second}]++;cout<<"class H="<<key.first<<" bits="<<key.second<<" strong="<<v.first<<" reciprocal="<<v.second<<"\n";}for(auto [key,v]:summary){auto[h,s,r]=key;cout<<"summary H="<<h<<" strong="<<s<<" reciprocal="<<r<<" classes="<<v<<"\n";}ofstream f(string("hist-base")+to_string(base)+"-n"+to_string(n)+".csv");f<<"H,count\n";for(auto[h,c]:hist)if(c)f<<h<<","<<c<<"\n";
}
