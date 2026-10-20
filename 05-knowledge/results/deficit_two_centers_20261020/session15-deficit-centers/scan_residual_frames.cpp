#include <bits/stdc++.h>
using namespace std;
using U=uint32_t;
struct Desc {int fixed; vector<U> b[2];};
static bool strong(const vector<U>&a,U rem){
    if(!rem)return false;
    U start=rem&-rem;
    U seen=start,front=start;
    while(front){int v=__builtin_ctz(front);front&=front-1;U add=a[v]&rem&~seen;seen|=add;front|=add;}
    if(seen!=rem)return false;
    seen=start;front=start;
    while(front){int v=__builtin_ctz(front);front&=front-1;U x=rem&~seen,add=0;while(x){int w=__builtin_ctz(x);x&=x-1;if(a[w]>>v&1)add|=U(1)<<w;}seen|=add;front|=add;}
    return seen==rem;
}
int main(int argc,char**argv){
    if(argc!=2)return 2;ifstream in(argv[1]);int n,q,K,F;in>>n>>q>>K>>F;
    vector<Desc>d(K);vector<int>freeidx;vector<U>a(n,0);
    for(int k=0;k<K;k++){in>>d[k].fixed;for(int b=0;b<2;b++){d[k].b[b].resize(n);for(auto&x:d[k].b[b])in>>x;}
        if(d[k].fixed<0)freeidx.push_back(k);int b=max(0,d[k].fixed);for(int i=0;i<n;i++)a[i]|=d[k].b[b][i];}
    assert((int)freeidx.size()==F);
    vector<array<int,4>> subs;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int k=j+1;k<n;k++)for(int l=k+1;l<n;l++)subs.push_back({i,j,k,l});
    U full=(U(1)<<n)-1;uint64_t N=uint64_t(1)<<F,passed=0;vector<uint64_t>survive;map<int,uint64_t>hist;
    ofstream cert(string(argv[1])+".survivors");
    ofstream witnesses(string(argv[1])+".witnesses.bin",ios::binary);
    auto t0=chrono::steady_clock::now();
    for(uint64_t it=0;it<N;it++){
        if(it){int changed=__builtin_ctzll(it);int k=freeidx[changed];for(int i=0;i<n;i++)a[i]^=d[k].b[0][i]^d[k].b[1][i];}
        bool ok=false;int pos=0;
        for(auto S:subs){pos++;U mask=0;for(int v:S)mask|=U(1)<<v;
            bool s4=true;for(int v:S){int deg=__builtin_popcount(a[v]&mask);if(deg==0||deg==3){s4=false;break;}}if(!s4)continue;
            U x=a[S[0]],y=a[S[1]],z=a[S[2]],w=a[S[3]];
            U terminal=((x&y&z)|(x&y&w)|(x&z&w)|(y&z&w))&~mask&full;if(!terminal)continue;
            if(!strong(a,full^mask))continue;
            ok=true;hist[pos]++; for(int v:S)witnesses.put(char(v));witnesses.put(char(__builtin_ctz(terminal)));break;
        }
        if(ok)passed++;else {uint64_t gray=it^(it>>1);survive.push_back(gray);for(int j=0;j<5;j++)witnesses.put(char(255));cert<<gray;for(U x:a)cert<<' '<<x;cert<<'\n';}
        if((it+1)%262144==0)cerr<<argv[1]<<" "<<it+1<<"/"<<N<<" survivors "<<survive.size()<<"\n";
    }
    double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();
    cout<<"{\"file\":\""<<argv[1]<<"\",\"n\":"<<n<<",\"q\":"<<q<<",\"presentations\":"<<N<<",\"prefix_excluded\":"<<passed<<",\"survivors\":"<<survive.size()<<",\"seconds\":"<<sec<<",\"distinct_subset_positions\":"<<hist.size()<<",\"last_subset_position\":"<<(hist.empty()?0:hist.rbegin()->first)<<",\"surviving_masks\":[";
    for(int i=0;i<(int)survive.size();i++){if(i)cout<<',';cout<<survive[i];}cout<<"]}"<<endl;
}
