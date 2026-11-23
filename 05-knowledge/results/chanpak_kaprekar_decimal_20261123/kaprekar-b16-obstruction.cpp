#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using namespace std;

static unsigned dp[1<<12][12];
int hp(const vector<int>& r){
    int n=(int)r.size(), F=(1<<n)-1;
    for(int s=1;s<=F;++s)for(unsigned m=s;m;m&=m-1){
        int i=__builtin_ctz(m), t=s^(1<<i); unsigned z=0;
        if(t)for(unsigned y=r[i]&t;y;y&=y-1)z+=dp[t][__builtin_ctz(y)];
        else z=1;
        dp[s][i]=z;
    }
    unsigned z=0;for(int i=0;i<n;++i)z+=dp[F][i];return (int)z;
}
bool strong(const vector<int>& r){
    int n=r.size(), F=(1<<n)-1;
    for(int i=0;i<n;++i){
        int s=1<<i;
        for(;;){
            int t=s;for(int j=0;j<n;++j)if((s>>j)&1)t|=r[j];
            if(s==t)break;
            s=t;
        }
        if(s!=F)return false;
    }
    return true;
}
int canonical(const vector<int>& r){
    vector<int> p(r.size());iota(p.begin(),p.end(),0);
    int best=1<<20;
    do{
        int z=0,k=0;
        for(int i=0;i<(int)p.size();++i)for(int j=i+1;j<(int)p.size();++j)
            z|=((r[p[i]]>>p[j])&1)<<k++;
        best=min(best,z);
    }while(next_permutation(p.begin(),p.end()));
    return best;
}
vector<vector<int>> cores(int n,int aim){
    vector<pair<int,int>> e;
    for(int i=0;i<n;++i)for(int j=i+2;j<n;++j)e.emplace_back(i,j);
    map<int,vector<int>> out;
    for(int m=0;m<(1<<(int)e.size());++m){
        vector<int> r(n);
        for(int i=0;i+1<n;++i)r[i]|=1<<(i+1);
        for(int j=0;j<(int)e.size();++j){
            auto [a,b]=e[j];if((m>>j)&1)r[a]|=1<<b;else r[b]|=1<<a;
        }
        if(hp(r)==aim&&strong(r))out.emplace(canonical(r),r);
    }
    vector<vector<int>> z;
    for(auto& [_,r]:out)z.push_back(r);
    return z;
}
vector<int> ordinal(const vector<vector<int>>& all){
    int n=0;for(auto& r:all)n+=(int)r.size();
    vector<int> z;
    for(auto& r:all){
        int a=z.size(), end=a+(int)r.size();
        for(int x:r)z.push_back((x<<a)|(((1<<n)-1)^((1<<end)-1)));
    }
    return z;
}
int main(){
    auto c3=cores(3,3), c5=cores(4,5), c15a=cores(5,15), c15b=cores(6,15);
    assert(c3.size()==1 && c5.size()==1 && c15a.size()==2 && c15b.size()==2);
    vector<vector<int>> c15=c15a;c15.insert(c15.end(),c15b.begin(),c15b.end());
    vector<vector<int>> test;
    constexpr int n=12;
    for(auto& r:c15)for(int a=0;a<=n-(int)r.size();++a){
        vector<vector<int>> B(a,{0});B.push_back(r);B.resize(n-r.size()+1,{0});
        test.push_back(ordinal(B));
    }
    for(int sign=0;sign<2;++sign)for(int a=0;a<=n-7;++a)for(int b=0;b<=n-7-a;++b){
        vector<vector<int>> B(a,{0});B.push_back(sign?c5[0]:c3[0]);
        B.insert(B.end(),b,{0});B.push_back(sign?c3[0]:c5[0]);
        B.resize(n-7+2,{0});test.push_back(ordinal(B));
    }
    assert(test.size()==72);
    long hits=0;
    for(auto A:test){
        assert(A.size()==n&&hp(A)==15);
        for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
            A[i]^=1<<j; A[j]^=1<<i;
            hits+=hp(A)==225;
            A[i]^=1<<j; A[j]^=1<<i;
        }
    }
    assert(hits==0);
    cout<<"72 SCC profiles; 4752 edge reversals; zero 15-to-225 seams\n";
}
