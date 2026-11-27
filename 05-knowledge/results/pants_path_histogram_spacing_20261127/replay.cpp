#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

using u32 = std::uint32_t;
using i64 = std::int64_t;
using Matrix = std::vector<std::vector<int>>;
using Hist = std::map<int,int>;

Matrix make(const std::vector<std::string>& rows,int start,int n) {
    Matrix A(n,std::vector<int>(n));
    for(int i=0;i<n-1;++i) {
        A[n-1][i]=1; A[i][n-1]=-1;
        for(int j=0;j<n-1;++j) if(i!=j) A[i][j]=rows[start+i][j]=='1'?1:-1;
    }
    for(int i=0;i<n;++i) for(int j=0;j<n;++j) {
        assert(A[i][j]==-A[j][i]);
        assert(std::inner_product(A[i].begin(),A[i].end(),A[j].begin(),0)==(i==j ? n-1 : 0));
    }
    return A;
}

int quartic(const Matrix& A,int a,int b,int c,int d) {
    int sum=0;
    for(int j=0;j<(int)A.size();++j) sum+=A[a][j]*A[b][j]*A[c][j]*A[d][j];
    return sum;
}

void both_triples(const Matrix& B,int r,int a,int b,int c,int& D0,int& D1) {
    int K=B[a][c]*B[b][c]+B[a][b]*B[c][b]+B[b][a]*B[c][a];
    int Q=quartic(B,r,a,b,c);
    D0=(int)B.size()-1+K+Q; D1=(int)B.size()-1+K-Q;
    assert(D0%8==0 && D1%8==0);
}

Hist verify(const Matrix& A,bool shell4) {
    int n=A.size();
    i64 L=n*(n-2);
    i64 s0=L*L, s1=s0-8*L*(n-4), s2=2*s1-s0+48*n*(n-4);
    i64 base3=s0-3*s1+3*s2, base4=3*s0-8*s1+6*s2;
    std::vector<u32> rows(n);
    std::vector<int> powers(n);
    for(int j=0;j<n;++j) {
        int m=n-1-2*j;
        powers[j]=(m*m-1)*(m*m-1);
    }
    for(int i=0;i<n;++i) for(int j=0;j<n;++j) if(A[i][j]<0) rows[i]|=1u<<j;
    auto score=[&](u32 cut){
        int s=0;
        for(int i=0;i<n;++i) s+=powers[std::popcount((rows[i]^cut)&~(1u<<i))];
        return s;
    };
    Hist H;
    for(int r=0;r<n;++r) {
        Matrix B=A;
        for(int i=0;i<n;++i) for(int j=0;j<n;++j)
            if(((rows[r]>>i)^(rows[r]>>j))&1u) B[i][j]=-B[i][j];
        std::vector<int> v;
        for(int i=0;i<n;++i) if(i!=r) v.push_back(i);
        std::vector<std::array<int,2>> d(n*n*n);
        auto index=[&](int a,int b,int c){return (a*n+b)*n+c;};
        int g=0, first=-1;
        for(int a:v) for(int b:v) if(a<b) for(int c:v) if(b<c) {
            auto& e=d[index(a,b,c)];
            both_triples(B,r,a,b,c,e[0],e[1]);
            for(int u=0;u<2;++u) {
                u32 cut=rows[r]^(1u<<a)^(1u<<b)^(1u<<c)^(u<<r);
                const int exact=score(cut);
                assert(exact==base3-192*e[u]);
                if(first<0) first=exact;
                g=std::gcd(g,exact-first);
            }
        }
        // The signed four-cut certificate has score difference 1536.
        assert(g==1536);
        if(!shell4) continue;
        for(int a:v) for(int b:v) if(a<b) for(int c:v) if(b<c) for(int e:v) if(c<e) {
            i64 q4=quartic(B,a,b,c,e);
            auto u=d[index(a,b,c)],v1=d[index(a,b,e)],w=d[index(a,c,e)],x=d[index(b,c,e)];
            u32 flips=rows[r]^(1u<<a)^(1u<<b)^(1u<<c)^(1u<<e);
            for(int k=0;k<2;++k) {
                int s=score(flips^(k<<r));
                assert(s==base4+384*q4-192*(u[k]+v1[k]+w[k]+x[k]));
                assert(s>=136192 && s<=159040);
                ++H[s];
            }
        }
    }
    return H;
}

int main(int argc,char** argv) {
    assert(argc==2);
    std::ifstream in(argv[1]);
    assert(in);
    std::map<int,std::vector<std::string>> rows;
    for(std::string s;std::getline(in,s);) {
        const auto a=s.find_first_not_of(" \t\r");
        if(a==std::string::npos) continue;
        s=s.substr(a,s.find_last_not_of(" \t\r")-a+1);
        if((s.size()==23 || s.size()==27) && s.find_first_not_of("01")==std::string::npos)
            rows[s.size()].push_back(s);
    }
    assert(rows[23].size()==368 && rows[27].size()==1755);
    for(int n:{8,12,20,32}) {
        int p=n-1;
        std::vector<std::string> paley(p,std::string(p,'0'));
        for(int i=0;i<p;++i) for(int j=1;j<p;++j) paley[i][(i+j*j)%p]='1';
        verify(make(paley,0,n),false);
    }
    for(int i=0;i<16;++i) verify(make(rows[23],i*23,24),false);
    std::vector<Hist> hist;
    std::set<Hist> signatures;
    for(int i=0;i<65;++i) {
        auto h=verify(make(rows[27],i*27,28),true);
        int total=0;
        for(const auto& [s,f]:h) total+=f;
        assert(total==982800);
        signatures.insert(h);
        hist.push_back(h);
    }
    assert(signatures.size()==41);
    std::cout<<"frequency at 8474100:";
    for(auto& h:hist) std::cout<<" "<<h.at(139840);
    std::cout<<"\ncomplete single path-count bins:";
    i64 falling=1;
    for(int j=23;j<=28;++j) falling*=j;
    std::set<int> all_scores;
    for(const auto& h:hist) for(const auto& [s,k]:h) all_scores.insert(s);
    for(int s:all_scores) {
        std::set<int> values;
        for(const auto& h:hist) {
            auto it=h.find(s);
            values.insert(it==h.end()?0:it->second);
        }
        if(values.size()==41) std::cout<<" "<<(falling+28*(3*27*27-2*27)-s-28*53)/32;
    }
    std::cout<<"\nchecked local gcd on every root and both formulas at 63882000 fourth-shell locations\n";
}
