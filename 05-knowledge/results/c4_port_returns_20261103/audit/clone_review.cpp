#include <algorithm>
#include <bit>
#include <cstdint>
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
using T=vector<unsigned>;
using I=int64_t;

void need(bool ok,const char*why){if(!ok)throw runtime_error(why);}
void valid(const T&t){int n=t.size();for(int i=0;i<n;i++){
    need(!(t[i]>>n)&&!(t[i]>>i&1),"row range or loop");
    for(int j=i+1;j<n;j++)need(((t[i]>>j&1)+(t[j]>>i&1))==1,"pair orientation");
}}
// Deletion recurrence, whereas the producer propagates paths forward.
I H(const T&t){int n=t.size(),N=1<<n;vector<I>d(N*n);vector<unsigned>in(n);
    for(int v=0;v<n;v++)for(int u=0;u<n;u++)if(t[u]>>v&1)in[v]|=1u<<u;
    for(int m=1;m<N;m++)for(unsigned b=m;b;b&=b-1){int v=countr_zero(b),rest=m^(1<<v);
        if(!rest){d[m*n+v]=1;continue;}
        for(unsigned q=rest&in[v];q;q&=q-1)d[m*n+v]+=d[rest*n+countr_zero(q)];
    }
    I answer=0;for(int v=0;v<n;v++)answer+=d[(N-1)*n+v];return answer;
}
// Directed searches from each vertex, not the producer's Floyd closure.
bool strong(const T&t){unsigned full=(1u<<t.size())-1;
    for(int root=0;root<(int)t.size();root++){unsigned seen=1u<<root,front=seen;
        while(front){int u=countr_zero(front);front&=front-1;unsigned fresh=t[u]&~seen;seen|=fresh;front|=fresh;}
        if(seen!=full)return false;
    }return true;
}
// Isomorphism by partial adjacency-preserving bijections: no bit-word keys.
bool iso(const T&a,const T&b){int n=a.size();if(b.size()!=a.size())return false;
    vector<int>da(n),db(n),order(n),map(n,-1);for(int i=0;i<n;i++){da[i]=popcount(a[i]);db[i]=popcount(b[i]);order[i]=i;}
    auto sa=da,sb=db;sort(sa.begin(),sa.end());sort(sb.begin(),sb.end());if(sa!=sb)return false;
    sort(order.begin(),order.end(),[&](int x,int y){int cx=count(da.begin(),da.end(),da[x]),cy=count(da.begin(),da.end(),da[y]);return cx!=cy?cx<cy:x<y;});
    function<bool(int,unsigned)>go=[&](int k,unsigned used){if(k==n)return true;int x=order[k];
        for(int y=0;y<n;y++)if(!(used>>y&1)&&da[x]==db[y]){bool ok=true;
            for(int j=0;j<k;j++){int z=order[j];if((a[x]>>z&1)!=(b[y]>>map[z]&1)){ok=false;break;}}
            if(ok){map[x]=y;if(go(k+1,used|(1u<<y)))return true;map[x]=-1;}
        }return false;
    };return go(0,0);
}
T clone(const T&t,int x,int count){int n=t.size();T out(n-1+count);vector<int>old;
    for(int i=0;i<n;i++)if(i!=x)old.push_back(i);
    for(int i=0;i<n-1;i++)for(int j=0;j<n-1;j++)if(t[old[i]]>>old[j]&1)out[i]|=1u<<j;
    for(int i=0;i<n-1;i++)for(int j=n-1;j<n-1+count;j++){
        if(t[old[i]]>>x&1)out[i]|=1u<<j;else out[j]|=1u<<i;
    }
    for(int i=n-1;i<n-1+count;i++)for(int j=i+1;j<n-1+count;j++)out[i]|=1u<<j;
    valid(out);return out;
}
vector<I>coefficients(const T&t,int x){int maxr=(t.size()-1)/2+1;vector<I>a;
    I stirling[10][10]{};stirling[0][0]=1;for(int i=1;i<=maxr;i++)for(int j=1;j<=i;j++)stirling[i][j]=stirling[i-1][j-1]+j*stirling[i-1][j];
    for(int k=1;k<=maxr;k++){I value=H(clone(t,x,k)),fact=1;
        for(int r=1;r<k;r++){fact*=r;value-=a[r-1]*fact*stirling[k][r];}
        fact*=k;need(value%fact==0,"nonintegral inferred run coefficient");a.push_back(value/fact);
    }return a;
}
int main(){int np;cin>>np;vector<T>reps;I births=0,sbirths=0,tbirths=0;
    for(int p=0;p<np;p++){T t(8);for(auto&r:t)cin>>r;valid(t);need(strong(t)&&H(t)<=97,"parent condition");
        for(unsigned cut=0;cut<256;cut++){T u=t;for(int i=0;i<8;i++)if(cut>>i&1)u[i]|=256;u.push_back(255^cut);births++;
            if(!strong(u))continue;sbirths++;I h=H(u);if(h!=91&&h!=97)continue;tbirths++;
            bool found=false;for(const auto&r:reps)if(iso(u,r)){found=true;break;}if(!found)reps.push_back(u);
        }
    }
    int ne;cin>>ne;set<int>matched;for(int k=0;k<ne;k++){T t(9);for(auto&r:t)cin>>r;valid(t);int found=-1;
        for(int i=0;i<(int)reps.size();i++)if(iso(t,reps[i])){need(found==-1,"duplicate independent classes");found=i;}
        need(found>=0&&matched.insert(found).second,"producer representative mismatch");
    }need(matched.size()==reps.size(),"producer omitted an independent class");
    cout<<"{\"parents\":"<<np<<",\"births\":"<<births<<",\"strong_births\":"<<sbirths<<",\"target_births\":"<<tbirths<<",\"class_counts\":[";
    for(int target:{91,97}){if(target==97)cout<<',';int c=0;for(auto&r:reps)c+=H(r)==target;cout<<c;}cout<<"],\"records\":[";
    int nr;cin>>nr;for(int k=0;k<nr;k++){int n,x;cin>>n>>x;T t(n);for(auto&r:t)cin>>r;valid(t);auto aa=coefficients(t,x);
        if(k)cout<<',';cout<<"{\"h\":"<<H(t)<<",\"strong\":"<<(strong(t)?"true":"false")<<",\"h2\":"<<H(clone(t,x,2))<<",\"A\":[";
        for(int i=0;i<(int)aa.size();i++){if(i)cout<<',';cout<<aa[i];}cout<<"]}";
    }cout<<"]}\n";
}
