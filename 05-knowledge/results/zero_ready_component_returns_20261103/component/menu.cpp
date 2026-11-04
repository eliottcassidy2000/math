#define main inherited_main
#include "../../triangle_return_170807_20261103/search.cpp"
#undef main
#include <algorithm>
#include <tuple>

int main(){
    int n;
    if(!(cin>>n)||n<3||n>18) return 2;
    T base(n);for(auto&x:base)cin>>x;
    vector<tuple<I,int,int,int,int,I,I>> nearest;
    int tested=0,strong_count=0,hits=0;
    cout<<"{\"status\":\"PASS\",\"hits\":[";
    for(int a=-1;a<n;a++)for(int b=a+1;b<n;b++){
        if(a==-1&&b!=0)continue;
        T t=base;if(a>=0){t[a]^=1<<b;t[b]^=1<<a;}
        ++tested;
        if(!strong(t))continue;
        ++strong_count;auto z=stats(t);
        for(int r=0;r<n;r++)for(int s=1;s<63;s++){
            if(z.E[r]>(512425-z.h)/(I(1)<<s)) break;
            I h=z.h+(I(1)<<s)*z.E[r];
            for(I target:{I(512423),I(512425)}){
                I delta=h-target;
                nearest.push_back({delta<0?-delta:delta,a,b,r,s,z.h,z.E[r]});
                if(h==target){
                    if(hits++)cout<<',';
                    cout<<"{\"arc\":["<<a<<','<<b<<"],\"root\":"<<r<<",\"prefix_size\":"<<s<<",\"target\":"<<target<<",\"h_B\":"<<z.h<<",\"F\":"<<z.E[r]<<",\"rows\":";
                    print(t);cout<<'}';
                }
            }
        }
    }
    sort(nearest.begin(),nearest.end());
    cout<<"],\"component_graphs\":"<<tested<<",\"strong_component_graphs\":"<<strong_count<<",\"nearest_below\":[";
    for(int i=0;i<min<int>(20,nearest.size());i++){
        if(i)cout<<',';auto [d,a,b,r,s,h,f]=nearest[i];
        cout<<'['<<d<<','<<a<<','<<b<<','<<r<<','<<s<<','<<h<<','<<f<<']';
    }
    cout<<"]}\n";
}
