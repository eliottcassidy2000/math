// Build: c++ -std=c++17 -O3 check_literal_boundary.cpp -o /tmp/literal_boundary
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <set>
#include <iterator>
#include <stdexcept>
#include <vector>
using namespace std;
using Rows=vector<int>;
void demand(bool b){if(!b)throw runtime_error("literal certificate mismatch");}
bool strong(const Rows&a){
    int n=a.size(), F=(1<<n)-1;
    for(int rev=0;rev<2;rev++){
        int S=1,todo=1;
        while(todo){
            int i=__builtin_ctz(unsigned(todo));todo&=todo-1;
            int b=rev?(F^a[i]^(1<<i)):a[i];
            todo|=b&~S; S|=b;
        }
        if(S!=F)return false;
    }
    return true;
}
struct Counts{
    int H[1024]{}, root[1024][10];
    Counts(const Rows&a,bool flip=false){
        int n=a.size(),F=(1<<n)-1, D[1024][10];
        H[0]=1;
        for(int S=1;S<=F;S++)for(int b=S;b;b&=b-1){
            int v=__builtin_ctz(unsigned(b)), T=S^(1<<v);
            int x=!T;
            for(int w=T&(flip?a[v]:~a[v]);w;w&=w-1)
                x+=D[T][__builtin_ctz(unsigned(w))];
            D[S][v]=x;H[S]+=x;
            root[S][v]=x;
        }
    }
};
int main(){
    const vector<Rows> input{
      {30,60,8,16,36,13},
      {158,252,8,16,36,13,189,60},
      {222,252,8,16,36,13,188,60}
    };
    const array<int,3> old{17,51,53};
    for(int c=0;c<3;c++){
        auto C=input[c];int m=C.size(), M=(1<<m)-1,q=old[c];
        demand(strong(C)&&Counts(C).H[M]==q);
        map<int,int> hits, full, exact;
        map<int,set<array<int,3>>> menus;
        for(int k=1;k<=2;k++){
            int n=m+k,F=(1<<n)-1;
            for(unsigned mask=0;mask<(1u<<(k*m));++mask){
                Rows a=C;
                for(int j=0;j<k;j++){
                    int U=(mask>>(m*j))&M;
                    a.push_back(U|(F^((1<<(m+j+1))-1)));  // old contacts + forward exterior
                    for(int b=M^U;b;b&=b-1)a[__builtin_ctz(unsigned(b))]|=1<<(m+j);
                }
                if(!strong(a))continue;
                Counts z(a);
                int h=z.H[F];
                if(h!=3*q&&h!=3*q+2&&h!=3*q+4&&!(c==1&&h==159))continue;
                ++hits[h];
                if(c>0){
                    Counts t(a,true);
                    for(int v=0;v<n;v++){
                        int E=0,O=0;
                        for(int S=1;S<=F;S++)if(S>>v&1)
                            (__builtin_popcount(unsigned(S))&1?O:E)+=t.root[S][v]*z.H[F^S];
                        int G=z.H[F^(1<<v)];
                        menus[h].insert({E,O,G});
                        if(E==81&&O==108&&G==27){
                            ++exact[h];
                            if(v==4&&h==3*q)++full[k];
                        }
                    }
                }
            }
        }
        cout<<"source "<<q<<":";
        for(int d:{0,2,4})cout<<" "<<3*q+d<<"="<<hits[3*q+d];
        if(c>0)cout<<" profile-matched zero at k=1,2: "<<full[1]<<","<<full[2];
        cout<<"\n";
        if(c==0)demand(hits[51]==10&&hits[53]==14&&hits[55]==0);
        if(c==1)demand(hits[153]==11&&hits[155]==2&&hits[157]==1&&full[1]==0&&full[2]==5);
        if(c==2)demand(hits[159]==10&&hits[161]==4&&hits[163]==0&&full[1]==2&&full[2]==4);
        if(c>0){
            for(int d:{0,2,4,6})cout<<"exact continuation "<<(3*q+d)<<"="<<exact[3*q+d]<<" ";
            cout<<"\n";
            demand(exact[3*q]==6&&exact[3*q+2]==0&&exact[3*q+4]==0);
            if(c==1)demand(exact[159]==3);
            for(int i:{0,2,4})for(int j:{0,2,4})if(i<j){
                vector<array<int,3>> both;
                set_intersection(menus[3*q+i].begin(),menus[3*q+i].end(),
                                 menus[3*q+j].begin(),menus[3*q+j].end(),back_inserter(both));
                cout<<"common arbitrary root profiles "<<3*q+i<<"/"<<3*q+j<<": "<<both.size()<<"\n";
                demand(both.empty());
            }
        }
    }
}
