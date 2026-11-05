#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using I=long long;
using G=vector<unsigned>;
void need(bool x,const string&s){
    if(!x)throw runtime_error(s);
}
struct Tab{
    vector<I> h;
    vector<vector<I>> ends;
};
Tab tab(const G&a){
    int n=a.size(),N=1<<n;
    Tab t{
        vector<I>(N),vector<vector<I>>(N,vector<I>(n))
    };
    t.h[0]=1;
    for(int mask=1;mask<N;mask++){
        for(int v=0;v<n;v++)if(mask>>v&1){
            int old=mask^(1<<v);
            I x=old?0:1;
            for(int u=0;u<n;u++)if((old>>u&1)&&(a[u]>>v&1))x+=t.ends[old][u];
            t.ends[mask][v]=x;
            t.h[mask]+=x;
        }
    }
    return t;
}
G converse(G a){
    unsigned f=(1u<<a.size())-1;
    for(unsigned i=0;i<a.size();i++)a[i]=f^a[i]^(1u<<i);
    return a;
}
struct P{
    I h,e,d,R,L;
    I U()const{
        return R-d;
    }
    I s()const{
        return R-h-d+e;
    }
};
vector<P> profiles(const G&a){
    auto t=tab(a),q=tab(converse(a));
    int n=a.size(),f=(1<<n)-1;
    vector<P>p(n);
    for(int r=0;r<n;r++){
        p[r]={
            t.h[f],t.ends[f][r],t.h[f^(1<<r)],0,0
        };
        for(int m=1;m<=f;m++)if(m>>r&1){
            p[r].R+=q.ends[m][r]*t.h[f^m];
            p[r].L+=t.ends[m][r]*t.h[f^m];
        }
    }
    return p;
}
I ham(const G&a){
    return tab(a).h.back();
}
G tour(int n,unsigned long long code){
    G a(n);
    int b=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,b++){
        if(code>>b&1)a[i]|=1u<<j;
        else a[j]|=1u<<i;
    }
    return a;
}
void edge(G&a,int i,int j){
    a[i]|=1u<<j;
}
G clone(const G&a,int r,int m){
    int n=a.size();
    vector<int>own;
    for(int i=0;i<n;i++)for(int k=0;k<(i==r?m:1);k++)own.push_back(i);
    G b(own.size());
    for(int i=0;i<(int)b.size();i++)for(int j=i+1;j<(int)b.size();j++){
        if(own[i]==own[j]||(a[own[i]]>>own[j]&1))edge(b,i,j);
        else edge(b,j,i);
    }
    return b;
}
G apex(unsigned q,int n){
    G a(n+1);
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++)edge(a,i,j);
        if(q>>(n-1-i)&1)edge(a,n,i);
        else edge(a,i,n);
    }
    return a;
}
G compiler(const G&y,int c,int m){
    int n=y.size(),a=n,b=n+1,z=n+2;
    G x=y;
    x.resize(n+2+m);
    edge(x,c,a);
    edge(x,c,b);
    edge(x,a,b);
    for(int v=0;v<n;v++)if(v!=c){
        edge(x,a,v);
        edge(x,b,v);
    }
    for(int i=z;i<(int)x.size();i++){
        edge(x,a,i);
        edge(x,b,i);
        edge(x,i,c);
        for(int v=0;v<n;v++)if(v!=c)edge(x,i,v);
        for(int j=i+1;j<(int)x.size();j++)edge(x,i,j);
    }
    return x;
}
bool strong(const G&a){
    for(int r=0;r<(int)a.size();r++){
        unsigned seen=1u<<r,old=0;
        while(old!=seen){
            old=seen;
            for(int v=0;v<(int)a.size();v++)if(seen>>v&1)seen|=a[v];
        }
        if(seen!=((1u<<a.size())-1))return false;
    }
    return true;
}
bool prime(I x){
    if(x<2)return false;
    for(I i=2;i*i<=x;i++)if(x%i==0)return false;
    return true;
}
I sigma(unsigned q){
    if(!q)return 0;
    if(!(q&1))return sigma(q>>1);
    return (q>>1)+sigma(q>>1)+((q>>1)&1);
}
I fh(unsigned q){
    if(!q)return 1;
    if(q&1)return fh(q>>1);
    return fh(q>>1)+(q>>1)+((q>>1)&1);
}
int match_rank(const G&a,int r){
    vector<int>o,in;
    for(int v=0;v<(int)a.size();v++)if(v!=r){
        if(a[r]>>v&1)o.push_back(v);
        else in.push_back(v);
    }
    int best=0;
    auto rec=[&](auto&&self,int i,unsigned used,int k)->void{
        best=max(best,k);
        if(i==(int)o.size())return;
        self(self,i+1,used,k);
        for(int j=0;j<(int)in.size();j++)if(!(used>>j&1)&&(a[o[i]]>>in[j]&1))self(self,i+1,used|(1u<<j),k+1);
    };
    rec(rec,0,0,0);
    return best;
}
// Literal path and two-cover counts do not use the subset recurrence.
P literal(const G&a,int r){
    int n=a.size();
    vector<int>v;
    for(int i=0;i<n;i++)v.push_back(i);
    P p{
        0,0,0,0,0
    };
    do{
        bool h=true,rr=true,ll=true;
        for(int i=0;i+1<n;i++){
            bool ok=a[v[i]]>>v[i+1]&1;
            h&=ok;
            if(v[i+1]!=r)rr&=ok;
            if(v[i]!=r)ll&=ok;
        }
        p.h+=h;
        p.e+=h&&v.back()==r;
        p.R+=rr;
        p.L+=ll;
    }
    while(next_permutation(v.begin(),v.end()));
    v.clear();
    for(int i=0;i<n;i++)if(i!=r)v.push_back(i);
    do{
        bool h=true;
        for(int i=0;i+1<(int)v.size();i++)h&=a[v[i]]>>v[i+1]&1;
        p.d+=h;
    }
    while(next_permutation(v.begin(),v.end()));
    return p;
}
int main(){
    try{
        I graphs=0,roots=0,up1=0,lone=0,words=0,neutral=0,compilers=0,compiledCounts=0,strongchecks=0,primepairs=0,evenprimepairs=0;
        set<array<I,5>>xs,ys;
        vector<pair<G,int>>unit;
        for(int n=1;n<=6;n++){
            for(unsigned long long code=0;code<(1ULL<<(n*(n-1)/2));code++){
                auto a=tour(n,code);
                auto pp=profiles(a);
                graphs++;
                for(int r=0;r<n;r++){
                    auto p=pp[r];
                    roots++;
                    if(n<=4){
                        auto l=literal(a,r);
                        need(p.h==l.h&&p.e==l.e&&p.d==l.d&&p.R==l.R&&p.L==l.L,"literal profile crosscheck");
                    }
                    need(p.e<=p.d&&p.d<=p.h&&p.h<=p.R,"physical inequality");
                    need(p.s()>=0&&p.s()%2==0,"slack parity");
                    if(p.U()==1){
                        up1++;
                        need(p.h==1&&p.e==0&&p.R==2&&p.d==1&&p.s()==0,"U1 lemma");
                    }
                    if(p.L==1){
                        lone++;
                        need(p.h==1&&p.d==1&&p.R>0,"L1 transitive");
                        need(a[r]==(((1u<<n)-1)^(1u<<r)),"L1 source");
                    }
                    xs.insert({
                        p.h,p.e,p.R,p.d,p.s()
                    }
                    );
                    if(n<=5)ys.insert({
                        p.h,p.d,p.L,0,0
                    }
                    );
                    if(n<=5&&match_rank(a,r)==1){
                        I b=2*p.h-ham(clone(a,r,2));
                        if(b==1||b==-1)unit.push_back({
                            a,r
                        }
                        );
                    }
                }
            }
            cerr<<"census n="<<n<<" done\n";
        }
        for(auto x:xs)for(auto y:ys){
            I U=x[2]-x[3],L=y[2],al=U*L,b=x[3]*y[0]-x[4]*y[1];
            if(prime(al)&&(b==1||b==-1)){
                if(al%2){
                    primepairs++;
                    need(U==al&&L==1&&y[0]==1&&y[1]==1&&x[4]==x[3]-b,"odd prime reduction");
                    if(x[3]==1)need(b==-1&&(al==3||al==5),"odd prime d1");
                }
                else{
                    evenprimepairs++;
                    need(b==1&&x[3]==1&&x[4]==0,"alpha2 classification");
                }
            }
        }
        for(int n=0;n<=10;n++)for(unsigned q=0;q<(1u<<n);q++){
            auto a=apex(q,n);
            auto p=profiles(a).back();
            words++;
            need(p.d==1&&p.h==fh(q)&&p.R==q+1&&p.U()==q&&p.e==1-(q%2)&&p.s()==sigma(q),"apex graph recurrence");
            if(n<=7){
                auto z=profiles(apex(q,n+2)).back();
                neutral++;
                need(p.h==z.h&&p.e==z.e&&p.d==z.d&&p.R==z.R,"leading neutral padding");
            }
        }
        for(unsigned q=0;q<(1u<<20);q++){
            I s=sigma(q);
            if(q){
                unsigned odd=q;
                while(!(odd&1))odd>>=1;
                need(odd<=2*s+1&&sigma(odd)==s,"fixed slack bound");
                bool pow2=!(q&(q-1));
                need((s==0)==pow2,"sigma0 iff");
                need((s==2)==(odd==3||odd==5),"sigma2 iff");
            }
        }
        // Exhaustive actual sign compilers on labelled inputs through n=4; one
        // representative of every (H,b,L0,L1) signature through n=5.
        set<array<I,4>>sig;
        for(auto [a,r]:unit){
            auto y=clone(a,r,2);
            auto yp=profiles(y);
            I h=ham(a),b=2*h-ham(y);
            array<I,4>key={
                h,b,yp[r].L,yp[r+1].L
            };
            bool fresh=sig.insert(key).second;
            if(a.size()>4&&!fresh)continue;
            for(int c:{
                r,r+1
            }
            ){
                auto py=yp[c];
                need(py.d==h&&py.h==2*h-b,"cloned Y deletion");
                bool st=strong(a);
                G f=compiler(y,c,1);
                int z=y.size()+2;
                need(f.size()==a.size()+4,"compiler vertex count");
                I N=6*py.L-b;
                need(N>=12*h-7*b,"compiler lower bound");
                need(ham(f)==N,"physical compiler m1");
                for(int m=1;m<=3;m++){
                    G out=compiler(y,c,m);
                    need(ham(out)==3*py.L*(1LL<<m)-b,"physical compiler formula");
                    if(m>1)need(out==clone(f,z,m),"literal clone graph equality");
                    if(st){
                        need(strong(out),"strong compiler");
                        strongchecks++;
                    }
                    compiledCounts++;
                }
                need(2*ham(f)-ham(clone(f,z,2))==-b,"new offset sign");
                need(match_rank(f,z)==1,"new star rank");
                need(3*py.L%3==0&&3*py.L>3,"composite amplitude");
                compilers++;
            }
        }
        // New exact odd-core exclusion and the larger-deletion positive contrast.
        vector<I> ten_certificate={0,2,2,6,4,8,8,14,8,14,12};
        for(int j=0;j<11;j++) need(sigma(2*j+1)==ten_certificate[j]&&ten_certificate[j]!=10,"slack10 certificate");
        P blocked=profiles(G{14,12,24,16,3})[4];
        need(blocked.h==11&&blocked.d==1&&blocked.L==13,"actual blocked partner");
        P contrast=profiles(G{16,1,3,7,14})[2];
        need(contrast.h==9&&contrast.e==2&&contrast.R==16&&contrast.d==5,"alpha11 input profile");
        G first_contrast={6,48,34,39,45,1};
        for(int m=1;m<=5;m++){
            G out=clone(first_contrast,5,m);
            need(ham(out)==11*(1LL<<m)+1&&strong(out),"alpha11 actual graph ray");
        }
        // U=0 is deliberately outside the finite positive-amplifier corollary.
        need(profiles(G{0})[0].U()==0,"singleton U0");
        for(int m=1;m<=5;m++)need(ham(apex(0,m))==1,"U0 neutral synthesis");
        // Two successive sign switches, with every choice of the marked Y clone.
        for(int first=0;first<2;first++)for(int second=0;second<2;second++){
            G a={
                2,4,1
            };
            int r=0;
            I b=1;
            for(int choice:{
                first,second
            }
            ){
                I h=ham(a);
                auto y=clone(a,r,2);
                int c=r+choice;
                I L=profiles(y)[c].L;
                auto f=compiler(y,c,1);
                r=y.size()+2;
                need(ham(f)==6*L-b&&ham(y)==2*h-b,"iterated compiler determinant");
                need(2*ham(f)-ham(clone(f,r,2))==-b,"iterated compiler sign");
                need(strong(f),"iterated compiler strongness");
                a=f;
                b=-b;
            }
        }
        cout<<"{\n  \"status\": \"PASS\",\n  \"labelled_graphs_through6\": "<<graphs<<",\n  \"rooted_graphs_through6\": "<<roots<<",\n  \"U1_roots\": "<<up1<<",\n  \"L1_roots\": "<<lone<<",\n  \"distinct_X_profiles\": "<<xs.size()<<",\n  \"distinct_Y_profiles_through5\": "<<ys.size()<<",\n  \"odd_prime_unit_profile_pairs\": "<<primepairs<<",\n  \"even_prime_unit_profile_pairs\": "<<evenprimepairs<<",\n  \"apex_words_through_length10\": "<<words<<",\n  \"leading_zero_padding_checks\": "<<neutral<<",\n  \"arithmetic_integers_checked\": 1048576,\n  \"slack10_exact_odd_core_checks\": 11,\n  \"alpha11_actual_graph_checks\": 5,\n  \"two_step_compiler_choice_sequences\": 4,\n  \"unit_input_roots_through5\": "<<unit.size()<<",\n  \"unit_signature_classes\": "<<sig.size()<<",\n  \"physical_compilers_both_clone_choices\": "<<compilers<<",\n  \"physical_compiler_count_checks\": "<<compiledCounts<<",\n  \"strongness_checks\": "<<strongchecks<<"\n}\n";
    }
    catch(const exception&e){
        cerr<<e.what()<<"\n";
        return 1;
    }
}
