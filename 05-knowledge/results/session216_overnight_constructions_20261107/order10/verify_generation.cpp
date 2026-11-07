#define main generator_main
#include "generate.cpp"
#undef main
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
vector<G>read(const string&path){ifstream in(path);int total;in>>total;vector<G>out;for(int k=0;k<total;k++){int n;in>>n;G a(n);for(int&x:a)in>>x;for(int i=0;i<n;i++){need((a[i]>>i&1)==0,"loop");for(int j=0;j<i;j++)need(((a[i]>>j&1)+(a[j]>>i&1))==1,"edge");}out.push_back(a);}return out;}
U brute(const G&a){vector<int>p(a.size());iota(p.begin(),p.end(),0);U best=~U(0);do{best=min(best,encode(a,p));}while(next_permutation(p.begin(),p.end()));return best;}
int main(){long long small=0,idempotent=0,relabels=0;for(int n=1;n<=5;n++){map<U,U>forward,backward;for(U bits=0;bits<(U(1)<<(n*(n-1)/2));bits++){auto a=decode(bits,n);U x=canonical(a),y=brute(a);need(!forward.count(x)||forward[x]==y,"canon conflation");need(!backward.count(y)||backward[y]==x,"canon split");forward[x]=y;backward[y]=x;small++;}}
auto old=read("dependencies/BANK_THROUGH8.txt"),fresh=read("BANK8.txt"),nine=read("BANK9.txt");set<U>A,B;for(auto&a:old)if(a.size()==8)A.insert(canonical(a));for(auto&a:fresh)B.insert(canonical(a));need(A==B&&A.size()==6880,"inherited order8 disagreement");
mt19937 rng(216);for(auto&a:nine){vector<int>p(a.size());iota(p.begin(),p.end(),0);U c=encode(a,p);need(canonical(a)==c,"idempotence");idempotent++;}
for(int k=0;k<2000;k++){auto&a=nine[(k*971)%nine.size()];vector<int>p(a.size());iota(p.begin(),p.end(),0);U c=canonical(a);shuffle(p.begin(),p.end(),rng);G g(a.size());for(int i=0;i<(int)a.size();i++)for(int j=0;j<(int)a.size();j++)if(a[p[i]]>>p[j]&1)g[i]|=1<<j;need(canonical(g)==c,"relabel invariance");relabels++;}
cout<<"{\"status\":\"PASS\",\"brute_force_graphs\":"<<small<<",\"inherited_order8_classes\":"<<A.size()<<",\"order9_idempotence\":"<<idempotent<<",\"order9_relabels\":"<<relabels<<"}\n";
}
