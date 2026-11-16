// Dependency-free Node.js certificate. Run: node check_ternary.js
function run(log){
  const C=[118,124,248,113,224,192,128,11], A=[
    [374,380,504,369,224,448,128,11,208],
    [118,372,504,115,224,192,128,11,249],
    [118,124,504,369,224,192,128,11,243]
  ], M0=6, root=7;
  const check=(f,msg)=>{if(!f)throw Error(msg);};
  const bit=x=>31-Math.clz32(x), pop=x=>{let p=0;while(x){p++;x&=x-1}return p};
  function validate(a){
    for(let i=0;i<a.length;i++){
      check(!(a[i]>>i&1)&&a[i]<(1<<a.length),"bad row");
      for(let j=0;j<i;j++)check((a[i]>>j&1)+(a[j]>>i&1)===1,"not a tournament");
    }
    const F=(1<<a.length)-1;
    for(let c=0;c<2;c++){let seen=1,b=1;
      while(b){const i=bit(b&-b);b&=b-1;const x=(c?(F^a[i]^(1<<i)):a[i])&~seen;seen|=x;b|=x;}
      check(seen===F,"non-strong");
    }
  }
  function dp(a,root=-1){
    const n=a.length,F=(1<<n)-1,D=new Int32Array((F+1)*n),H=new Int32Array(F+1);H[0]=1;
    for(let i=0;i<n;i++)if(root<0||root===i)D[(1<<i)*n+i]=1;
    for(let S=1;S<=F;S++)for(let b=S;b;b&=b-1){
      const i=bit(b&-b),x=D[S*n+i]; if(!x)continue; H[S]+=x;
      for(let c=a[i]&~S;c;c&=c-1){const v=c&-c;D[(S|v)*n+bit(v)]+=x;}
    }
    return H;
  }
  function state(a,r){
    const H=dp(a),P=dp(a,r),F=H.length-1;
    let E=0,O=0;
    for(let S=1;S<=F;S++){const z=P[S]*H[F^S];if(pop(S)&1)O+=z;else E+=z;}
    return [H[F],E,O,H[F^(1<<r)]];
  }
  function skeletons(a,removed=0,first=-1){
    const F=((1<<a.length)-1)^((1<<M0)|removed);
    const d=new Int32Array((1<<a.length)*3);
    if(first<0)d[0]=1;
    function visit(S,k,v){
      d[3*S+k]++;
      for(let b=a[v];b;b&=b-1){
        const j=bit(b&-b);
        if(j===M0){if(k<2)visit(S,k+1,j);}
        else if((F>>j&1)&&!(S>>j&1))visit(S|(1<<j),k,j);
      }
    }
    for(let v=0;v<a.length;v++)if(first<0||v===first){
      if(v===M0)visit(0,1,v);else if(F>>v&1)visit(1<<v,0,v);
    }
    return [F,d];
  }
  function table(a){
    const [F,H]=skeletons(a),[K,G]=skeletons(a,1<<root),[,P]=skeletons(a,0,root);
    let E=[0,0,0],O=[0,0,0];
    for(let S=0;S<=F;S++)if((F&S)===S)for(let i=0;i<3;i++)if(P[3*S+i])for(let j=0;j<3;j++)if(H[3*(F^S)+j]&&i+j){
      check(i+j<=2&&i<=1&&j<=1,"unanticipated rooted pattern");
      const z=P[3*S+i]*H[3*(F^S)+j], idx=i+j===2?2:i;
      E[idx]+=z;O[idx]+=(pop(S)&1?-1:1)*z;
    }
    return [[H[3*F+1],H[3*F+2]],[G[3*K+1],G[3*K+2]],E,O];
  }
  function substitute(a,M){
    let w=[];
    for(let i=0;i<a.length;i++)for(let j=0;j<(i===M0?M.length:1);j++)w.push([i,j]);
    return [w.map(([i,k])=>w.reduce((r,[j,l],t)=>r+(((a[i]>>j&1)||(i===M0&&j===M0&&(M[k]>>l&1)))?1<<t:0),0)),w.findIndex(([i])=>i===root)];
  }
  const TT=n=>Array.from({length:n},(_,i)=>((1<<n)-1)^((1<<(i+1))-1));
  const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
  const flipped=C.map((v,i)=>v^(i===1?8:i===3?2:0));
  check(eq(table(flipped),[[169,88],[5,0],[93,88,88],[-5,0,0]]),"preparation coefficients");
  for(let s=0;s<=1;s++){
    const [a,r]=substitute(C,TT(s)),F=(1<<a.length)-1,q=state(a,r)[0],hits=[];
    for(let c=0;c<=F;c++){
      const b=a.map((v,i)=>v|((1-(c>>i&1))<<a.length)).concat(c);
      const p=state(b,r);if(p[0]===3*q+2)hits.push([c,p]);
    }
    check(eq(hits,s?[]:[[104,[239,127,142,15]]]),"unchanged-source single-append census");
  }
  const profiles=[];
  for(let i=0;i<4;i++){
    const a=[C,...A][i],f=i===0?1:3,d=i===0?0:2*(i-1);
    validate(a);check(a[M0]===(1<<root),"marked module has another exit");
    if(i>0)check(eq(a.slice(0,8).map((m,v)=>(m&255)^C[v]),i===2?[0,8,0,2,0,0,0,0]:[0,0,0,0,0,0,0,0]),"induced/rewired status");
    const tab=table(a);
    check(eq(tab,[[163*f+d,84*f],[5*f,0],[89,84,84].map(x=>f*x),[-5*f,0,0]]),"run coefficients");
    const cases=[];
    for(const M of [TT(0),TT(1),TT(2),TT(3),TT(4),[2,4,1],[6,12,8,1]]){
      const DH=dp(M),F=DH.length-1,u=DH[F],R=DH.reduce((sum,h,S)=>sum+h*DH[F^S],0);
      const expected=[f*(84*R-5*u)+d*u,42*f*R,f*(42*R+5*u),5*f*u];
      const [b,r]=substitute(a,M); validate(b);
      const got=state(b,r); check(eq(expected,got),"independent direct subset DP");
      cases.push(got);
    }
    profiles.push({i,tab,cases});
  }
  log(profiles);
  return profiles;
}
run(x=>console.log(JSON.stringify(x,null,2)));
