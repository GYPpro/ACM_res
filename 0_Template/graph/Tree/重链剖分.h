/*
n个节点树，节点有权值，q次询问，四种操作
op1：x到y最短路径上所有节点的值加上val
op2：输出x到y最短路径上所有节点的权值和 
op3：将以x为根节点的子树内所有节点值都加上val 
op4：输出以x为根节点的子树内所有节点的权值和 
op5：输入x,y 输出lca(x,y) 
*/
const int N=5e5+10;
vector<int> e[N];
//如果开vector 有爆空间的可能 
int fa[N],son[N],sz[N],dep[N],dfn[N],rkn[N],top[N];  
int fi[N],la[N];
int idx=0;
void dfs1(int id,int u){  
	fa[id]=u;
	sz[id]=1;
	dep[id]=dep[u]+1;
	for(auto x:e[id]){
		if(x==u) continue;
		dfs1(x,id);
		sz[id]+=sz[x];
		if(sz[x]>sz[son[id]]) son[id]=x;
	}
}
void dfs2(int id,int tp){
	top[id]=tp;
	dfn[id]=++idx;
	fi[id]=la[id]=idx;
	rkn[idx]=id;
	if(son[id]) dfs2(son[id],tp);
	for(auto x:e[id]){
		if(x==fa[id]||x==son[id]) continue;
		dfs2(x,x);
	}
	for(auto x:e[id]){
		if(x==fa[id]) continue;
		la[id]=max(la[id],la[x]);
	}
}
vector<int> a(N);
#define lc p<<1
#define rc p<<1|1
struct node{
	int l,r,sum,mx;
	int lz;
}tr[4*N];
void pushup(int p){
	tr[p].sum=tr[lc].sum+tr[rc].sum;
	tr[p].mx=max(tr[lc].mx,tr[rc].mx);
}
void build(int p,int l,int r){  // (1,1,n)
	tr[p].l=l;tr[p].r=r;tr[p].lz=0;
	if(l==r){
		tr[p].sum=tr[p].mx=a[rkn[l]];
		return;
	}
	int m=l+r>>1;
	build(lc,l,m);
	build(rc,m+1,r);
	pushup(p);
}
void pushdown(int p){  
	if(tr[p].lz){
		tr[lc].mx+=tr[p].lz;
		tr[rc].mx+=tr[p].lz;
		tr[lc].sum+=(tr[lc].r-tr[lc].l+1)*tr[p].lz;
		tr[rc].sum+=(tr[rc].r-tr[rc].l+1)*tr[p].lz;
		tr[lc].lz+=tr[p].lz;
		tr[rc].lz+=tr[p].lz;
		tr[p].lz=0;
	}
}
void update1(int p,int x,int val){  //单点修改  传(1，dfn[id],val) 
	if(tr[p].l==tr[p].r){ 
		// tr[p].mx=tr[p].sum=val;
		tr[p].mx+=val;
		tr[p].sum+=(tr[p].r-tr[p].l+1)*val;
		return; 
	}
	pushdown(p);
	int m=tr[p].l+tr[p].r>>1;
	if(x<=m) update1(lc,x,val);
	else update1(rc,x,val);
	pushup(p);
}
void update2(int p,int l,int r,int val){  //区间修改  传(1,dfn[x],dfn[y],val) 
	if(l<=tr[p].l&&tr[p].r<=r){
		tr[p].sum+=(tr[p].r-tr[p].l+1)*val;
		tr[p].mx+=val;
		tr[p].lz+=val;
		return;
	}
	pushdown(p);
	int m=tr[p].l+tr[p].r>>1;
	if(l<=m) update2(lc,l,r,val);
	if(m<r) update2(rc,l,r,val);
	pushup(p);
}
int query1(int p,int l,int r){  // 查询sum  传(1,dfn[x],dfn[y]) 
	if(l<=tr[p].l&&tr[p].r<=r) return tr[p].sum;
	pushdown(p);
	int sum=0;
	int m=tr[p].l+tr[p].r>>1;
	if(l<=m) sum+=query1(lc,l,r);
	if(m<r) sum+=query1(rc,l,r);
	return sum;
}
int query2(int p,int l,int r){  //查询mx 传(1,dfn[x],dfn[y]) 
	if(l<=tr[p].l&&tr[p].r<=r) return tr[p].mx;
	pushdown(p);
	int mx=-1e9;
	int m=tr[p].l+tr[p].r>>1;
	if(l<=m) mx=max(mx,query2(lc,l,r));
	if(m<r) mx=max(mx,query2(rc,l,r));
	return mx;
}
int query_sum(int x,int y){  //查询两点路径的sum  传(x,y) 
	int ans=0,fx=top[x],fy=top[y];
	while(fx!=fy){
		if(dep[fx]>=dep[fy]) ans+=query1(1,dfn[fx],dfn[x]),x=fa[fx];
		else ans+=query1(1,dfn[fy],dfn[y]),y=fa[fy];
		fx=top[x];fy=top[y];
	}
	if(dfn[x]<dfn[y]) ans+=query1(1,dfn[x],dfn[y]);
	else ans+=query1(1,dfn[y],dfn[x]);
	return ans;
}
int query_mx(int x,int y){  //查询两点路径的mx  传(x,y) 
	int ans=-1e9,fx=top[x],fy=top[y];
	while(fx!=fy){
		if(dep[fx]>=dep[fy]) ans=max(ans,query2(1,dfn[fx],dfn[x])),x=fa[fx];
		else ans=max(ans,query2(1,dfn[fy],dfn[y])),y=fa[fy];
		fx=top[x];fy=top[y];
	}
	if(dfn[x]<dfn[y]) ans=max(ans,query2(1,dfn[x],dfn[y]));
	else ans=max(ans,query2(1,dfn[y],dfn[x]));
	return ans;
}
void update3(int x,int y,int val){  //区间更新两点路径的值  传(x,y) 
	int fx=top[x],fy=top[y];
	while(fx!=fy){
		if(dep[fx]>=dep[fy]) update2(1,dfn[fx],dfn[x],val),x=fa[fx];
		else update2(1,dfn[fy],dfn[y],val),y=fa[fy];
		fx=top[x];fy=top[y];
	}
	if(dfn[x]<dfn[y]) update2(1,dfn[x],dfn[y],val);
	else update2(1,dfn[y],dfn[x],val);
}
int lca(int x,int y){  //最近公共祖先 
	while(top[x]!=top[y]){
		if(dep[top[x]]>dep[top[y]]) x=fa[top[x]];
		else y=fa[top[y]];
	}
	return dep[x]>dep[y] ? y : x; 
}
void solve() {
	int n,q,root;
	cin>>n>>q>>root;
	for(int i=1;i<=n;i++) cin>>a[i];
	for(int i=1;i<n;i++){
		int u,v;
		cin>>u>>v;
		e[u].push_back(v);
		e[v].push_back(u);
	}
	dfs1(root,0);
	dfs2(root,root);
	build(1,1,n);
	while(q--){
		int op;cin>>op;
		if(op==1){
			int x,y,val;
			cin>>x>>y>>val;
			update3(x,y,val);
		}
		else if(op==2){
			int x,y;
			cin>>x>>y;
			cout<<query_sum(x,y)%mod<<"\n";
		}
		else if(op==3){
			int x,val;
			cin>>x>>val;
			update2(1,fi[x],la[x],val);
		}
		else if(op==4){
			int x;
			cin>>x;
			cout<<query1(1,fi[x],la[x])%mod<<"\n";
		}
		else if(op==5){
			int x,y;
			cin>>x>>y;
			cout<<lca(x,y)<<"\n";
		}else cout<<"FUCK YOU!\n";
	}
}
