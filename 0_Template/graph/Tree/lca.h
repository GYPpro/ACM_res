
class LCA{ 
public:
    vector<vector<pii>> cnj;
    vector<int> lg,dep;
    vector<array<int,32>> fa,wei;
    int n;

    LCA(int _n) {
        n = _n;
        cnj.resize(n+1);
        lg.resize(n+1),fa.resize(n+1),dep.resize(n+1),wei.resize(n+1);
        for(int i = 1; i <= n; i ++)
            lg[i] = lg[i-1] + (1 << lg[i-1] == i);
    }

    void addEdge(int u,int v,int w) {
        cnj[u].push_back({v,w});
        cnj[v].push_back({u,w});
    }

    void build(int rt = 1) {
        using itf = function<void(int,int)>;
        itf dfs = [&](int p,int f) -> void {
            fa[p][0] = f,dep[p] = dep[f] + 1;
            // wei[p][0] = 0;
            for(int i = 1;i <= lg[dep[p]];i ++) fa[p][i] = fa[fa[p][i-1]][i-1];
            for(int i = 1;i <= lg[dep[p]];i ++) wei[p][i] = max(wei[p][i-1],wei[fa[p][i-1]][i-1]);
            for(auto [x,w]:cnj[p]) if(x == f) continue;
            else wei[x][0] = w,dfs(x,p);
        };
        dfs(rt,0);
    }

    int get(int x,int y) {
        if(dep[x] < dep[y]) swap(x,y);
        while(dep[x] > dep[y]) x = fa[x][lg[dep[x] - dep[y]] - 1];
        if(x == y) return x;
        for(int k = lg[dep[x]]-1;k >= 0;k --) if(fa[x][k] != fa[y][k]) x = fa[x][k],y = fa[y][k];
        return fa[x][0];
    }

    int getmaxw(int x,int y) {
        int curmx = 0;
        if(dep[x] < dep[y]) swap(x,y);
        while(dep[x] > dep[y]) curmx = max(curmx,wei[x][lg[dep[x] - dep[y]] - 1]), x = fa[x][lg[dep[x] - dep[y]] - 1];
        if(x == y) return curmx;
        for(int k = lg[dep[x]]-1;k >= 0;k --) 
            if(fa[x][k] != fa[y][k]) 
                curmx = max(curmx,wei[x][k]),x = fa[x][k],
                curmx = max(curmx,wei[y][k]),y = fa[y][k];
        return max({curmx,wei[x][0],wei[y][0]});
    } 
};

//----------VERSION 2-----------


void solve(){
	int n,m,root;
	cin>>n>>m>>root;
	vector<int> e[n+1],dep(n+1);
	vector<vector<int>> fa(n+1,vector<int>(21));
	for(int i=1;i<n;i++){
		int u,v;
		cin>>u>>v;
		e[u].push_back(v);
		e[v].push_back(u);
	}
	function<void(int,int)> dfs=[&](int id,int u){
		fa[id][0]=u;
		dep[id]=dep[u]+1;
		for(int i=1;i<=20;i++) fa[id][i]=fa[fa[id][i-1]][i-1];
		for(auto x:e[id]){
			if(x==u) continue;
			dfs(x,id);
		}
	};
	dfs(root,root); 
	function<int(int,int)> lca=[&](int x,int y){
		if(dep[x]<dep[y]) swap(x,y);
		int tmp=dep[x]-dep[y];
		for(int i=0;i<=20;i++){
			if(tmp>>i&1) x=fa[x][i];
		}
		if(x==y) return x;
		for(int i=20;i>=0;i--){
			if(fa[x][i]!=fa[y][i]){
				x=fa[x][i];
				y=fa[y][i];
			}
		}
		return fa[x][0];
	};
	while(m--){
		int a,b;
		cin>>a>>b;
		cout<<lca(a,b)<<"\n";
	}
	return;
} 
