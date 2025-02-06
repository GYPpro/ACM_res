/*
n座城市,m条双向边,每条路有限重,问从x到y最多能运输多重的货物 
*/
void solve(){
	int n,m;
	cin>>n>>m;
	vector<int> e[n*2+1];
	priority_queue<array<int,3>> q;
	for(int i=1;i<=m;i++){
		int u,v,w;
		cin>>u>>v>>w;
		q.push({w,u,v});
	}
	vector<int> fa(n*2+10),w(n*2+10);
	iota(fa.begin(),fa.end(),0);
	function<int(int)> find=[&](int x){
		return fa[x]==x ? x : fa[x]=find(fa[x]);
	};
	auto merge=[&](int x,int y){
		int fx=find(x),fy=find(y);
		fa[fx]=fy;
	};
	int cnt;
	auto Kruskal=[&]()->void {
		cnt=n;
		while(!q.empty()){
			auto a=q.top();q.pop();
			int fx=find(a[1]),fy=find(a[2]);
			if(fx==fy) continue;
			cnt++;
			merge(fx,cnt);
			merge(fy,cnt);
			w[cnt]=a[0];
			e[cnt].push_back(fy);
			e[cnt].push_back(fx);
		}
	};
	Kruskal();
	vector<vector<int>> faa(n*2+10,vector<int>(19));
	vector<int> dep(n*2+10);
	vector<bool> vis(n*2+10);
	function<void(int,int)> dfs=[&](int id,int u){
		faa[id][0]=u;dep[id]=dep[u]+1;vis[id]=1;
		for(int i=1;i<=19;i++){
			faa[id][i]=faa[faa[id][i-1]][i-1];
		}
		for(auto x:e[id]) dfs(x,id);
	};
	for(int i=cnt;i>=1;i--){
		if(!vis[i]) dfs(i,i);
	}
	auto lca=[&](int x,int y)->int{
		if(find(x)!=find(y)) return 0;
		if(dep[x]<dep[y]) swap(x,y);
		int tmp=dep[x]-dep[y];
		for(int i=0;i<19;i++){
			if((tmp>>i&1)) x=faa[x][i];
		}
		if(x==y) return x;
		for(int i=18;i>=0;i--){
			if(faa[x][i]!=faa[y][i]){
				x=faa[x][i];
				y=faa[y][i];
			}
		}
		return faa[x][0];
	};
	w[0]=-1;
	int qq;cin>>qq;
	while(qq--){
		int x,y;
		cin>>x>>y;
		int l=lca(x,y);
		cout<<w[l]<<"\n";
	}
}

