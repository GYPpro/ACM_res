vector<int> dep(n+1);
int mx=0,tmp=-1;
function<void(int,int)> dfs=[&](int id,int fa){
	dep[id]=dep[fa]+1;
	for(auto x:e[id]){
		if(x==fa) continue;
		dfs(x,id);
	}
	if(dep[id]>mx){
		mx=dep[id];
		tmp=id;
	}
};
dfs(1,0);
int d1=tmp;
mx=0;
dfs(d1,0);
int d2=tmp; //d1,d2
