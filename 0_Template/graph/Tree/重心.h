int sz[N];  //�ڵ��С 
int weight[N]; //���������С 
int centroid[2]; //��������
void dfs(int id,int fa){
	sz[id]=1;
	weight[id]=0;
	for(auto x:e[id]){
		if(x==fa) continue;
		dfs(x,id);
		sz[id]+=sz[x];
		weight[id]=max(weight[id],sz[x]);
	}
	weight[id]=max(weight[id],n-sz[id]);
	if(weight[id]<=n/2){
		centroid[centroid[0]!=0]=id;
	}
} 
