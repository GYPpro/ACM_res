int dfn[N];
int h[N],m,a[N],len;
bool cmp(int x,int y){
	return dfn[x]<dfn[y];
}
void build(){
	sort(h+1,h+1+m,cmp);  //dfn排序 
	for(int i=1;i<m;i++){
		a[++len]=h[i];
		a[++len]=lca(h[i],h[i+1]);  //插入LCA 
	}
	a[++len]=h[m];
	sort(a+1,a+1+len,cmp);  //DFN排序 
	len=unique(a+1,a+len+1)-a+1;  //去重 
	for(int i=1,lc;i<len;i++){
		lc=lca(a[i],a[i+1]);
		conn(lc,a[i+1]);  //连边 
	}
}
