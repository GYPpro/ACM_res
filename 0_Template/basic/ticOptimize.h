//时间戳优化：对付多组数据很常见的技巧。 
int tag[N], t[N], Tag;
int lowbit(int x){
	return x&-x;
} 
void reset(){
	++Tag;
}
void add(int x,int val){
	while(x<=n){
		if(tag[x]!=Tag) t[x]=0;
		t[x]+=val;tag[x]=Tag;
		x+=lowbit(x);
	}
}
int getsum(int x){
	int ans=0;
	while(x){
		if(tag[x]==Tag) ans+=t[x];
		x-=lowbit(x);
	}
	return ans;
}
