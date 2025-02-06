// sz =int(sqrt(n)) 注意考虑sqrtl(n) 
struct node{
	int l,r,id;
	bool operator<(const node &x) const{
		if(l/sz!=x.l/sz) return l<x.l;
		if((l/sz)&1) return r<x.r;
		else return r>x.r;	
	}
};
