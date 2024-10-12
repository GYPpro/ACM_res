#include<bitset>
int n;
auto bin = bitset<32>(n);

cout << bin;
//输出二进制位

cout << bin.to_ullong();
//输出十进制

bin[1] = 1
//随机访问

bin = !bin ^ (bin & bin | bitset<32>(1))
//位运算

bin != bin
//比较运算符

bin.count()
//1的数量

bin.test(i)
//随机访问，类似std::vector::pos()

bin.any()
//有一位1就true

bin.none()
//全0就返回true

bin.all()
//全1就返回true

bin.flip()
//翻转全部

bin.flip(i)
//a[i] = !a[i]

bin._Find_first()
//第一个1的下标

bin._Find_next(i)
//从下标n往后第一个1的下标