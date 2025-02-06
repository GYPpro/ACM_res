from decimal import *
import sys
setcontext(Context(prec=2000000, Emax=2000000, Emin=0)) 
a = sys.stdin.readline()
b = sys.stdin.readline()
print(Decimal(a) * Decimal(b))