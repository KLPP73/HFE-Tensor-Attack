

import itertools
import random
pubk=load("public_key.sobj")
Fq=pubk[0][0][0].base_ring()
q=Fq.order()
# ~ n1=len(pubk[0])
# ~ n2=len(pubk[0][0])
# ~ n3=len(pubk)


# ~ CorrectOrderTensor=[[[pubk[k][i][j] for k in range(n3)] for j in range(n2)] for i in range(n1)]

n1=len(pubk)
n2=len(pubk[0])
n3=len(pubk[0][0])

CorrectOrderTensor=[[[pubk[i][j][k] for k in range(n3)] for j in range(n2)] for i in range(n1)]

print(n1,n2,n3)

R=PolynomialRing(GF(q),'lambda',n1)
v=R.gens()

Mt=[[0 for i in range(n3)] for j in range(n2)]

for i in range(n1):
	for j in range(n2):
		for k in range(n3):
			Mt[j][k]+=v[i]*CorrectOrderTensor[i][j][k]

Mt=Matrix(Mt)

d=4
I=[]
for _ in range(pow(n1,d)):
	S1=random.sample(range(n2),d)
	S2=random.sample(range(n3),d)
	I.append(Mt[S1,S2].det())

I=Ideal(I+[v[0]-1])

g=I.groebner_basis()

print(g)
print(len(list(g)))

exit(0)
V=I.variety()

print(V)
print(len(V))

for a in V:
	print("R",Mt.subs(a).rank())
