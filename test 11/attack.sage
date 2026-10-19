import random

load('library/tensor_attack.sage')

pubk=load("public_key.sobj")


n=pubk.nrows()

pubk=[[pubk[i][j] for j in range(pubk.ncols())] for i in range(pubk.nrows())]
for i in range(n):
	for j in range(i,n):
		pubk[j][i]+=pubk[i][j]

pubk=Matrix(pubk)
Fe=pubk.base_ring()

tensor=[[pubk[i][j].list() for j in range(pubk.ncols())] for i in range(pubk.nrows())]


Fq=tensor[0][0][0].base_ring()


tensor=switch_dimensions(tensor,0,1)
n1,n2,n3=dimensions(tensor)
# ~ tensor=tensor[n2:]

print(dimensions(tensor))
n1,n2,n3=dimensions(tensor)

n1=30

q=len(tensor[0][0][0].base_ring())
print(q)

R=PolynomialRing(GF(q**n2),'lambda',n1)
v=R.gens()

Mt=[[0 for i in range(n3)] for j in range(n2)]

for i in range(n1):
	for j in range(n2):
		for k in range(n3):
			Mt[j][k]+=v[i]*tensor[i][j][k]

Mt=Matrix(Mt)

d=2
I=[]
for _ in range(pow(n1,d)):
	S1=random.sample(range(n2),d)
	S2=random.sample(range(n3),d)
	I.append(Mt[S1,S2].det())

I=Ideal(I+[v[0]-1])

g=I.groebner_basis()

print(g)
print(len(list(g)))
V=I.variety()

print(V)
print(len(V))
