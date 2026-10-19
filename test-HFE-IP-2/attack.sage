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


tensor=switch_dimensions(tensor,0,2)
n1,n2,n3=dimensions(tensor)
# ~ tensor=tensor[n2:]

print(dimensions(tensor))
n1,n2,n3=dimensions(tensor)

p=4
print(compute_stats(tensor,p,2*p+1))

# ~ print(construct_attack(tensor,p,2*p+1))
