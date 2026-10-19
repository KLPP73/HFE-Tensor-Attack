import secrets
import random

q=65537
n=10

# ~ q=q**n

F=GF(q)

V1=[(F**n).random_element() for i in range(n)]
V2=[(F**n).random_element() for i in range(n)]


# ~ St=random.sample(range(n**2),n*u)
# ~ St=[(a%n,a//n) for a in St]

# ~ St=[(i,j) for j in range(n) for i in range(n)]
# ~ k=6
# ~ St=[(i,j) for (i,j) in St if max(i,j)>=k]

# ~ Mt=[[int((i,j) in St) for j in range(n)] for i in range(n)]

# ~ print(Matrix(Mt))


u=4
while True:
	pi=[Permutations(n).random_element() for i in range(u)]
	W=True
	for i in range(n):
		S=[pi[j][i] for j in range(u)]
		if len(set(S))!=u:
			W=False
	if not W:
		continue
	St=[(i,p[i]-1) for i in range(n) for p in pi]
	break


print(St)
tensor=[]

for i in range(2*n):
	Mt=Matrix(F,n,n)
	
	for (j,k) in St:
		Mt+=Matrix(V1[j]).transpose()*Matrix(V2[k])*F.random_element()
	
	# ~ for j in range(n):
		# ~ for k in range(n):
			# ~ Mt+=Matrix(V1[j]).transpose()*Matrix(V2[k])*F.random_element()
	Mt=[[Mt[i][j] for j in range(n)] for i in range(n)]
	
	tensor.append(Mt)
	

# ~ tensor=[[[F.random_element() for k in range(n)] for j in range(n)] for i in range(2*n)]


save(tensor,"public_key")
