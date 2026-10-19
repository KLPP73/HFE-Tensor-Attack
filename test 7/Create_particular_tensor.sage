import secrets
import random
import itertools

q=65537
n1=10
n2=n1

# ~ q=q**n

F=GF(q)

# ~ V1=[(F**n).random_element() for i in range(n)]
# ~ V2=[(F**n).random_element() for i in range(n)]
V1=[identity_matrix(n1)[i] for i in range(n1)]
V2=[identity_matrix(n2)[i] for i in range(n2)]

# ~ St=random.sample(range(n**2),n*u)
# ~ St=[(a%n,a//n) for a in St]

# ~ St=[(i,j) for j in range(n) for i in range(n)]
# ~ k=6
# ~ St=[(i,j) for (i,j) in St if max(i,j)>=k]

# ~ Mt=[[int((i,j) in St) for j in range(n)] for i in range(n)]

# ~ print(Matrix(Mt))

# ~ St=[(i,j) for j in range(n) for i in range(n) if i!=j]
# ~ u=2
# ~ while True:
	# ~ pi=[Permutations(n1).random_element() for i in range(u)]
	# ~ W=True
	# ~ for i in range(n1):
		# ~ S=[pi[j][i] for j in range(u)]
		# ~ if len(set(S))!=u:
			# ~ W=False
	# ~ if not W:
		# ~ continue
	# ~ St=[(i,p[i]-1) for i in range(n1) for p in pi]
	# ~ break


# ~ print(St)

u=3
St=[]
r=random.sample(range(n1),u)
r=range(u)
for i in range(n1):
	St+=[((i+r[j])%n1,i) for j in range(u)]

# ~ u=3

# ~ St=[]
# ~ for i in range(n2):
	# ~ r=random.sample(range(n1),u)
	# ~ St+=[(a,i) for a in r]

# ~ print(St)

tensor=[]

for i in range(100):
	Mt=Matrix(F,n1,n2)
	
	for (j,k) in St:
		Mt+=Matrix(V1[j]).transpose()*Matrix(V2[k])*F.random_element()
	
	# ~ for j in range(n):
		# ~ for k in range(n):
			# ~ Mt+=Matrix(V1[j]).transpose()*Matrix(V2[k])*F.random_element()
	
	Mt=Mt.transpose()
	Mt=[[Mt[i][j] for j in range(Mt.ncols())] for i in range(Mt.nrows())]
	tensor.append(Mt)
	

# ~ tensor=[[[F.random_element() for k in range(n1)] for j in range(n1)] for i in range(2*n1)]


save(tensor,"public_key")
# ~ print(tensor)
# ~ tot=Matrix(F,n,n)
# ~ for per in Permutations(n):
	# ~ res=identity_matrix(F,n)
	# ~ for a in per:
		# ~ res*=tensor[a-1]
	# ~ tot+=res*per.sign()

# ~ print(tot)
