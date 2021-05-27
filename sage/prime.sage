


def prime(x):
	return 2*x^12 - 1

xs = [6334792777, 8176556533, 8426067021]
ps = [prime(x) for x in xs]
pms = [(p-1).factor() for p in ps]
pps = [(p+1).factor() for p in ps]

factors = [list(pms[i]*pps[i]) for i in [0..len(xs)-1]]

def qf_from_Gram(gram):
	n = gram.ncols()
	H = [0]*(n*(n+1)//2)
	k = 0
	for i in range(n):
		for j in range(i,n):
			if i == j:
				H[k] = gram[i][j]
			else:
				H[k] = 2*gram[i][j]
			k += 1
	return QuadraticForm(ZZ, n, H);

# constraints: 
# N1, N2 coprime and odd (or at least N2 odd)
# N2 > N1^2
# log N1 ~ 2 lambda = 256
# log N2 ~ 4 lambda = 512
# d = N2 mod N1^2
# D = (N2^2 - d^2)/N1^2
# D is square mod p
# -D not square mod every prime divisor of N1


def check(N1,N2,p,q):
	if (N2 < N1^2):
		print("Failure: the constraint N2 > N1^2 is not satisfied")
		return false
	if GCD(N1,N2) != 1:
		print("Failure: N1 and N2 have a common factor")
		return false
	if N2 % 2 == 0:
		print("Failure: N2 is even")
		return false
	if (p^2-1) % (N1*N2) != 0:
		print("Failure: N1*N2 does not divide p^2-1")
		return false

	d = N2 % (N1^2)
	D1 = (N2 - d)
	D2 = (N2 + d)
	D = D1*D2 // N1^2
	D1 = D // GCD(D1,D);
	D2 = D // D1;


	if not Zmod(p)(q*D).is_square(): 
		print("Failure: q*D is not a square mod p")
		return false

	if not Zmod(q)(p*D).is_square(): 
		print("Failure: p*D is not a square mod q")
		return false
	fac = N1.factor()
	for ell,n in fac:
		if Zmod(ell)(-D).is_square():
			print("Failure: -D is a square mod", ell)
			return false
	print("Success: all tests passed for the parameters")
	print("\t p =", p)
	print("\t N1 =", N1.factor())
	print("\t N2 =", N2.factor())
	return true

# p = prime(6334792777)
# N1 = 1066643^12
# N2 = 3^3 * 5 * 7 * 13 * 17 * 53 * 199 * 271 * 277 * 433 * 461 * 547 * 3373 * 5939^12 * 10181 * 16273 * 104173 * 600109 * 2229307 * 3254137 * 3290561 * 5790487 * 6109057 * 9918889 * 12011149 * 20650789 * 24685489

# check(N1,N2,p)

p = prime(8426067021)
# N1 = 43^12 * 84719^11
# N2 = 3^21 * 5^1 * 7^1 * 13^1 * 17^1 * 19^1 * 23^1 * 41^0 * 73^1 * 257^12 * 313^1 * 1009^1 * 2857^1 * 3733 * 5519 * 6961 * 53113 * 499957 * 763369 * 2101657 * 2616791 * 7045009 * 11959093 * 17499277 * 20157451 * 33475999 * 39617833 *45932333


# BEST CANDIDATE

N1 = 43^12 * 84719^10
N2 = 3^22 * 5^1 * 7^1 * 13^1 * 17^1 * 19^1 * 23^1 * 41^1 * 73^1 * 257^12 * 313^1 * 1009^1 * 2857^1 * 3733 * 5519 * 6961 * 53113 * 499957 * 763369 * 2101657 * 2616791 * 7045009 * 11959093 * 17499277 * 20157451 * 33475999 * 39617833 *30298249


# SMALLER VALUE GOOD FOR TESTS

# N1 = 43^11 * 84719^0
# N2 = 3^23 * 5^1 * 7^1 * 13^1 * 17^1 * 19^1 * 23^1 * 41^1 * 73^1 * 257^12 * 313^1 * 1009^1 * 2857^1 * 3733 * 5519 * 6961 * 53113



check(N1,N2,p,3)


