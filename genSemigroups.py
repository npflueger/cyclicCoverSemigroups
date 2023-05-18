####
#### This file contains various functions to generate the semigroups
#### considered in "Weierstrass semigroups from cyclic covers of hyperelliptic curves."
####


"""
Class to represent a numerical semigroup and its various statistics.
Semigroups are created using their standard basis e_0, ..., e_(m-1),
where m = multiplicity = minimum positive element of S, and
e_i = minimum positive element of S in i + m Z.
"""
class Semigroup():
    """
    Initialize a semigroup from a "standard basis" e, where e[0] = multiplicity.
    The optional "memo" argument allows the user to add a string with any data 
    about how the semigroup was constructed.
    """
    def __init__(self, e, memo=None):
        m = len(e)
        assert(e[0] == m)
        g = sum([(e[i]-1) // m for i in range(m)]) 
        gaps = [n for n in range(1,2*g) if n < e[n % m] ]
        assert(g == len(gaps))
        wt = sum(gaps) - sum(range(1,g+1))
        
        gens = []
        for n in range(1,2*g+m):
            if n in gaps: continue
            isGen = True
            for a in gens:
                if n - a not in gaps:
                    isGen = False
            if isGen:
                gens += [n]

        ewt = 0
        for a in gens:
            for b in gaps:
                if a < b: ewt += 1
        self.e = e
        self.m = m
        self.g = g
        self.gaps = gaps
        self.gens = gens
        self.wt = wt
        self.ewt = ewt
        self.memo = memo

    def showStats(self):
        if self.memo: print(self.memo)
        print("Semigroup statistics:")
        print("\te =",self.e)
        print(f"\tg = {self.g}\n\twt = {self.wt}\n\tewt = {self.ewt} = {self.ewt/self.g} g\n\tgens = {self.gens}")
        print("-"*40)



"""
Constructs a semigroup from the data N, gamma, d, and the 
"a-vector" defined in Question 1, i.e. (a_1, ..., a_{N-1}), where
a_j = ceiling of  min(e_{N-j}, e_{2N-j}) / N.
Note that a[0] is never accessed, since it is not defined in the paper's notation.
"""
def sgFromAVector(N, gamma, d, a, memo = None):
    if memo == None: memo = f"Constructing a semigroup from these data:\n\tN={N}    gamma={gamma}    d={d}"
    memo += f"\n\ta-vector        a = {a[1:]} (beginning with a_1)"
    m = [None] + [a[j]-j*d for j in range(1,N)]
    memo += f"\n\tMult. profile   m = {m[1:]} (beginning with m_1)"
    g = (d * (N-1) * N + 2*N*(gamma-1) + 2) // 2
    e = [None]*(2*N)
    e[0] = 2*N
    e[N] = N*(2*gamma+1)
    badIneq = False
    for j in range(1,N):
        # min of e_{N-j} and e_{2N-j}
        minE = N*a[j] - j
        sumE = 2*j*N*d + (2*gamma+1)*N - 2*j
        maxE = sumE - minE
        assert(minE <= maxE)

        if a[j] % 2 == 0:
            e[N-j], e[2*N-j] = maxE, minE
        else:
            e[N-j], e[2*N-j] = minE, maxE

    return Semigroup(e,memo)

"""
Constructs a semigroup from data N, gamma, d, and a given
multiplication profile m_1, ..., m_{N-1}. (m_0 never accessed)
"""
def sgFromMP(N, gamma, d, m, memo = None):
    a = [None]*N # Initialize
    for j in range(1,N):
        a[j] = m[j] + j*d
    # Note: a[0] = None still, but a_0 is undefined so this is never accessed.
    return sgFromAVector(N,gamma,d,a,memo)

"""
Constructs a semigroup S_N from the "staircase" multiplication profile
of given N, d, and mult. profile as specified in Example 5.11.
"""
def staircaseSG(N, gamma, d):
    # NOTE: this m begins with m[0] = 0, although the paper writes mult. profiles
    # beginning at 1. this doesn't matter, since m_0 is undefined.
    m = list(range(gamma)) + (N-2*gamma)*[gamma] + list(range(gamma-1,-1,-1)) 
    memo = f"Staircase semigroup with data:\n\tN = {N}    gamma = {gamma}    d = {d}"
    return sgFromMP(N,gamma,d,m,memo)


       
"""
Constructs a semigroup from its generators.
Doesn't try to do this too cleverly or efficiently.
Method is: increase a variable n one at at time, and at each step check
if it belongs to the semigroup. At the same time, maintain a partial 
list of the standard basis, including only elements less than n (with
elements set to None otherwise). Using this partial list, and the 
generators, check whether n is in the semigroup, and add it to the 
standard basis if its the first in its class mod m. Stop when the full
standard basis has been seen.
"""
def semigroupFromGens(gens, memo = None):
    gens.sort()
    if memo == None: memo = f"Constructing a semigroup from generators {gens}"
    m = gens[0]
    e = [None]*m
    e[0] = m
    remainingE = m-1
    n = 0
    while remainingE > 0:
        n += 1
        # Check whether or not n is in the semigroup and minimal in its class
        if e[n%m] != None: continue
        inS = False
        for a in gens:
            b = n - a
            r = b % m
            if b==0 or (e[r] != None and e[r] <= b):
                inS = True
        if inS:
            e[n%m] = n
            remainingE -= 1
    return Semigroup(e, memo)
        
"""
Generator for semigroups with N=3 and given gamma,d.
Proceeds by finding all a-vectors in the feasible set.
"""
def nThreeSemigroups(N, gamma, d):
    N = 3
    g = (d * (N-1) * N + 2*N*(gamma-1) + 2) // 2
    assert(N == 3) # For now, just consider this case
    for a in feasibleSet(N,gamma,d):
        if a[2]%2 == 0 or 4*d+2*gamma-2*a[1]+1 <= a[2] <= 2*a[1]:
            yield sgFromAVector(N,gamma,d,a)

"""
Generator for feasible set
For convenience in looking up entries, the yielded lists are
[None, a_1, a_2, ..., a_{N-1}] so that a_i is a[i].
"""
def feasibleSet(N,gamma,d):
    for a in feasibleSetTilde(N,gamma,d):
        if feasibleIneqs(a,N):
            yield a
def feasibleSetTilde(N,gamma, d):
    if N == 1:
        yield [None]
        return
    for a in feasibleSet(N-1,gamma,d):
        for ai in range((N-1)*d, (N-1)*d+gamma+1):
            yield a + [ai]
def feasibleIneqs(a,N):
    for i in range(1,N):
        for j in range(1,N-i):
            if a[i+j] > a[i] + a[j]:
                return False
    return True


"""
The Carvalho-Torres semigroup 
"""
def carvalhoTorres(N,gamma,d):
    memo = f"Carvalho-Torres semigroup for data\n\tN = {N}    gamma = {gamma}    d = {d}"
    gens = [2*N, (2*gamma+1)*N, N*d-1]
    return semigroupFromGens(gens,memo)

"""
Generate all semigroups produced by our construction,
for N = 3 and given genus.
Return value is number of semigroups found
"""
def listByGenus(g):
    count = 0
    N = 3
    for gamma in range(1, g // (2*N-1) + 1):
        d = 0
        while True:
            d += 1
            thisg = (d * (N-1) * N + 2*N*(gamma-1) + 2) // 2
            if thisg > g: break
            if thisg < g: continue
            for S in nThreeSemigroups(N, gamma, d):
                S.showStats()
                print("")
                count += 1
    return count
