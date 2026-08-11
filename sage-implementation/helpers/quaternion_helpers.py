import itertools
from sage.all import *
from random import randint

def ideal_generator(I):
    p = I.quaternion_algebra().ramified_primes()[0]
    bound = max(ceil(10*log(p,2)), 100)
    while True:
        alpha = sum([b * randint(-bound, bound) for b in I.basis()])
        if gcd(alpha.reduced_norm(), I.norm()**2) == I.norm():
            return alpha
        
def quaternion_order_basis(alpha, O):
    assert alpha in O
    Obasis = O.basis()
    M_O = Matrix(QQ, [b.coefficient_tuple() for b in Obasis]).transpose()
    vec_alpha = vector(alpha.coefficient_tuple())
    coeffs = M_O.inverse() * vec_alpha
    coeffs = [ZZ(c) for c in coeffs]
    #assert alpha == sum([c * beta for c, beta in zip(coeffs, Obasis)])
    return coeffs

def make_primitive(alpha, O):
    v_alpha = quaternion_order_basis(alpha, O)
    d = gcd(v_alpha)
    assert alpha/d in O
    return alpha/d
        
def cornacchia(QF, m):
    m_prime = prod([l**e for l, e in factor(m, limit=1000) if l < 1000])
    if not is_pseudoprime(m/m_prime):
        return None, None, False
    sol = QF.solve_integer(m)
    if not sol:
        return None, None, False
    return sol[0], sol[1], True

def represent_integer(O0, M):
    B = O0.quaternion_algebra()
    QF = BinaryQF([1,0,1])
    i,j,k = B.gens()
    p = B.ramified_primes()[0]
    sol = None
    for _ in range(10000):
        m1 = floor(sqrt(round((4*M)/p, 5)))
        z = randint(-m1, m1)
        m2 = floor(sqrt(round((4*M-z**2)/p, 5)))
        w = randint(-m2, m2)
        Mm = M - p*QF(z,w)
        x, y, found = cornacchia(QF, Mm)
        if not found:
            continue
        else:
            gamma = x + y*i + j*z + k*w
            gamma = make_primitive(gamma,O0)
            if gamma.reduced_norm() == M:
                return gamma
    assert False

def heuristic_random_ideal(O0, M):
    p = O0.quaternion_algebra().ramified_primes()[0]
    N = next_prime(randint(p, 2*p))
    while gcd(N, M) > 1:
        N = next_prime(randint(p, p**2))
    alpha = represent_integer(O0, N*M)
    I = O0*alpha + O0*M
    assert I.norm() == M
    return I

def order_isomorphism(O, alpha):
    B = O.quaternion_algebra()
    return B.quaternion_order([alpha * g * ~alpha for g in O.basis()])

def randomise_order(O):
    B = O.quaternion_algebra()
    p = B.ramified_primes()[0]
    while True:
        r = sum(randrange(-p,p)*g for g in B.basis())
        if r:
            break
    return order_isomorphism(O, r), r

def connecting_ideal(O1, O2):
    I = O1*O2
    I *= I.norm().denominator()
    return I

def gram_matrix(basis):
    M = []
    for a in basis:
        M.append([QQ(2) * a.pair(b) for b in basis])
    return Matrix(QQ, M)

def reduced_basis(I):
    B = I.basis()
    G = gram_matrix(I.basis())
    U = G.LLL_gram().transpose()
    return [sum(c*beta for c, beta in zip(row, B)) for row in U]

def reduced_ideal(I, return_elt=False):
    reduced_basis_elements = reduced_basis(I)
    beta = reduced_basis_elements[0]
    J = I*(beta.conjugate()/I.norm())
    # assert J.conjugate().is_equivalent(I.conjugate())
    if return_elt:
        return J, beta
    return J

def reduced_ideal_odd(I):
    reduced_basis_elements = reduced_basis(I)
    c = 0
    beta = reduced_basis_elements[c]
    while beta.reduced_norm()/I.norm() % 2 == 0:
        c += 1
        if c > 10000:
            assert False
        if c > 3:
            beta = sum(randint(1,10)*gen for gen in reduced_basis_elements)
        else:
            beta = reduced_basis_elements[c]

    J = I*(beta.conjugate()/I.norm())
    assert J.conjugate().is_equivalent(I.conjugate())
    return J

def reduced_ideal_prime(I):
    reduced_basis_elements = reduced_basis(I)
    c = 0
    beta = reduced_basis_elements[c]
    while not is_pseudoprime(beta.reduced_norm()/I.norm()):
        c += 1
        if c > 10000:
            assert False
        if c > 3:
            beta = sum(randint(1,10)*gen for gen in reduced_basis_elements)
        else:
            beta = reduced_basis_elements[c]

    J = I*(beta.conjugate()/I.norm())
    assert J.conjugate().is_equivalent(I.conjugate())
    return J

def reduced_basis_from_gens(gens, Quat):   
    """
    More generally reduces the basis of any (not necessarily full rank) lattice
    """
    M = matrix(list(g) for g in gens)
    big = Integer(10)**99
    scal = diagonal_matrix(ZZ, [ceil(sqrt(g.reduced_norm())*big) for g in Quat.basis()])
    M = (M * scal).LLL() * ~scal
    return [sum(c*g for c,g in zip(v, Quat.basis())) for v in M]

def _roots_mod_prime_pow(d, l, e):
    """ This is extremely slow with sage, so do it custom """
    if not is_square(Integers(l**e)(d)):
        return []
    if d % l**e == 0:
        return [ZZ(r) for r in Integers(l**e)(0).sqrt(all=True)]

    k_start = 1
    while d % l**k_start == 0:
        k_start += 1
    assert k_start <= e, "Replace this hacky code with sage 10.3 upgrade"
    rts = Integers(l**k_start)(d).sqrt(all=True)
    rts = [ZZ(r) for r in rts]
    if l < 1000:
        for k in range(k_start+1,e+1):
            lift_rts = []
            for r in rts:
                Z_le = Integers(l**k)
                for m in range(l):
                    if Z_le(r + m*l**(k-1))**2 == Z_le(d):
                        lift_rts.append(r + m*l**(k-1))
            rts = list(set(lift_rts))
            if len(rts) == 0:
                return rts
    else:
        assert e <= 2, "e was too small"
        lift_rts = []
        if e == 2:
            for r in rts:
                x1 = (((d-r**2)/l)*pow(2*r, -1, l)) % l

                lift_r = (r + x1*l) % l**2
                assert Integers(l**2)(lift_r)**2 == Integers(l**2)(d), "lift was incorrect"
                lift_rts.append(lift_r)
            rts = lift_rts
    return rts
            

def all_roots_mod(d, M):
    fac_M = factor(M)
    
    all_rts = []
    mods = []
    for l, e in fac_M:
        #_, X = PolynomialRing(Integers(l**e), "X").objgen()
        rts = _roots_mod_prime_pow(d, l, e)
        #print(f"for l, e: {l, e}")
        #print(f"     > {sorted((X**2 - d).roots(multiplicities=False)) == sorted(rts)}")
        if len(rts) == 0:
            return []
        all_rts.append(rts)
        mods.append(l**e)

    sols = []
    for crt_in in itertools.product(*all_rts): 
        sols.append(crt(list(crt_in), mods))
    
    return sols

def all_cornacchia(d, m, only_mod_2 = False):
    # Find all solutions to x^2 + dy^2 = m
    if m < 0:
        return []
    gfullist = [[l**k for k in range((e//2) + 1)] for l, e in factor(m/m.squarefree_part())]
    usedgs = []
    sols = []
    for glist in itertools.product(*gfullist):
        if not glist:
            g = 1
        else:
            g = prod(glist)
        if g in usedgs:
            continue
        usedgs.append(g)
        tempm = ZZ(m/(g**2))
        #_, X = PolynomialRing(Integers(tempm), "X").objgen()
        rs = all_roots_mod(-ZZ(d), tempm)
        #print(sorted((X**2 + d).roots(multiplicities=False)) == sorted(rs))
        bound = round(tempm**(1/2),5)
        for r in rs:
            n = tempm
            while r > bound:
                n, r = r, n%r
            s = sqrt((tempm - r**2)/d)
            if s in ZZ:
                yield g*r, g*s
                if not only_mod_2:
                    yield g*r, -g*s
                    yield -g*r, g*s
                    yield g*r, -g*s

                if d == 1:
                    yield g*s, g*r
                    if not only_mod_2:
                        yield g*s, -g*r
                        yield -g*s, g*r
                        yield g*s, -g*r

def order_isomorphism(O, alpha):
    B = O.quaternion_algebra()
    return B.quaternion_order([alpha * g * ~alpha for g in O.basis()])

def reduced_trace_0_basis(I):
    T0 = (QQ**4).submodule((QQ**4).basis()[1:])
    M = I.free_module().intersection(T0).basis_matrix()
    return reduced_basis_from_gens([sum(c*g for c,g in zip(v, I.quaternion_algebra().basis())) for v in M], I.quaternion_algebra())


def make_cyclic(I):
    O = I.left_order()
    d = gcd([gcd(quaternion_order_basis(beta, O)) for beta in I.basis()])
    return I*(1/d)

def is_isomorphic(O1, O2):
    I = connecting_ideal(O1, O2)
    for alpha in reduced_basis(I):
        if alpha.reduced_norm() == I.norm():
            return True
    return False

def compute_isomorphism(O1, O2):
    I = connecting_ideal(O1, O2)
    for alpha in reduced_basis(I):
        if alpha.reduced_norm() == I.norm():
            return alpha
    assert False, "Not isomorphic"

def basis_matrix(O):
    M_O = Matrix(QQ, [ai.coefficient_tuple() for ai in O.gens()])
    return M_O

def successive_minima(O):
    return [b.reduced_norm() for b in reduced_basis(O)]

def pullback(I, J):
    r"""Compute the pullback of I by J"""
    assert I.left_order() == J.right_order()
    assert gcd(I.norm(), J.norm()) == 1

    O = J.left_order()
    return J*I + O*I.norm()

def pushforward(I, J):
    r"""Compute the pushforward of I by J"""
    assert I.left_order() == J.left_order()
    assert gcd(I.norm(), J.norm()) == 1
    
    O = J.right_order()
    return J.conjugate()*I + O*I.norm()



##############################
#
#       N-torsion stuff
#
##############################

class O_mod_N():
    def __init__(self, O, N):
        assert is_prime(N), "Not implemented for non-prime N :(("
        self.N = N
        self.B = O.quaternion_algebra()
        self.O = O
        self.basis_quat = self.B.basis()
        mat_i, mat_j, mat_k = self.B.modp_splitting_data(N)
        self.Z_N = Integers(N)
        mat_1 = Matrix(self.Z_N, [[1, 0], [0, 1]])
        self.mat_basis_quat = [mat_1, mat_i, mat_j, mat_k]

        self.basis_order = O.basis()
        self.mat_basis_order = []
        for beta in self.basis_order:
            self.mat_basis_order.append(sum(self.Z_N(c)*mat for c, mat in zip(beta.coefficient_tuple(), self.mat_basis_quat)))
        
    def project(self, alpha):
        assert alpha in self.O
        return sum(self.Z_N(c)*mat for c, mat in zip(alpha.coefficient_tuple(), self.mat_basis_quat))
    
    def lift(self, mat_alpha):
        system = []
        a_vec = []
        for ii in range(2):
            for jj in range(2):
                system.append([m[ii, jj] for m in self.mat_basis_order])
                a_vec.append(mat_alpha[ii, jj])
        system = Matrix(self.Z_N, system)
        a_vec = vector(self.Z_N, a_vec)
        v = system.solve_right(a_vec)
        alpha = sum(self.B(c)*b for c, b in zip(v, self.basis_order))
        assert alpha in self.O
        assert self.project(alpha) == mat_alpha
        return alpha
    
    def random_ideal(self):
        "Returns a random ideal of O of norm N"
        a = randint(0, self.N-1)
        b = randint(0, self.N-1)
        assert not a == b == 0

        alpha = self.lift(Matrix(self.Z_N, [[a, b], [0, 0]]))
        assert alpha.reduced_norm() % self.N == 0

        return self.O*self.N + self.O*alpha


def ideal_mod_N(gamma, alpha, N, O):
    "Find quaternion beta_1 such that beta_1*alpha = gamma mod N"
    "Assumes this is a left O0-ideal"
    print(N)
    assert is_prime(N), "Not implemented for non-prime N :(("
    OmodN = O_mod_N(O, N)
    mat_alpha = OmodN.project(alpha)
    mat_gamma = OmodN.project(gamma)
    b = mat_alpha.solve_left(mat_gamma) #Actually this one is always a row vector
    beta_1 = OmodN.lift(b)
    assert beta_1 in O
    return beta_1

def decompose_in_ideal(gamma, alpha, N, O):
    "Write gamma as beta_1*alpha + beta_2*N (assumes this is possible)"
    beta_1 = ideal_mod_N(gamma, alpha, N, O)
    beta_2 = (gamma - beta_1*alpha)/N
    assert beta_1 in O
    assert beta_2 in O
    assert beta_1*alpha + beta_2*N == gamma
    return beta_1, beta_2

def find_correct_pushforward_endo(I, J):
    r"""
    given two left O-ideals I, J, find theta \in O so that J is the pushforward of I by alpha
    """
    O = I.left_order()
    assert O == J.left_order()
    N = ZZ(I.norm())
    assert N == ZZ(J.norm()), "Otherwise this makes no sense"

    OmodN = O_mod_N(O, N)
    alpha_I = ideal_generator(I)
    alpha_J = ideal_generator(J)
    mat_I = OmodN.project(alpha_I)
    mat_J = OmodN.project(alpha_J)

    K_I = mat_I.right_kernel().basis()[0]
    K_J = mat_J.right_kernel().basis()[0]

    system = Matrix(OmodN.Z_N, [[K_I[0], K_I[1], 0, 0], [0, 0, K_I[0], K_I[1]]])
    theta_0 = system.solve_right(K_J) #This one is typically singular, meaning its norm is divisible by N

    theta = Matrix(OmodN.Z_N, [[0,0], [0,0]])
    while theta.is_singular():
        theta = theta_0 + system.right_kernel().random_element()
        theta = Matrix(OmodN.Z_N, [[theta[0], theta[1]], [theta[2], theta[3]]])

    theta = OmodN.lift(theta)
    J_prime = pushforward(I, O*theta)
    alpha = compute_isomorphism(J_prime.left_order(), J.left_order())
    J_prime = alpha**(-1)*J_prime*alpha
    assert J_prime.conjugate().is_equivalent(J.conjugate())
    return theta

def find_correct_pushforward(I, J):
    assert I.norm() == J.norm()
    K = connecting_ideal(I.left_order(), J.left_order())
    while gcd(K.norm(), I.norm()) != 1:
        print("Retrying norm...")
        Kbasis = reduced_basis(K)
        beta = sum(randint(-10, 10)*b for b in Kbasis)
        K = K*(beta.conjugate()/K.norm())

    I_prime = pushforward(I, K)
    alpha = compute_isomorphism(I_prime.left_order(), J.left_order())
    I_prime = alpha**(-1)*I_prime*alpha
    theta = find_correct_pushforward_endo(I_prime, J)
    K = K*(alpha*theta*alpha**(-1))
    return K


###############################################
#                                             #
#    Jumping between quaternion algebras      #
#                                             #
###############################################

def isomorphism_gamma(B_old, B):
    r"""
    Used for computing the isomorphism between B_old and B
    See Lemma 10 [EPSV23]
    """
    if B_old == B:
        return 1
    i_old, j_old, k_old = B_old.gens()
    q_old = -ZZ(i_old**2)
    i, j, k = B.gens()
    q = -ZZ(i**2) 
    p = -ZZ(j**2)
    x, y = DiagonalQuadraticForm(QQ, [1,p]).solve(q_old/q)
    return x + j*y, (x + j_old*y)**(-1)

def eval_isomorphism(alpha, B, gamma):
    r"""
    Given alpha \in B_old, and gamma deteremining the isomorphism from B_old to B,
    returns alpha \in B
    """
    i, j, k = B.gens()
    return sum([coeff*b for coeff, b in zip(alpha.coefficient_tuple(), [1, i*gamma, j, k*gamma])]) 


######## Connecting ########

def find_isomorphism(O1, omega1, omega2):
    #alpha*omega1 - omega2*alpha = 0
    B = O1.quaternion_algebra()
    i,j,k = B.gens()
    a1, b1, c1, d1 = omega1.coefficient_tuple()
    a2, b2, c2, d2 = omega2.coefficient_tuple()

    System = [[a1-a2, i**2*(b1-b2), j**2*(c1 - c2), k**2*(d1 - d2)],
              [(b1*i - i*b2)/i, (i*a1 - a2*i)/i, (j*k*d1 - k*d2*j)/i, (k*c1*j - c2*j*k)/i],
              [(c1*j - c2*j)/j, (i*d1*k - d2*k*i)/j, (c1*j - c2*j)/j, (k*b1*i - b2*i*k)/j],
              [(d1*k - d2*k)/k, (i*c1*j - c2*j*i)/k, (j*b1*i - b2*i*j)/k, (k*a1 - a2*k)/k]]
    alpha = list(Matrix(QQ, System).right_kernel().basis_matrix())[0]
    assert alpha != 0

    return sum(c*g for c,g in zip(alpha, B.basis()))
