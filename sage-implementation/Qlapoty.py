from sage.all import *
from helpers.quaternion_helpers import *
from helpers.sum_of_squares import sum_of_squares
from time import time



#--------------------------------------------------------------------------------------------------
#----------------------------------------- Qlapoti-------------------------------------------------
#--------------------------------------------------------------------------------------------------

def smallish_gen_old(I, N, I_basis=None):
    if not I_basis:
        print("SHOULDNT GET CALLED")
        I_basis = reduced_basis(I)
    #N = I.norm() This is slow, pass instead

    while True:
        alpha = sum(randint(1,100000)*gen for gen in I_basis)

        a_alpha = alpha.coefficient_tuple()[0]
        if gcd(2*a_alpha, N) == 1 and gcd(alpha.reduced_norm(), N*N) == N:
            break

    #assert gcd(alpha.reduced_norm(), N**2) == N
    assert alpha in I
    return alpha



def _succ_min(L, p):
    fourth_root_p = round(p**(1/4), 10)
    lam1 = round(L.row(0).norm()/fourth_root_p, 10)
    lam2 = round(L.row(1).norm()/fourth_root_p, 10)
        
    return lam1, lam2

def _succ_min_v2(L, N):
    sqrtN = round(N.sqrt(), 10)
    lam1 = round(L.row(0).norm()/sqrtN, 10)
    lam2 = round(L.row(1).norm()/sqrtN, 10)
        
    return (lam1 + lam2)/2


def predict_endtype(tlist):
    #Pain and Misery in concentrated form
    #transformlist = {1: [[2, 0, 1, 2], [2, 2, 3, 0], [2, 2, 1, 0], [2, 0, 3, 2]], 2: [[0, 2, 2, 1], [0, 2, 2, 3], [2, 2, 0, 1], [2, 2, 0, 3]]}
    # t1, t2, t3, t4 = tlist
    # Again, disagreeing on what j and k is
    t1, t2, t4, t3 = tlist

    if t1 == 2:
        if t2 == 2:
            if t3 == 1 and t4 == 0:
                return 1
            elif t3 == 3 and t4 == 2:
                return 1
            elif t3 == 0:
                if t4 == 1:
                    return 2
                elif t4 == 3:
                    return 2
                else:
                    return 0
            else:
                return 0
        elif t2 == 0:
            if t3 == 1 and t4 == 2:
                return 1
            elif t3 == 3 and t4 == 2:
                return 1
            else:
                return 0
        else:
            return 0
    elif t1 == 0:
        if t2 == 2 and t3 == 2:
            if t4 == 1:
                return 2
            elif t4 == 3:
                return 2
            else:
                return 0
        else:
            return 0
    else:
        return 0

def qlapoti(J, e, allow_diags = True, odd_norm_output = False, stats = False):
    I, betaij = reduced_ideal(J, return_elt=True) #The smallest ideal in the class
    p = I.quaternion_algebra().ramified_primes()[0]
    num_tries = 0
    num_cornacchia = 0

    # The algorithm...
    N = ZZ(I.norm())

    solved = False
    I_basis = reduced_basis(I)

    keep_alpha = False
    lam = 1


    while not solved:
        num_tries += 1
        if num_tries % 20000 == 0:
            lam = 1
            keep_alpha = False
        #print("\n\nTrying new gen....")
        if keep_alpha:
            alpha += alpha_0
            lam += 1
            a_alpha = list(alpha.coefficient_tuple())[0]
            #multiply alpha by an invertible lambda, and check that 2*alpha is still invertible
            while gcd(lam, N) > 1:
                alpha += alpha_0
                lam += 1
                a_alpha = list(alpha.coefficient_tuple())[0]
        else:
            alpha_0 = smallish_gen_old(I, N, I_basis=I_basis)
            alpha = alpha_0
            #print(f"alpha = {alpha}")
            alpha_0_norm = alpha.reduced_norm()
            #print(f"norm = {factor(ZZ(alpha_0_norm), limit=100)}")

        a_alpha = list(alpha.coefficient_tuple())[0]
        b_alpha = list(alpha.coefficient_tuple())[1]

        M = ZZ(2**e - 2*lam**2*alpha_0_norm/N)
        if M < 0:
            print("M already too small! NB! Shouldnt happen")
            continue
        #N(a1**2 + a2**2 + b1**2 + b2**2) + 2*a_alpha1*a1 + 2*a_alpha2*a2 + 2*balpha1*b1 + 2*balpha2*b2 = M
        #b1 + balpha1^-1*balpha2*b2) = M*(2*balpha1)^-1 (mod N)
        #A + (2*b_alpha/2*a_alpha)B = M/(2*a_alpha) (mod N)

        if keep_alpha:
            a_alpha_inv2 = inverse_mod(ZZ(2*a_alpha), N)
            T = (M*a_alpha_inv2) % N
            v_target = vector(ZZ, [-T, 0])
        else:
            a_alpha_inv2 = inverse_mod(ZZ(2*a_alpha), N)
            x = ZZ((2*b_alpha*a_alpha_inv2) % N)
            L = Matrix(ZZ, [[N-x, 1], [N, 0]])
            L = L.LLL() # Make custom gauss reduction if necessary for speedup...

            T = ZZ((M*a_alpha_inv2) % N)
            v_target = vector(ZZ, [-T, 0])
            L_inv = Matrix(QQ, L).inverse()
            if _succ_min(L, p)[1] < 1 and gcd(2*a_alpha, 2*b_alpha) == 1: #Just make sure our lattice isnt completely unreasonable
                if not stats:
                    keep_alpha = True

        AB_vec = vector([round(c) for c in v_target*L_inv])
        v_close = AB_vec*L
    
        A, B = v_close - v_target
        
        assert (2*(a_alpha*A + b_alpha*B)) % N == M % N

        M2 = M - 2*a_alpha*A - 2*b_alpha*B
        M2 = ZZ(M2/N)

        # Complete the square
        # (a2 = A - a1, b2 = B - b1)
        #2*a1^2 + A^2 - 2Aa1 + 2*b1^2 + B^2 - 2Bb1 = M2
        #M3 = M2 - A**2 - B**2

        #(2*a1)^2 + (2*b1)^2 - 4*Aa1 - 4*Bb1 = (2*a1 - A)^2 + (2*b1 - B)^2 = 4*M3 + A^2 + B^2
        # 4*M3 = (2*a1)^2 + (2*b1)^2 - 4*A*a1 - 4*B*b1 = (2*a1 - A)^2 + (2*b1 - B)^2 - (A^2 + B^2)

        M4 = 2*M2 - A**2 - B**2
        if not allow_diags and M4 < 0:
            continue

        #Unsolvable cases
        if M4 % 8 == 0:
            continue
        if A % 2 == B % 2 == 0:
            if M4 % 4 != 0:
                continue
        elif A % 2 == B % 2 == 1:
            if M4 % 4 != 2:
                continue
        else:
            if M4 % 4 != 1:
                continue

        #### Can we predict already here?
        
        num_cornacchia += 1

        ab = sum_of_squares(M4)

        if not ab:
            continue

        solved = True
        ad1, bd1 = ab

        if ad1 % 2 != A % 2:
            temp = ad1
            ad1 = bd1
            bd1 = temp

        assert ad1 % 2 == A % 2 and bd1 % 2 == B % 2
        solved = True
        a1 = ZZ((ad1 + A)/2)
        b1 = ZZ((bd1 + B)/2)
        a2 = A - a1
        b2 = B - b1

        i, j, k = I.quaternion_algebra().gens()
        gamma1 = a1 + i*b1
        gamma2 = a2 + i*b2

        mu1 = N*gamma1 + alpha
        mu2 = N*gamma2 + alpha

        θ = mu2 * mu1.conjugate() / N

        if allow_diags:
            d1 = mu1.reduced_norm()/N
            d2 = mu2.reduced_norm()/N
            if not odd_norm_output or d1 % 2 == d2 % 2 == 1:
                break
            else:
                solved = False
                keep_alpha = False #Some congruence condition it seems like, so can get stuck
                lam = 1
                continue
        
        if θ.denominator() == 1:
            # If this should be possible, it must come from a "transformation"
            theta_profile = [c % 4 for c in θ]
            res = predict_endtype(theta_profile)
            if res == 0:
                solved = False
                keep_alpha = False #Some congruence condition it seems like, so can get stuck
                lam = 1
                continue
            elif res == 1:
                temp = (mu1 + mu2)/2
                mu2 = (mu1 - mu2)/2
                mu1 = temp
            else:
                temp = (mu1 + i*mu2)/2
                mu2 = (mu1 - i*mu2)/2
                mu1 = temp
        
            e = e-1 
            θ = mu2 * mu1.conjugate() / N
        else:
            # Check odd norm
            theta_profile = [c % 4 for c in 2*θ]
            tx, ty, tz, tw = 2*θ 
            if not ((tx % 2 != tz % 2) and ([c % 4 for c in 2*θ].count(2) == 1)):
                solved = False
                keep_alpha = False #Some congruence condition it seems like, so can get stuck
                lam = 1
                continue
                
        d1 = mu1.reduced_norm()/N
        d2 = mu2.reduced_norm()/N

    assert allow_diags or θ.denominator() == 2 
    
    assert solved
    
    assert N*(a1**2 + a2**2 + b1**2 + b2**2) + 2*a_alpha*(a1 + a2) + 2*b_alpha*(b1 + b2) == M

    assert mu1 in I
    assert mu2 in I
    assert mu1.reduced_norm() + mu2.reduced_norm() == 2**e*N
    assert θ == mu2 * mu1.conjugate() / N

    #J1 = I*(mu1.conjugate()/N)
    #J2 = I*(mu2.conjugate()/N)
    #assert J1.norm() + J2.norm() == 2**e

    if stats:
        return mu1, mu2, d1, d2, θ, N, I, betaij, num_tries, num_cornacchia
    return mu1, mu2, d1, d2, θ, N, I, betaij, e 


#--------------------------------------------------------------------------------------------------
#----------------------------------------- Qlapoty V1----------------------------------------------
#--------------------------------------------------------------------------------------------------


#Function called to sample a generator in the first version of the algorithm. 
def Sample_alpha0_regular(HNF_x,HNF_y, N, Bound): 
    Sampling_bound = Bound >> 1

    while True: 
        c3 = randint(-Sampling_bound, Sampling_bound)
        c4 = randint(-Sampling_bound, Sampling_bound)
        c3 = ((c3 >> 1) << 1)     #Ensure it is even
        c4 = ((c4 >> 1) << 1) + 1 #Ensure it is odd


        #reduce mod N the first 2 coefficient. 
        z_I = (HNF_x + quat_i*HNF_y)
        z_c = (c3 + quat_i*c4)
        z_cI = z_c*z_I

        #Apply reduction mod N
        a_alpha0 = ZZ((z_cI[0] % (2*N)))/2  
        b_alpha0 = ZZ((z_cI[1] % (2*N)))/2

        if gcd(2*a_alpha0, 2*N) != 1:
            continue

        assert(2*a_alpha0 % 2 == 1)
        assert(2*b_alpha0 % 2 == 0)

        alpha_0 = (a_alpha0 + quat_i*b_alpha0) + (z_c*quat_j)/2
        
        #Compute the lattice
        a_alpha_inv2 = inverse_mod(ZZ(2*a_alpha0), (4*N))
        x = ZZ((2*b_alpha0*a_alpha_inv2) % (4*N))
        L = Matrix(ZZ, [[-x, 1], [4*N, 0]])
        L = L.LLL()
        
        if _succ_min_v2(L, N) <= 3: #Just make sure our lattice isnt completely unreasonable
            good_generator = True 
            L_inv = Matrix(QQ, L).inverse()
            return alpha_0, L,  L_inv
        
        
        

        
def qlapoty(J, e, stats = False):
    I, betaij = reduced_ideal(J, return_elt=True) #The smallest ideal in the class
    p = I.quaternion_algebra().ramified_primes()[0]
    num_tries_alpha = 0
    num_tries_st = 0
    num_cornacchia = 0
    num_parity_st = 0

    # The algorithm...
    N = ZZ(I.norm())
    
    good_generator = False
    solved = False
    HNF_x,HNF_y  = HNF_coeffs(I)
    
    #Fist loop
    alpha_0, L,  L_inv = Sample_alpha0_regular(HNF_x,HNF_y, N, 128)
    
    alpha_0_norm = alpha_0.reduced_norm()
    r0 = alpha_0_norm/N
    bitrev_r0_mod2 = 1 - (r0 % 2)
    
    alpha = -alpha_0
    lam = -1

    while not solved: 
        alpha += 2*alpha_0
        lam += 2
        if gcd(lam, N) != 1:
            continue
        num_tries_st += 1

        a_alpha = list(alpha.coefficient_tuple())[0]
        b_alpha = list(alpha.coefficient_tuple())[1]

        r = lam**2*alpha_0_norm/N
        M = ZZ(2**e - 2*r - N - 2*N*bitrev_r0_mod2)
        if M < 0:
            print("M already too small! NB! Shouldnt happen")
            continue
        
        a_alpha_inv2 = inverse_mod(ZZ(2*a_alpha), 4*N)
        T = (M*a_alpha_inv2) % (4*N)
        v_target = vector(ZZ, [-T, 0])

        AB_vec = vector([round(c) for c in v_target*L_inv])
        v_close = AB_vec*L
    
        A, B = v_close - v_target
        
        assert (2*(a_alpha*A + b_alpha*B)) % (4*N) == M % (4*N)

        M2 = M - 2*a_alpha*A - 2*b_alpha*B
        M2 = ZZ(M2/N) + 1 + 2*bitrev_r0_mod2
        
        M4 = 2*M2 - A**2 - B**2
        
        #Diagonal isogeny 
        if  M4 % 4 != 1:
            num_parity_st += 1
            continue
                 
        #Early reject: 
        if M4 < 0:
            continue
        
                      
        #Ensure odd solutions
        assert (r + ((M2 % 4) >> 1)) % 2 == 1


        num_cornacchia += 1
        ab = sum_of_squares(M4)

        if not ab:
            continue
        else: 
            solved = True

    ad1, bd1 = ab
    if ad1 % 2 != A % 2:
        temp = ad1
        ad1 = bd1
        bd1 = temp


    assert ad1 % 2 == A % 2 and bd1 % 2 == B % 2
    solved = True
    a1 = ZZ((ad1 + A)/2)
    b1 = ZZ((bd1 + B)/2)
    a2 = A - a1
    b2 = B - b1
    
    
    i, j, k = I.quaternion_algebra().gens()
    gamma1 = a1 + i*b1
    gamma2 = a2 + i*b2

    mu1 = N*gamma1 + alpha
    mu2 = N*gamma2 + alpha

    θ = mu2 * mu1.conjugate() / N
        
    d1 = mu1.reduced_norm()/N
    d2 = mu2.reduced_norm()/N
                
    assert solved    
    assert N*(a1**2 + a2**2 + b1**2 + b2**2) + 2*a_alpha*(a1 + a2) + 2*b_alpha*(b1 + b2) == M + N + 2*N*bitrev_r0_mod2

    assert d1 % 2 == 1
    assert d2 % 2 == 1
    assert mu1 in I
    assert mu2 in I
    assert mu1.reduced_norm() + mu2.reduced_norm() == 2**e*N 
    assert θ == mu2 * mu1.conjugate() / N
    assert 2*(θ[0] + θ[1]) % 2 == 1
        
    if stats:
        return mu1, mu2, d1, d2, θ, N, I, betaij, num_tries_alpha, num_tries_st, num_cornacchia, num_parity_st
    return mu1, mu2, d1, d2, θ, N, I, betaij





#--------------------------------------------------------------------------------------------------
#------------------------------------------ V2 ----------------------------------------------------
#--------------------------------------------------------------------------------------------------



    
def basis_to_good_HNF(mi):
    #Reorder the HNF to be in the right shape
    mt = mi.transpose()
    mt = matrix(QQ,[mt[3],mt[2],mt[1],mt[0]])
    mmixed = mt.transpose()
    mmixed = mmixed*2
    mmixed = mmixed.change_ring(ZZ)
    me = mmixed.echelon_form()
    me = matrix(QQ,[me[3],me[2],me[1],me[0]])
    mt = me.transpose()
    mt = matrix(QQ,[mt[3],mt[2],mt[1],mt[0]])
    reorder = mt.transpose()
    return reorder


def HNF_coeffs(I):
    """Given a left O0 ideal, computes its x, y HNF coefficient

    Args:
        I (O0-ideal): an ideal 

    Returns:
        QQ: x, y coeffient
    """
    M = basis_to_good_HNF(I.basis_matrix())
    x =  M[2][0] 
    y =  M[2][1] 

    return x, y



def Gaussian_int_GCD(a, b): 
    # Takes as input a, b in Z[i] over the Quaternion algebra and 
    # returns gcd_{Z[i]}(a, b) over the quaternion algebra.         
    a_2 = ZI([a[0], a[1]])
    b_2 = ZI([b[0], b[1]])
    c_2 = gcd(a_2, b_2)
    return c_2[0] + quat_i*c_2[1]
    


        
def find_generator_minima(omega, alpha, N):
    #find x,y such that omega = x*b1 + y*b2 for lattice of det 4N
    a_alpha = 2*alpha[0]
    b_alpha = 2*alpha[1]
    
    div_a = gcd(a_alpha, 4*N)
    N_mod = ZZ(4*N / div_a)
    a_inv = inverse_mod(ZZ(a_alpha/ div_a), N_mod)
    const = ZZ((b_alpha*a_inv) % N_mod)

    # x = (d*omega_0 - b*omega_1)/det(B) = omega_1*(N_mod/4N)
    x = omega[1]/div_a
    #y = (a*omega_1 - c*omega_0)/det(B) = (const*omega_1 + div_x*omega_0)/(4N)
    y = (const*omega[1] + div_a*omega[0])/(4*N)
    
    return x, y, const, div_a

def Final_lattice_reduction(L):
    #Take a lattice whose first number is already a minima

    v1 =  L[0][0]**2 + L[0][1]**2
    v2 =  L[1][0]**2 + L[1][1]**2
    m = L[0][0]*L[1][0] + L[0][1]*L[1][1]
    if v1 < v2:
        #first minima
        q = round(m/v1)
        L_ret = Matrix(ZZ, [[L[0][0], L[0][1]], [L[1][0] - q*L[0][0], L[1][1] - q*L[0][1]]])
        
    else:
        #second minima
        q = round(m/v2)
        L_ret = Matrix(ZZ, [[ L[0][0] - q*L[1][0],L[0][1] - q*L[1][1]], [L[1][0], L[1][1]]])
        
    return L_ret



def sample_generator(bound_a, bound_b, omega_I, N, HNF_x, HNF_y, z_c0):
    """Implementation of the first loop of SQIsign NIST 1. 

    Args:
        bound_a (int): lower bound of sampling
        bound_b (int): upper bound of sampling
        omega_I (Z[i]): minima of z_c0*(HNF_x+iHNF_y)
        N (int): Norm of the ideal
        HNF_x (int): First component of the HNF
        HNF_y (int): Second component of the HNF
        z_c0 (Z[i]): A ideal of norm 5^f that complete z_I. 

    Returns:
        _type_: _description_
    """
    Big_enough = (bound_b - bound_a) >= (log(N,2)**2)
    ell = gcd(omega_I[0], omega_I[1])
    
    while True:
        
        #Random sampling
        
        c3 = randint(ceil(bound_a/sqrt(2)), bound_b)
        c4 = randint(floor(sqrt(max(0, bound_a**2 - c3**2))), floor(sqrt(bound_b**2 - c3**2)))
                
        if ((c3 + c4) % 2) == 0:
            continue
        
        
        if gcd(c3, c4) != 1:
            continue

        #c3, c4 = first_const, sec_const
        
        #Condition to succeed if b - a is big 
        if Big_enough:
            if gcd(c3**2+c4**2, N) != 1:
                continue
        else:
            if gcd(c3**2+c4**2, ell) != 1:
                continue

        
        #flip sign
        if randint(0,1):
            c3 = -c3
        if randint(0,1):
            c4 = -c4
        
        
            
        z_c = (c3 + quat_i*c4)
        #compute GCD over Z[i]
        gamma = Gaussian_int_GCD(z_c.conjugate(), omega_I) 
        Gaussian_coprime =  ((gamma[0]==0 and (gamma[1]==1 or gamma[1]==-1)) or (gamma[1]==0 and (gamma[0]==1 or gamma[0]==-1)))
        while  not Gaussian_coprime: 
            z_c = (z_c * gamma**2) / (gamma[0]**2 + gamma[1]**2)
            gamma = Gaussian_int_GCD(z_c.conjugate(), omega_I) 
            Gaussian_coprime =  ((gamma[0]==0 and (gamma[1]==1 or gamma[1]==-1)) or (gamma[1]==0 and (gamma[0]==1 or gamma[0]==-1)))

        #compute alpha
        alpha_0 = z_c * z_c0 * (HNF_x + quat_i * HNF_y + quat_j)/2
        
        #Reduce mod N the first 2 coeff
        a_alpha0 = ZZ(2*alpha_0[0] % (4*N))/2
        b_alpha0 = ZZ(2*alpha_0[1] % (4*N))/2
        alpha_0 = (a_alpha0 + quat_i*b_alpha0) + (z_c * z_c0 * quat_j)/2

        # assert geometric conditions
        assert((2*(a_alpha0 + b_alpha0 )) % 2 == 1)


        #compute stuff for target rank 2 lattice
        g, u, v = xgcd(2*a_alpha0, 2*b_alpha0)
        vv, ss, _ = xgcd(g, 4*N)
        
        # assert we are working with a generator
        assert(vv == 1)
        
        Tar1, Tar2 = (u*ss), (ss*v)
        assert((2*a_alpha0*Tar1 + 2*b_alpha0*Tar2) % (4*N) == 1)
        
        
        
        # Compute small minima method: 
        minima_alpha_0 = z_c * omega_I

        # check that our special minima is indeed in the lattice. 
        assert(2*(minima_alpha_0.conjugate()*alpha_0)[0] % (4*N) == 0)
        
        #Given canonical basis of L, compute a1, a2 such to express our minima in it.  
        a1, a2, target, div_a = find_generator_minima(minima_alpha_0, alpha_0, N)       
        vv, a3, a4 = xgcd(a1, a2)   #a1*a3 + a2*a4 = 1
        assert(vv == 1)
        a3, a4 = -a4, a3            #a1*a4 - a2*a3 = 1
        
                                
        L = Matrix(ZZ, [[-target, div_a], [(4*N)/div_a, 0]])
        L_transf = Matrix(ZZ, [[a1, a2], [a3, a4]])
        L = L_transf*L
        
        #Fast LLL using the fact that minima_alpha_0 is a minima.  
        L = Final_lattice_reduction(L)
        
        #As we can ensure omega^I is a minima, then can drop the LLL      
        #L, Transf = L.LLL(transformation=True)
        #print(Transf)

        L_inv = Matrix(QQ, L).inverse()

        return alpha_0, L, L_inv, Tar1, Tar2

    
    

def new_qlapoty(J, e, stats = False):
    
    I, betaij = reduced_ideal(J, return_elt=True) #The smallest ideal in the class
    num_tries_alpha = 0
    num_tries_st = 0
    num_cornacchia = 0


    # The algorithm...
    N = ZZ(I.norm())
    
    good_generator = False
    solved = False

    #Compute the HNF lattice
    HNF_x,HNF_y  = HNF_coeffs(I)
    
    
    #Compute the minima of lattice given by the HNF.
    div_x = gcd(HNF_x, 4*N)
    N_mod = ZZ(4*N / div_x)
    x_inv = inverse_mod(ZZ(HNF_x/ div_x), N_mod)
    const = ZZ((HNF_y*x_inv) % N_mod)
    Ideal_lattice = Matrix(ZZ, [[-const, div_x], [N_mod, 0]])

    Ideal_lattice_reduced = Ideal_lattice.LLL() # reduce the lattice
    s_0 = Ideal_lattice_reduced.row(0)[0]       #First minima
    t_0 = Ideal_lattice_reduced.row(0)[1]
    
    omega_I = (s_0 + quat_i*t_0)    

    first_minima_0 = sqrt(s_0**2 + t_0**2)
    Bound_alpha = floor(sqrt(float(N*(2**(e+3) -16*N))/p))
        
    # Setup Bound 
    const_1 = 0.7 #Ensure that b - a >= 2
    const_2 = 1
    Bound_1 = (const_1*sqrt(N))
    Bound_2 = (const_2*sqrt(8*N))


    #Following the theorem, if n_I, if we cannot find a solution and must therefore reject.   
    if (first_minima_0 < sqrt(2)*Bound_1 / Bound_alpha):
        print("eject. Minima too small")
        return 

    #Set our sampling bounds. 
    bound_a = (Bound_1/ first_minima_0)
    bound_b = min((Bound_alpha/sqrt(2)),(Bound_2/first_minima_0))



    #In case our sampling bound are too small, 
    if (bound_b - bound_a <= ceil(log(N,2)**2)):
        ell = gcd(omega_I[0], omega_I[1])
        g = ceil(log(log(N,2)**4, 5))
        pow_5 = gcd(5**g, ell)
        f = log(pow_5, 5)
        omega_I = omega_I/pow_5
                
        # Test that (2 - i) does not divide omega_I
        if ((2*omega_I[0] - omega_I[1]) % 5 != 0):
            z_c0 = (2 + quat_i)**f
            omega_I = z_c0 * omega_I
        else:   #it does divide (2 - i). 
            z_c0 = (2 - quat_i)**f
            omega_I = z_c0 * omega_I

        #update bounds
        bound_a = ceil(bound_a*sqrt(pow_5))
        bound_b = floor(bound_b*sqrt(pow_5))
        
    
    #Sample alpha_0
    alpha_0, L, L_inv, Tar1, Tar2 = sample_generator(bound_a, bound_b, omega_I, N, HNF_x, HNF_y, z_c0)

    
    
    alpha_0_norm = alpha_0.reduced_norm()        
    r_0 = alpha_0_norm / N
    delta_odd_sol = 1 + 2*(1 - (r_0 % 2))
    
    alpha = -alpha_0
    lam = -1
    
    bound_lambda = floor(sqrt(2**(e-2) / r_0))
    
    while lam < bound_lambda: 
        num_tries_st += 1
        alpha += 2*alpha_0
        lam += 2
        if gcd(lam, N) != 1:
            continue 

        a_alpha = list(alpha.coefficient_tuple())[0]
        b_alpha = list(alpha.coefficient_tuple())[1]

        r = lam**2*r_0
        M = ZZ(2**e - 2*r - N*delta_odd_sol)
        if M < 0:
            print("M already too small! NB! Shouldnt happen")
            continue
        

        #compute target
        lamb_inv = inverse_mod(ZZ(lam), 4*N)
        MS = (lamb_inv * M) % (4*N) # M * lamb^{-1}
        V1, V2 = (- MS * Tar1) % (4*N), (- MS * Tar2) % (4*N)
        
        assert((2*a_alpha*V1 + 2*b_alpha*V2 )%(4*N) == -M %(4*N))        
        v_target = vector(ZZ, [V1, V2])



        AB_vec = vector([round(c) for c in v_target*L_inv])
        v_close = AB_vec*L
    
        # solution of modular equation
        A, B = v_close - v_target
        
        assert (2*(a_alpha*A + b_alpha*B) % (4*N) == M % (4*N))

        M2 = M - 2*a_alpha*A - 2*b_alpha*B + delta_odd_sol*N
        M2 = ZZ(M2/N) 
        
        M4 = 2*M2 - A**2 - B**2
                                 
        #Early reject: 
        if M4 < 0:
            continue
        
        #Diagonal isogeny 
        if  M4 % 4 != 1:
            continue

                      
        #Ensure odd solutions
        assert (r + ((M2 % 4) >> 1)) % 2 == 1


        num_cornacchia += 1
        ab = sum_of_squares(M4)

        if not ab:
            continue
        else: 
            solved = True
            break

    #TODO: implement a loop back to first line
    if (not solved):
        print("Error: Did not find valid s, t ")

    ad1, bd1 = ab
    if ad1 % 2 != A % 2:
        temp = ad1
        ad1 = bd1
        bd1 = temp


    assert ad1 % 2 == A % 2 and bd1 % 2 == B % 2
    solved = True
    a1 = ZZ((ad1 + A)/2)
    b1 = ZZ((bd1 + B)/2)
    a2 = A - a1
    b2 = B - b1
    
    
    i, j, k = I.quaternion_algebra().gens()
    gamma1 = a1 + i*b1
    gamma2 = a2 + i*b2

    mu1 = N*gamma1 + alpha
    mu2 = N*gamma2 + alpha

    θ = mu2 * mu1.conjugate() / N
        
    d1 = mu1.reduced_norm()/N
    d2 = mu2.reduced_norm()/N
                
    assert solved    
    assert N*(a1**2 + a2**2 + b1**2 + b2**2) + 2*a_alpha*(a1 + a2) + 2*b_alpha*(b1 + b2) == M + N*delta_odd_sol

    assert d1 % 2 == 1
    assert d2 % 2 == 1
    assert mu1 in I
    assert mu2 in I
    assert mu1.reduced_norm() + mu2.reduced_norm() == 2**e*N 
    assert θ == mu2 * mu1.conjugate() / N
    assert 2*(θ[0] + θ[1]) % 2 == 1
    
        
    return mu1, mu2, d1, d2, θ, N, I, betaij




if __name__ == "__main__":
    proof.all(False)
    param = "NIST-I"

    if param == "NIST-I":
        p = 2**248*5 - 1
        e = 246
        f = 5
    elif param == "NIST-III":
        p = 2**376*65 - 1
        e = 374
        f = 65 
    elif param == "NIST-V":
        p = 2**500*27 - 1
        e = 498
        f = 27
    assert is_prime(p)
    assert p%4 == 3
    
    B = QuaternionAlgebra(-1, -p)
    #O0 = B.maximal_order()
    # Sage's default is not the isogeny gringos default O0
    quat_i, quat_j, quat_k = B.gens()
    O0 = B.quaternion_order([1, quat_i, (quat_i + quat_j) / 2, (1 + quat_k) / 2])
    ZI = QuadraticField(-1).ring_of_integers()
    
    
    
    
    N = next_prime(randint(p, p**2))

    O_N = O_mod_N(O0, N)

    timings = []
    n_runs = 250
    
    
#--------------------------------------------------------------------------------------------------
#------------------------------------------ MAIN --------------------------------------------------
#--------------------------------------------------------------------------------------------------

    from tqdm import tqdm
    two_torsion = e 

    
    for _ in tqdm(range(n_runs)):        
        I = O_N.random_ideal()
        t_start = time()
        qlapoti(I, e, allow_diags= False , odd_norm_output=True)
        timings.append(time()-t_start)
    print(f"All runs done! Qlapoti without diagonal isogeny took on average {round(sum(timings)/n_runs, 5)}")
    timings = []
    for _ in tqdm(range(n_runs)):        
        I = O_N.random_ideal()
        t_start = time()
        qlapoti(I, e, allow_diags= True , odd_norm_output=True)
        timings.append(time()-t_start)
    print(f"All runs done! Qlapoti with diagonal isogeny took on average {round(sum(timings)/n_runs, 5)}")
    timings = []
    for _ in tqdm(range(n_runs)):        
        I = O_N.random_ideal()
        t_start = time()
        qlapoty(I, two_torsion , stats=False)
        timings.append(time()-t_start)
    print(f"All runs done! Qlapoty Loop 2 took on average {round(sum(timings)/n_runs, 5)}")
    
    timings = []
    for _ in tqdm(range(n_runs)):        
        I = O_N.random_ideal()
        t_start = time()
        new_qlapoty(I, two_torsion)
        timings.append(time()-t_start)
    print(f"All runs done! full Qlapoty algorithm took on average {round(sum(timings)/n_runs, 5)}")
