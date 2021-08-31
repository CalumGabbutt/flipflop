# cython: profile=False

import numpy as np
import scipy.linalg as linalg  

cimport numpy as np
cimport scipy.linalg as linalg  
cimport cython
from libc.math cimport lgamma, log, log1p, pi

cpdef beta_lpdf(double[:] y, double[:] alpha, double[:] beta):
    cdef int N = len(y)
    cdef int Z = len(alpha)
    cdef Py_ssize_t z
    cdef Py_ssize_t n
    cdef double lgamma_alpha
    cdef double lgamma_beta
    cdef double lgamma_alphaplusbeta
    cdef double[:,:] lpk = np.empty((len(alpha), len(y)))

    # Check that alpha and beta are the same shape, if not raise exception
    if len(alpha)!= len(beta):
        raise Exception("alph and beta must be the same shape")

    # Loop through each of the possible Z=2S+1 peaks and calculate the
    # log-pdf of the nth beta value in y[N] for the zth beta distribution
    with cython.boundscheck(False):
        with cython.wraparound(False):
            with cython.nonecheck(False):
                for z in range(Z):

                    # Precompute the log-gamma constants
                    lgamma_alpha = lgamma(alpha[z])
                    lgamma_beta = lgamma(beta[z])
                    lgamma_alphaplusbeta = lgamma(alpha[z] + beta[z])

                    for n in range(N):
                        # Calculate the log pdf of the beta distribution for the nth 
                        # datapoint drawn from the zth peak
                        lpk[z, n] = ((alpha[z] - 1)*log(y[n]) + (beta[z] - 1) * log1p(-y[n]) - 
                                    lgamma_alpha - lgamma_beta + lgamma_alphaplusbeta)


    return lpk

cdef generateRateMatrix(int S, double lam, double mu, double gamma, int[:,:] stateVar, int num_states):
    cdef np.ndarray[np.float_t, ndim=2] RateMatrix 
    cdef Py_ssize_t down
    cdef Py_ssize_t across
    cdef int k_down
    cdef int k_across
    cdef int k
    cdef int m

    # Preallocate a matrix containing all zeros
    RateMatrix = np.zeros((stateVar.shape[0], stateVar.shape[0]))

    # Loop through the state variables (k - the number of cells with one 
    # methylated allele, m - the number of cells with both alleles methylated)
    # to construct the transition matrix, which defines the rate at which the 
    # system transitions from state (k, m) to state (k', m').
    # For a derivation of these equations, please see Gabbutt et al. (2020)
    with cython.boundscheck(False):
        with cython.wraparound(False):
            with cython.nonecheck(False):
                for down in range(num_states):
                    k_down = stateVar[down, 0]
                    m_down = stateVar[down, 1]
                    for across in range(num_states):
                        k = stateVar[across, 0]
                        m = stateVar[across, 1]

                        if k == k_down-1 and m == m_down :
                            RateMatrix[down, across] = (S-m-k) * (k*lam/(S-1) + 2*mu)
                        elif k == k_down and m == m_down-1:
                            RateMatrix[down, across] = m * (S-m-k) * lam / (S-1)
                        elif k == k_down+1 and m == m_down-1:
                            RateMatrix[down, across] = k * (m*lam/(S-1)+mu)
                        elif k == k_down+1 and m == m_down:
                            RateMatrix[down, across] = k * ((S-m-k)*lam/(S-1)+gamma)
                        elif k == k_down and m == m_down+1:
                            RateMatrix[down, across] = m * (S-m-k) * lam / (S-1)
                        elif k == k_down-1 and m == m_down+1:
                            RateMatrix[down, across] = m * (k*lam/(S-1)+2*gamma)
                        elif k == k_down and m == m_down:
                            RateMatrix[down, across] = -(2*((k+m)*(S-m-k)+k*m)*lam/(S-1) + (k+2*m)*gamma + (2*S-(k+2*m))*mu) 

    return RateMatrix

cdef generateStateVar(int S, int num_states):
    cdef int[:,:] stateVar = np.empty((num_states, 2), dtype=np.int32)
    cdef int i
    cdef int k
    cdef int m

    # Loop through each combination of k and m and ensure that only valid states
    # (that satifsy 0<= k+m <= S) are included. There should be 1/2 * (S+1) * (S+2)
    # states.
    i = 0
    with cython.boundscheck(False):
        with cython.wraparound(False):
            with cython.nonecheck(False):
                for m in range(S+1):
                    for k in range(S+1):
                        if k+m <=S :
                            stateVar[i, 0] = k
                            stateVar[i, 1] = m

                            i += 1

    return stateVar


cdef findProbDist(np.ndarray[np.float_t, ndim=2] RateMatrix, double[:] InitialConditions, int[:,:] stateVar, double age, int S, int num_states):
    cdef double[:] ProbStates
    cdef double[:] ProbDist
    cdef Py_ssize_t index
    cdef int k
    cdef int m

    # Having constructed the transition matrix, use matrix exponentiation 
    # to solve the system of linear ordinary differential equations
    ProbStates = linalg.expm(RateMatrix * age) @ InitialConditions 

    # The ProbStates contains the probability that the system contains
    # k singly methylated cells and m doubly methylated cells, however 
    # we actually measure the population level methylation, hence we 
    # marginalise over the combination of k and m states that give rise 
    # to a given z state (e.g. P(z=0) = P(k=0, m=0), 
    # P(z=2) = P(k=2, m=0) + P(k=0, m=1) etc.)
    ProbDist = np.zeros(2*S+1)
    with cython.boundscheck(False):
        with cython.wraparound(False):
            with cython.nonecheck(False):
                for index in range(num_states):
                    k = stateVar[index, 0]
                    m = stateVar[index, 1]
                    ProbDist[k+2*m] += ProbStates[index]

    return ProbDist

cpdef runModel(int S, double lam, double mu, double gamma, double age):
    cdef int num_states = int(0.5 * ((S+1) * (S+2)))
    cdef int[:,:] stateVar 
    cdef np.ndarray[np.float_t, ndim=2] RateMatrix 
    cdef double[:] InitialConditions
    cdef double[:] ProbDist

    # Generate the possible (k, m) states for a given S value
    stateVar = generateStateVar(S, num_states)

    # For a given set of stem cell dynamics parameters, generate the transition 
    # matrix 
    RateMatrix = generateRateMatrix(S, lam, mu, gamma, stateVar, num_states)

    # Define the initial conditions. Given no other information, we assume that 
    # the population is initially clonal with an equal probability of all the alleles
    # being methylated or unmethylated 
    InitialConditions = np.zeros(num_states)
    InitialConditions[0] = 0.5 
    InitialConditions[-1] = 0.5

    # Given the above initial conditions, calculate the probability distribution at 
    # time age
    ProbDist = findProbDist(RateMatrix, InitialConditions, stateVar, age, S, num_states)

    return ProbDist