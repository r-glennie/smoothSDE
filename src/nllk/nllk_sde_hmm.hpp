
#ifndef _SDE_
#define _SDE_

#include "tr_dens.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla; 
using namespace density; 
using namespace Eigen; 

//' Penalised negative log-likelihood for SDE model
template <class Type>
Type nllk_sde_hmm(objective_function<Type>* obj) { 
    //======//
    // DATA //
    //======//
    DATA_STRING(type); // Model type for each response variable
    DATA_VECTOR(ID); // Time series ID
    DATA_VECTOR(times); // Observation times 
    DATA_VECTOR(coarse_time); // Coarse times 
    DATA_INTEGER(n_coarse); // Number of coarse time points 
    DATA_INTEGER(n_states); // Number of coarse states 
    DATA_MATRIX(obs); // Response variables
    DATA_SPARSE_MATRIX(X_fe); // Design matrix for fixed effects
    DATA_SPARSE_MATRIX(X_re); // Design matrix for random effects
    DATA_SPARSE_MATRIX(S); // Penalty matrix
    DATA_IVECTOR(ncol_re); // Number of columns of S and X_re for each random effect
    DATA_VECTOR(other_data); // Optional extra data needed to evaluate the likelihood
    
    // Number of observations
    int n = obs.rows();
    // Time intervals
    vector<Type> dtimes = diff(times);
    
    //============//
    // PARAMETERS //
    //============//
    PARAMETER_VECTOR(coeff_fe); // Fixed effect parameters
    PARAMETER_VECTOR(log_lambda); // Smoothness parameters
    PARAMETER_VECTOR(coeff_re); // Random effect parameters
    PARAMETER_VECTOR(log_tpm); // Transition probability matrix parameters 
    
    // Derived parameters (linear predictors)
    vector<Type> par_vec = X_fe * coeff_fe + X_re * coeff_re;
    matrix<Type> par_mat(n * n_states, par_vec.size()/(n * n_states));
    for(int i = 0; i < par_mat.cols(); i++) {
        // Matrix with one row for each time step and
        // one column for each parameter
        par_mat.col(i) = par_vec.segment(i * n * n_states, n * n_states);
    }
    
    // Transition probability matrix
    // Transition probability matrix  
    matrix<Type> tpm(n_states, n_states);
    int cur = 0;
    for (int i = 0; i < n_states; ++i) {
        tpm(i, i) = 1;
        for (int j = 0; j < n_states; ++j) {
            if (i != j) {
                tpm(i, j) = exp(log_tpm(cur)); 
                ++cur;
            }
        }
        tpm.row(i) /= tpm.row(i).sum();
    }
    
    // Compute stationary distribution
    matrix<Type> delta(1, n_states);
    matrix<Type> I = matrix<Type>::Identity(n_states, n_states);
    matrix<Type> tpminv = I;
    tpminv -= tpm;
    tpminv = (tpminv.array() + 1).matrix();
    matrix<Type> ivec(1, n_states); for (int i = 0; i < n_states; ++i) ivec(0, i) = 1;
    // if tpm is ill-conditioned then just use uniform initial distribution
    try {
        tpminv = tpminv.inverse();
        delta = ivec * tpminv;
    } catch(...) {
        for (int i = 0; i < n_states; ++i) delta(0, i) = 1.0 / n_states;
    }
    
    //========================//
    // Fine Scale Likelihoods //
    //========================//
    // Initialise log-likelihood
    matrix<Type> fine_llks(n_coarse, n_states); 
    fine_llks.setZero(); 
    // Loop over states
    for (int b = 0; b < n_states; ++b) {
        // Loop over observations
        for(int i = 1; i < n; i ++) {
            // No contribution if first observation of the track
            if(ID(i-1) == ID(i)) {
                fine_llks(coarse_time(i) - 1, b) = fine_llks(coarse_time(i) - 1, b) + tr_dens<Type>(obs.row(i), obs.row(i-1), dtimes(i-1),
                                          par_mat.row(i-1 + b * n), true, type, other_data);
            }
        } 
    }
    
    //=========================//
    // Coarse Scale Likelihood //
    //=========================//
    Type llk = 0; 
    vector<Type> phi(delta); 
    for (int i = 0; i < n_coarse; ++i) {
        // Re-initialise phi at first observation of each time series
        if(i == 0 || ID(i-1) != ID(i)) {
            phi = delta;
        }
        phi = (phi.array() * fine_llks.row(i).array()).matrix();
        phi = phi * tpm;
        sumphi = phi.sum();
        llk = llk + log(sumphi);
        phi = phi / sumphi;
    }
    
    //===================//
    // Smoothing penalty //
    //===================//
    Type nllk = -llk;
    // Are there random effects?
    if(ncol_re(0) > 0) {
        // Index in matrix S
        int S_start = 0;
        
        // Loop over smooths
        for(int i = 0; i < ncol_re.size(); i++) {
            // Size of penalty matrix for this smooth
            int Sn = ncol_re(i);
            
            // Penalty matrix for this smooth
            Eigen::SparseMatrix<Type> this_S = S.block(S_start, S_start, Sn, Sn);
            
            // Coefficients for this smooth
            vector<Type> this_coeff_re = coeff_re.segment(S_start, Sn);
            
            // Add penalty
            nllk = nllk -
                Type(0.5) * Sn * log_lambda(i) +
                Type(0.5) * exp(log_lambda(i)) * 
                density::GMRF(this_S).Quadform(this_coeff_re);
            
            // Increase index
            S_start = S_start + Sn;
        }
    }
    
    return nllk;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif