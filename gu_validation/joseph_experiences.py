#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:04:54 2020

@author: Joseph Muré
"""

from Classe_Gu import *

input_sample = ot.Sample(np.array([[0.1,0.6],[0.3,0.4],[0.5,0.1],[0.2,0.8],[0.9,1.0]]))
output_sample = ot.Sample(np.array([[1.0],[5.0],[4.2],[7.3],[4.5]]))

DIMENSION_SPATIALE = 2
AMPLITUDE = [1.0]
SCALE = [0.55,0.77]
NU = 2.5

#BASE = ot.ConstantBasisFactory(DIMENSION_SPATIALE).build()
BASE = ot.Basis(0)

noyau_tensorise = ot.ProductCovarianceModel([ot.MaternModel([SCALE[0]], AMPLITUDE, NU), 
                                             ot.MaternModel([SCALE[1]], AMPLITUDE, NU)])


g = Gu(input_sample, output_sample, BASE, noyau_tensorise, prior='reference', parametrization='standard')
g.set_scale(SCALE)
print("JOSEPH log likelihood = ", g._current_log_likelihood)
print("JOSEPH log prior = ", g._current_log_prior)
print("JOSEPH optimized function = ", g._current_log_likelihood + g._current_log_prior)
print("JOSEPH info fisher = ")
print(g.info_fisher_matern())
print("JOSEPH correlation matrix = ")
print(g.compute_matern_correlation_matrix())
print("JOSEPH inverse correlation matrix = ")
print(np.linalg.inv(g.compute_matern_correlation_matrix()))
print("JOSEPH derivative correlation matrix 0 = ")
print(g.compute_matern_derivative_matrix(0))
print("JOSEPH derivative correlation matrix 1 = ")
print(g.compute_matern_derivative_matrix(1))

noyau_tensorise = ot.ProductCovarianceModel([ot.MaternModel([SCALE[0]], AMPLITUDE, NU), 
                                             ot.MaternModel([SCALE[1]], AMPLITUDE, NU)])
noyau_tensorise.setScaleParametrization(noyau_tensorise.STANDARD)
algo = ot.KrigingAlgorithm(input_sample, output_sample, noyau_tensorise, BASE, False)
algo.setScalePrior(ot.GeneralLinearModelAlgorithm.REFERENCE)
fonction_optimisee = algo.getReducedLogLikelihoodFunction()
print("OT fonction optimisee = ", fonction_optimisee(SCALE))
