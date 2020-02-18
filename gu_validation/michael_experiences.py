#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 11:37:14 2020

@author: osboxes
"""

import openturns as ot

sampleSize = 6
dimension = 1

f = ot.SymbolicFunction(['x0'], ['x0 * sin(x0)'])

X = ot.Sample(sampleSize, dimension)
X2 = ot.Sample(sampleSize, dimension)
for i in range(sampleSize):
    X[i, 0] = 3.0 + i
    X2[i, 0] = 2.5 + i
X[0, 0] = 1.0
X[1, 0] = 3.0
Y = f(X)
print(X)

# Version OpenTURNS/Gu
prior = ot.GeneralLinearModelAlgorithm.REFERENCE
parametrization = ot.CovarianceModelImplementation.STANDARD
basis = ot.ConstantBasisFactory(dimension).build()
covarianceModel = ot.MaternModel([1], [4.123456], 0.5)
covarianceModelCopy = ot.MaternModel(covarianceModel)
covarianceModel.setScaleParametrization(parametrization)

algo = ot.KrigingAlgorithm(X, Y, covarianceModel, basis)
algo.setScalePrior(prior)
algo.run()

# perform an evaluation
result = algo.getResult()
print(' ')
print('prior = ', prior)
print('parametrization = ', parametrization)
scaleOTGU = result.getCovarianceModel().getParameter()
print('OpenTURNS/Gu : parameters=', scaleOTGU)
rllfunction = algo.getReducedLogLikelihoodFunction()
objective_value = rllfunction(scaleOTGU)
print("Objectif value=", objective_value)

# Version Joseph/Gu

essai = Gu(X, Y, basis, covarianceModelCopy, 'reference', 'standard')
result = essai.optimize_scale()
scaleJOGU = result.x
print("Joseph/Gu : scale=", scaleJOGU)
objective_value = essai.set_scale(scaleJOGU)
print("Objectif value sur Optimum Joseph/Gu=", objective_value)

# Validation de la fonction objectif
scale = ot.Point([0.5])
print("Objectif value sur Optimum OpenTURNS/Gu=", rllfunction(scale))
print("Objectif value sur Optimum Joseph/Gu=", essai.set_scale(scale))
