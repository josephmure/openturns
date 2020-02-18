#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 11:35:01 2019

@author: C39575
"""

import numpy as np
import openturns as ot

import re
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
from  Classe_Gu import *
ro.numpy2ri.activate()
importr('nloptr')
importr('RobustGaSP')


def fait_des_trucs():
    design = ro.r.matrix(np.identity(4))
    ro.r.assign('design', design)

class rgaspFromOT:
    def __init__(self, input_Sample, output_Sample, trend_Basis, initial_CovarianceModel, prior_choice='ref_gamma'):
        if type(input_Sample)!=ot.typ.Sample:
            raise TypeError("input_Sample should be an OpenTURNS Sample, not be of {}.".format(type(input_Sample)))
        if type(output_Sample)!=ot.typ.Sample:
            raise TypeError("output_Sample should be an OpenTURNS Sample, not be of {}.".format(type(input_Sample)))
        if output_Sample.getDimension()!=1:
            raise ValueError("output_Sample should be of dimension 1.")
        if output_Sample.getSize()!=input_Sample.getSize():
            message = "input_Sample has size {} and ".format(input_Sample.getSize())
            message += "output_Sample has size {}, ".format(output_Sample.getSize())
            message += "but they should have the same size."
            raise ValueError(message)
        if type(trend_Basis)!=ot.func.Basis:
            raise TypeError("trend_Basis should be an OpenTURNS Basis, not be of {}.".format(type(trend_Basis)))
        
        if type(initial_CovarianceModel)==ot.statistics.ProductCovarianceModel:
            collection = initial_CovarianceModel.getCollection()

            for model in collection:
                if model.getImplementation().getClassName()!='MaternModel':
                    raise TypeError("In a ProductCovarianceModel, all models should be Matern, not be of type {}.".format(model.getImplementation().getClassName()))
                if model.getFullParameter()[-1]!=collection[0].getFullParameter()[-1]:
                    raise TypeError("In a ProductCovarianceModel, all models should have the same nu.")
            
            self._nu = collection[0].getFullParameter()[-1]
        elif type(initial_CovarianceModel)!=ot.statistics.MaternModel:
            raise TypeError("initial_CovarianceModel should be an OpenTURNS MaternModel, not be of {}.".format(type(initial_CovarianceModel)))
        else:
            self._nu = initial_CovarianceModel.getNu()
         
        if initial_CovarianceModel.getScale().getSize()!=input_Sample.getDimension():
            initial_scale = str(list(np.array(initial_CovarianceModel.getScale())))
            message = "The scale parameter of initial_CovarianceModel is {}, ".format(initial_scale)
            message += "but input_Sample is of dimension {}. ".format(input_Sample.getDimension())
            message += "The dimensions do not match."
            raise ValueError(message)
            
        
        self._input_Sample = input_Sample # Should it be a deepcopy?
        self._output_Sample = output_Sample
        self._trend_Basis = trend_Basis
        self._current_CovarianceModel = initial_CovarianceModel
        
        
        self._prior_choice = prior_choice

        if self._prior_choice=='ref_approx':
            self._python_equivalent = Gu(self._input_Sample, 
                                         self._output_Sample, 
                                         self._trend_Basis, 
                                         self._current_CovarianceModel, 
                                         'jointly-robust', 'inverse')
        
        self._W, self._H = self.compute_W()

        
        self.assign_r_objects()
        
        
    def compute_W(self):
        if self._trend_Basis.getSize()==0:
            return np.identity(self._input_Sample.getSize()), 0
        else:
            if self._trend_Basis.getSize()>=self._input_Sample.getSize():
                message = "Size of trend_Basis {} >= ".format(self._trend_Basis.getSize())
                message += "size of input/output_Sample {}.".format(self._input_Sample.getSize())
                raise ValueError(message)
                
            matrice_H = np.array(self._trend_Basis.build(0)(self._input_Sample))
            for dim in range(1, self._trend_Basis.getSize()):
                matrice_H = np.concatenate((matrice_H, self._trend_Basis.build(dim)(self._input_Sample)), axis=1)
            
            u, s, vh = np.linalg.svd(matrice_H)
            return u[:, self._trend_Basis.getSize():], matrice_H
        
    def assign_r_objects(self):
        design = ro.r.matrix(np.array(self._input_Sample), nrow=self._input_Sample.getSize(), ncol=self._input_Sample.getDimension())
        ro.r.assign('design', design)
        
        response = ro.r.matrix(np.array(self._output_Sample), nrow=self._input_Sample.getSize(), ncol=1) #on a besoin d'une matrice colonne en entree
        ro.r.assign('response', response)
        
        if 'MaternModel' in str(self._current_CovarianceModel):
            nu = float(re.findall(r"[-+]?\d*\.\d+|\d+", str(self._current_CovarianceModel).split("nu=")[1])[0])
            if nu!=2.5 and nu!=1.5:
                raise ValueError("La régularité {!r} n'est pas gérée par le paquet R RobustGaSP.".format(self._current_CovarianceModel.getNu()))
            elif nu==2.5:
                kernel_type = 'matern_5_2'
            else:
                kernel_type = 'matern_3_2'
            alpha = 1 #c'est un placeholder : alpha n'est pas utilisé par le paquet R en cas de noyau de Matérn
        #elif 'GeneralizedExponential' in str(initial_CovarianceModel): # CANNOT HAPPEN BECAUSE WE ENSURED MATERN
        #    kernel_type = 'pow_exp'
        #    alpha = float(re.findall(r"[-+]?\d*\.\d+|\d+", str(initial_CovarianceModel).split("p=")[1])[0])
        #    alpha = np.array([alpha] * initial_CovarianceModel.getScale().getSize())
        else:
            raise TypeError("Ce noyau n'est pas géré par le paquet R RobustGaSP.")
        
        
        ro.r.assign('kernel_type', kernel_type)
        
        ro.r.assign('alpha', alpha)
        ro.r('alpha <- c(alpha)') #alpha doit etre un vecteur et non une matrice colonne        
        
        trend = ro.r.matrix(self._H, nrow=self._input_Sample.getSize(), ncol=self._trend_Basis.getSize())
        ro.r.assign('trend', trend)
        
        if self._trend_Basis.getSize()==0:
            zero_mean = 'Yes'
        else:
            zero_mean = 'No'
            
        ro.r.assign('zero.mean', zero_mean)
        
        ro.r.assign('prior_choice', self._prior_choice)
        
        current_range = ro.vectors.FloatVector(np.array(self._current_CovarianceModel.getScale()))
        ro.r.assign('current_range', current_range)
        
    def use_rgasp(self, fixed_scale=True):
        if fixed_scale:
            ro.r("range_par <- current_range")
        else:
            ro.r("range_par <- NA")
        commande_rgasp = """
        m1<- rgasp(design=design, 
                   response=response, 
                   trend=trend, #une matrice est specifiee meme en krigeage simple, et est ignoree si zero.mean==TRUE
                   zero.mean=zero.mean, #ce parametre decide du krigeage simple ou universel. Ce n'est PAS un booleen mais un string.
                   nugget=0, #version sans nugget effect de la programmation
                   nugget.est=F, #sans nugget : pas besoin de l'estimer
                   range.par=range_par, #on estime ou pas les longueurs de correlation.
                   prior_choice=prior_choice, #Attention, prior_choice specifie choix du prior (reference) ET parmetrisation (standard)
                   kernel_type=kernel_type,
                   alpha=alpha, #alpha est ignore en cas de noyau de Matern
                   num_initial_values=1) #pas de multi-start pour diminuer le temps d'execution
        """
        ro.r(commande_rgasp)
        new_scale = 1. / np.array(ro.r('m1@beta_hat'))
        self._current_CovarianceModel.setScale(new_scale)
        
    # Get slots of rgasp object
        
    def get_rgasp_p(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@p')[0]
    
    def get_rgasp_num_obs(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@num_obs')[0]
    
    def get_rgasp_input(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@input'))
    
    def get_rgasp_output(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@output'))    
    
    def get_rgasp_X(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@X'))
    
    def get_rgasp_zero_mean(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@zero_mean')[0]
    
    def get_rgasp_q(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@q')[0]
    
    def get_rgasp_LB(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@LB'))
    
    def get_rgasp_beta_initial(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@beta_initial'))

    def get_rgasp_beta_hat(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@beta_hat'))
    
    def get_rgasp_log_post(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@log_post')[0]    
    
    def get_rgasp_R0(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@R0')
    
    def get_rgasp_theta_hat(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@theta_hat'))
    
    def get_rgasp_L(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@L'))
    
    def get_rgasp_sigma2_hat(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@sigma2_hat')[0]        
    
    def get_rgasp_LX(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@LX'))    

    def get_rgasp_CL(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return np.array(ro.r('m1@CL'))
    
    def get_rgasp_nugget(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@nugget')[0]
    
    def get_rgasp_nugget_est(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@nugget.est')[0]
    
    def get_rgasp_kernel_type(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@kernel_type')[0]
    
    def get_rgasp_alpha(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@alpha')[0]    
    
    def get_rgasp_lower_bound(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@lower_bound')[0]    
    
    def get_rgasp_max_eval(self):
        self.assign_r_objects()
        #self.use_rgasp(True)
        return ro.r('m1@max_eval')[0] 
    
    # Get scale parameter directly from rgasp object
    def get_rgasp_scale(self):
        return 1./ self.get_rgasp_beta_hat()
    
    # Get jointly robust prior only
    def compute_rgasp_log_jointly_robust_prior(self):
        if self._prior_choice!='ref_approx':
            raise TypeError("compute_rgasp_log_jointly_robust_prior necessitates prior_choice=='ref_approx'.")
        self._python_equivalent.set_scale(self.get_rgasp_scale())
        return self.get_rgasp_log_post() - self._python_equivalent._current_log_likelihood
    
    # Get likelihood only (if jointly robust prior is used)
    def compute_rgasp_log_likelihood(self):
        if self._prior_choice!='ref_approx':
            raise TypeError("compute_rgasp_log_jointly_robust_prior necessitates prior_choice=='ref_approx'.")
        self._python_equivalent.set_scale(self.get_rgasp_scale())
        return self._python_equivalent._current_log_likelihood
  
    def use_log_approx_ref_prior(self):
        if self._prior_choice!='ref_approx':
            raise TypeError("compute_rgasp_log_jointly_robust_prior necessitates prior_choice=='ref_approx'.")
        ro.r.assign('b1', self._python_equivalent._b1)
        ro.r.assign('b', self._python_equivalent._b)
        self.assign_r_objects()
        return ro.r('log_approx_ref_prior(param=log(m1@beta_hat), nugget=m1@nugget, nugget_est=F, CL=m1@CL, a=b1, b=b)')[0]


NB_PTS_PLAN_EXPERIENCE = 5
DIMENSION_SPATIALE = 2

AMPLITUDE = [1.0]
SCALE = [.1,0.2]

ot.RandomGenerator.SetSeed(0)


loi_ech = ot.ComposedDistribution([ot.Uniform(0,1)] * DIMENSION_SPATIALE)
points = loi_ech.getSample(NB_PTS_PLAN_EXPERIENCE)
TREND = ot.SymbolicFunction(['x1', 'x2'],['5'])
BASE = ot.ConstantBasisFactory(DIMENSION_SPATIALE).build()

noyau_covariance = ot.MaternModel(SCALE, AMPLITUDE, 2.5)

realisation = ot.GaussianProcess(ot.TrendTransform(TREND,ot.Mesh(points)), noyau_covariance,ot.Mesh(points)).getRealization() 

noyau_tensorise = ot.ProductCovarianceModel([ot.MaternModel([1.0], AMPLITUDE, 2.5), ot.MaternModel([1.0], AMPLITUDE, 2.5)])

g = Gu(points, realisation.getValues(), BASE, noyau_tensorise, prior='reference', parametrization='standard')
g.set_scale([1.0,1.0])
g.optimize_scale()

noyau_tensorise = ot.ProductCovarianceModel([ot.MaternModel([1.0], AMPLITUDE, 2.5), ot.MaternModel([1.0], AMPLITUDE, 2.5)])
        
r = rgaspFromOT(points, realisation.getValues(), BASE, noyau_tensorise, prior_choice='ref_gamma')
r.use_rgasp(False)

noyau_tensorise = ot.ProductCovarianceModel([ot.MaternModel([1.0], AMPLITUDE, 2.5), ot.MaternModel([1.0], AMPLITUDE, 2.5)])
noyau_tensorise.setScaleParametrization(noyau_tensorise.STANDARD)

algo = ot.KrigingAlgorithm(points, realisation.getValues(), noyau_tensorise, BASE)
algo.setScalePrior(ot.GeneralLinearModelAlgorithm.REFERENCE)
algo.run()
model = algo.getResult().getCovarianceModel()
model.setScaleParametrization(noyau_tensorise.STANDARD)

print("Résultat formulé dans la param standard OT :", model.getScale())
print("Résultat formulé dans la param standard Gu_Python :", g._current_CovarianceModel.getScale())
print("Résultat formulé dans la param standard RobustGaSP_R :", r.get_rgasp_scale())


#NB_PTS_PLAN_EXPERIENCE = 20
#DIMENSION_SPATIALE = 1
#
#AMPLITUDE = [1.0]
#SCALE = [0.2]
#
#
#ot.RandomGenerator.SetSeed(3)
#
#loi_ech = ot.ComposedDistribution([ot.Uniform(0,1)] * DIMENSION_SPATIALE)
#points = loi_ech.getSample(NB_PTS_PLAN_EXPERIENCE)
#TREND = ot.SymbolicFunction(['x1'],['5'])
#BASE = ot.ConstantBasisFactory(DIMENSION_SPATIALE).build()
#
#noyau_covariance = ot.MaternModel(SCALE, AMPLITUDE, 2.5)
#realisation = ot.GaussianProcess(ot.TrendTransform(TREND,ot.Mesh(points)), noyau_covariance,ot.Mesh(points)).getRealization() 
#
#
#g = Gu(points, realisation.getValues(), ot.Basis(0), noyau_covariance)
#g.set_scale([0.5])
#res = g.optimize_scale()
#
#r = rgaspFromOT(points, realisation.getValues(), ot.Basis(0), noyau_covariance)

#NB_PTS_PLAN_EXPERIENCE = 20
#DIMENSION_SPATIALE = 3
#
#TREND = ot.SymbolicFunction(['x1','x2','x3'],['5'])
##TREND = ot.SymbolicFunction(['x'], ['0'])
#NU = 1.5
#
#AMPLITUDE = [2.0]
#SCALE = [0.1, 0.1, 0.1]
##SCALE = [.3]
#
#ot.RandomGenerator.SetSeed(1)
#
#loi_ech = ot.ComposedDistribution([ot.Uniform(0,1)] * DIMENSION_SPATIALE)
#points = loi_ech.getSample(NB_PTS_PLAN_EXPERIENCE)
#
#modele_covariance = ot.MaternModel(SCALE, AMPLITUDE, NU)
#sorties= ot.GaussianProcess(ot.TrendTransform(TREND ,ot.Mesh(points)), modele_covariance,ot.Mesh(points)).getRealization()
#
#g = Gu(points, sorties.getValues(), ot.Basis(0), modele_covariance)
#
#g._current_log_likelihood
