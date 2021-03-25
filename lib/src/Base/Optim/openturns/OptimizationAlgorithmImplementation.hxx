//                                               -*- C++ -*-
/**
 *  @brief OptimizationAlgorithmImplementation implements an algorithm for solving an optimization problem
 *
 *  Copyright 2005-2021 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef OPENTURNS_OPTIMIZATIONALGORITHMIMPLEMENTATION_HXX
#define OPENTURNS_OPTIMIZATIONALGORITHMIMPLEMENTATION_HXX

#include "openturns/OTprivate.hxx"
#include "openturns/PersistentObject.hxx"
#include "openturns/OptimizationProblem.hxx"
#include "openturns/OptimizationResult.hxx"
#include "openturns/LowDiscrepancySequence.hxx"
#include "openturns/SobolSequence.hxx"
#include "openturns/Experiment.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class OptimizationAlgorithmImplementation
 * OptimizationAlgorithmImplementation implements an algorithm for solving an optimization problem
 */

class OT_API OptimizationAlgorithmImplementation
  : public PersistentObject
{

  CLASSNAME
public:

  typedef OT::Collection<OT::OptimizationResult>           OptimizationResultCollection;
  typedef OT::PersistentCollection<OT::OptimizationResult> OptimizationResultPersistentCollection;

  /** Default constructor */
  OptimizationAlgorithmImplementation();

  /** Constructor with parameters */
  explicit OptimizationAlgorithmImplementation(const OptimizationProblem & problem);

  /** Virtual constructor */
  OptimizationAlgorithmImplementation * clone() const override;

  /** Performs the actual computation */
  virtual void run();

  /** Starting sample accessor */
  virtual Sample getStartingSample() const;
  virtual void setStartingSample(const Sample & startingSample, const Bool resetProblem);

  /** Append starting sample */
  virtual void appendStartingSample(const Point & startingPoint);
  virtual void appendStartingSample(const Sample & startingSample);

  /** Automatically select starting points */
  virtual void generateAdditionalStartingPoints(const UnsignedInteger nbStartingSampleToBeSelected, const LowDiscrepancySequence & generator = SobolSequence());
  virtual void generateAdditionalStartingPoints(Experiment & experiment);

  /** Problem accessor */
  virtual OptimizationProblem getProblem() const;
  virtual void setProblem(const OptimizationProblem & problem);

  /** Result accessor */
  virtual OptimizationResult getResult() const;

  /** Flag for results management accessors */
  virtual Bool getKeepResults() const;
  virtual void setKeepResults(const Bool keepResults);
  virtual OptimizationResultCollection getResultCollection() const;

  /** Result accessor */
  virtual void setResult(const OptimizationResult & result);

  /** Maximum iterations number accessor */
  virtual void setMaximumIterationNumber(const UnsignedInteger maximumIterationNumber);
  virtual UnsignedInteger getMaximumIterationNumber() const;

  /** Maximum evaluations number accessor */
  virtual void setMaximumEvaluationNumber(const UnsignedInteger maximumEvaluationNumber);
  virtual UnsignedInteger getMaximumEvaluationNumber() const;

  /** Maximum absolute error accessor */
  virtual Scalar getMaximumAbsoluteError() const;

  /** Maximum absolute error accessor */
  virtual void setMaximumAbsoluteError(const Scalar maximumAbsoluteError);

  /** Maximum relative error accessor */
  virtual Scalar getMaximumRelativeError() const;

  /** Maximum relative error accessor */
  virtual void setMaximumRelativeError(const Scalar maximumRelativeError);

  /** Maximum residual error accessor */
  virtual Scalar getMaximumResidualError() const;

  /** Maximum residual error accessor */
  virtual void setMaximumResidualError(const Scalar maximumResidualError);

  /** Maximum constraint error accessor */
  virtual Scalar getMaximumConstraintError() const;

  /** Maximum constraint error accessor */
  virtual void setMaximumConstraintError(const Scalar maximumConstraintError);

  /** String converter */
  String __repr__() const override;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;

  /** Verbose accessor */
  virtual Bool getVerbose() const;
  virtual void setVerbose(const Bool verbose);

  /** Progress callback */
  typedef void (*ProgressCallback)(Scalar, void * state);
  virtual void setProgressCallback(ProgressCallback callBack, void * state = 0);

  /** Stop callback */
  typedef Bool (*StopCallback)(void * state);
  virtual void setStopCallback(StopCallback callBack, void * state = 0);

protected:
  /** Check whether this problem can be solved by this solver.  Must be overloaded by the actual optimisation algorithm */
  virtual void checkProblem(const OptimizationProblem & problem) const;

  /** Run optimization algorithm from a specific point */
  virtual OptimizationResult runFromStartingPoint(const Point & startingPoint, const UnsignedInteger remainingEval);

  /** The result of the algorithm */
  OptimizationResult result_;

  // callbacks
  std::pair< ProgressCallback, void *> progressCallback_;
  std::pair< StopCallback, void *> stopCallback_;

private:
  Point startingSample_;
  OptimizationProblem problem_;

  /** Number of outermost iterations (in case of nested iterations) */
  UnsignedInteger maximumIterationNumber_;

  /** Maximum function calls */
  UnsignedInteger maximumEvaluationNumber_;

  Scalar maximumAbsoluteError_;    /**< Value of ||x_n - x_{n-1}|| */
  Scalar maximumRelativeError_;    /**< Value of ||x_n - x_{n-1}|| / ||x_n|| */
  Scalar maximumResidualError_;    /**< Value of ||objectiveFunction(x_n) - objectiveFunction(x_{n-1})|| */
  Scalar maximumConstraintError_;  /**< Value of ||constraints(x_n)|| for the active constraints */
  Bool verbose_ = false;

  /** Flag to tell if the collection of optimization results have to be kept */
  Bool keepResults_;
  OptimizationResultPersistentCollection resultCollection_;

} ; /* class OptimizationAlgorithmImplementation */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_OPTIMIZATIONALGORITHMIMPLEMENTATION_HXX */
