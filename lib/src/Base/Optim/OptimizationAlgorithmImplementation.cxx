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

#include "openturns/OptimizationAlgorithmImplementation.hxx"
#include "openturns/PersistentObjectFactory.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(OptimizationAlgorithmImplementation)

static const Factory<OptimizationAlgorithmImplementation> Factory_OptimizationAlgorithmImplementation;

/* Default constructor */
OptimizationAlgorithmImplementation::OptimizationAlgorithmImplementation()
  : PersistentObject()
  , progressCallback_(std::make_pair<ProgressCallback, void *>(0, 0))
  , stopCallback_(std::make_pair<StopCallback, void *>(0, 0))
  , startingSample_(Sample())
  , maximumIterationNumber_(ResourceMap::GetAsUnsignedInteger("OptimizationAlgorithm-DefaultMaximumIterationNumber"))
  , maximumEvaluationNumber_(ResourceMap::GetAsUnsignedInteger("OptimizationAlgorithm-DefaultMaximumEvaluationNumber"))
  , maximumAbsoluteError_(ResourceMap::GetAsScalar("OptimizationAlgorithm-DefaultMaximumAbsoluteError"))
  , maximumRelativeError_(ResourceMap::GetAsScalar("OptimizationAlgorithm-DefaultMaximumRelativeError"))
  , maximumResidualError_(ResourceMap::GetAsScalar("OptimizationAlgorithm-DefaultMaximumResidualError"))
  , maximumConstraintError_(ResourceMap::GetAsScalar("OptimizationAlgorithm-DefaultMaximumConstraintError"))
{
  // Nothing to do
}

/*
 * @brief Standard constructor: the optimization problem is managed by the optimization solver, and the actual solver is in charge to check if it is able to solve it.
 */
OptimizationAlgorithmImplementation::OptimizationAlgorithmImplementation(const OptimizationProblem & problem)
  : PersistentObject()
  , progressCallback_(std::make_pair<ProgressCallback, void *>(0, 0))
  , stopCallback_(std::make_pair<StopCallback, void *>(0, 0))
  , startingSample_(Sample())
  , problem_(problem)
  , maximumIterationNumber_(ResourceMap::GetAsUnsignedInteger("OptimizationAlgorithm-DefaultMaximumIterationNumber"))
  , maximumEvaluationNumber_(ResourceMap::GetAsUnsignedInteger("OptimizationAlgorithm-DefaultMaximumEvaluationNumber"))
  , maximumAbsoluteError_(ResourceMap::GetAsScalar("OptimizationAlgorithm-DefaultMaximumAbsoluteError"))
  , maximumRelativeError_(ResourceMap::GetAsScalar("OptimizationAlgorithm-DefaultMaximumRelativeError"))
  , maximumResidualError_(ResourceMap::GetAsScalar("OptimizationAlgorithm-DefaultMaximumResidualError"))
  , maximumConstraintError_(ResourceMap::GetAsScalar("OptimizationAlgorithm-DefaultMaximumConstraintError"))
{
  // Nothing to do
}

/* Result accessor */
OptimizationResult OptimizationAlgorithmImplementation::getResult() const
{
  return result_;
}

/* Result accessor */
void OptimizationAlgorithmImplementation::setResult(const OptimizationResult & result)
{
  result_ = result;
}

/* Maximum iterations number accessor */
UnsignedInteger OptimizationAlgorithmImplementation::getMaximumIterationNumber() const
{
  return maximumIterationNumber_;
}

/* Maximum iterations number accessor */
void OptimizationAlgorithmImplementation::setMaximumIterationNumber(const UnsignedInteger maximumIterationNumber)
{
  maximumIterationNumber_ = maximumIterationNumber;
}

void OptimizationAlgorithmImplementation::setMaximumEvaluationNumber(const UnsignedInteger maximumEvaluationNumber)
{
  maximumEvaluationNumber_ = maximumEvaluationNumber;
}

UnsignedInteger OptimizationAlgorithmImplementation::getMaximumEvaluationNumber() const
{
  return maximumEvaluationNumber_;
}

/* Maximum absolute error accessor */
Scalar OptimizationAlgorithmImplementation::getMaximumAbsoluteError() const
{
  return maximumAbsoluteError_;
}

/* Maximum absolute error accessor */
void OptimizationAlgorithmImplementation::setMaximumAbsoluteError(const Scalar maximumAbsoluteError)
{
  maximumAbsoluteError_ = maximumAbsoluteError;
}

/* Maximum relative error accessor */
Scalar OptimizationAlgorithmImplementation::getMaximumRelativeError() const
{
  return maximumRelativeError_;
}

/* Maximum relative error accessor */
void OptimizationAlgorithmImplementation::setMaximumRelativeError(const Scalar maximumRelativeError)
{
  maximumRelativeError_ = maximumRelativeError;
}

/* Maximum residual error accessor */
Scalar OptimizationAlgorithmImplementation::getMaximumResidualError() const
{
  return maximumResidualError_;
}

/* Maximum residual error accessor */
void OptimizationAlgorithmImplementation::setMaximumResidualError(const Scalar maximumResidualError)
{
  maximumResidualError_ = maximumResidualError;
}

/* Maximum constraint error accessor */
Scalar OptimizationAlgorithmImplementation::getMaximumConstraintError() const
{
  return maximumConstraintError_;
}

/* Maximum constraint error accessor */
void OptimizationAlgorithmImplementation::setMaximumConstraintError(const Scalar maximumConstraintError)
{
  maximumConstraintError_ = maximumConstraintError;
}

/* String converter */
String OptimizationAlgorithmImplementation::__repr__() const
{
  OSS oss;
  oss << "class=" << OptimizationAlgorithmImplementation::GetClassName()
      << " problem=" << problem_
      << " startingSample=" << startingSample_
      << " maximumIterationNumber=" << maximumIterationNumber_
      << " maximumEvaluationNumber=" << maximumEvaluationNumber_
      << " maximumAbsoluteError=" << maximumAbsoluteError_
      << " maximumRelativeError=" << maximumRelativeError_
      << " maximumResidualError=" << maximumResidualError_
      << " maximumConstraintError=" << maximumConstraintError_
      << " verbose=" << verbose_;
      << " keepResults=" << keepResults_;
  return oss;
}

/* Problem accessor */
OptimizationProblem OptimizationAlgorithmImplementation::getProblem() const
{
  return problem_;
}

void OptimizationAlgorithmImplementation::setProblem(const OptimizationProblem & problem)
{
  const UnsignedInteger problemDimension = problem.getDimension();
  if ( (problemDimension > 0) && (startingSample_.getSize() > 0) ) // only perform check if problem is initalized and starting points are already defined
  {
    if (problemDimension != startingSample_.getDimension())
      throw InvalidArgumentException(HERE) << "Problem dimension (" << problemDimension
                                           << ") and starting points dimension (" << startingSample_.getDimension() << ") do not match.";
  } // checks

  checkProblem(problem);
  problem_ = problem;
}

/* Performs the actual checks. Must be overloaded by the actual optimisation algorithm */
void OptimizationAlgorithmImplementation::checkProblem(const OptimizationProblem & ) const
{
  throw NotYetImplementedException(HERE) << "In OptimizationAlgorithmImplementation::checkProblem()";
}

/* Performs the actual computation. Must be overloaded by the actual optimisation algorithm */
OptimizationResult OptimizationAlgorithmImplementation::runFromStartingPoint(const Point & startingPoint, const UnsignedInteger remainingEval)
{
  throw NotYetImplementedException(HERE) << "In OptimizationAlgorithmImplementation::runFromStartingPoint()";
}

/* Virtual constructor */
OptimizationAlgorithmImplementation * OptimizationAlgorithmImplementation::clone() const
{
  return new OptimizationAlgorithmImplementation(*this);
}

/* Verbose accessor */
Bool OptimizationAlgorithmImplementation::getVerbose() const
{
  return verbose_;
}

/* Verbose accessor */
void OptimizationAlgorithmImplementation::setVerbose(const Bool verbose)
{
  verbose_ = verbose;
}

void OptimizationAlgorithmImplementation::appendStartingSample(const Point & startingPoint)
{
  const UnsignedInteger problemDimension = getProblem().getDimension();
  if (problemDimension > 0) // only perform the check if the problem has been set
  {
    if (problemDimension != startingPoint.getDimension())
      throw InvalidArgumentException(HERE) << "Proposed starting point has dimension " << startingPoint.getDimension()
                                           << ", but the optimization problem has dimension " << problemDimension;
  }

  if (startingSample_.getSize() == 0) startingSample_ = Sample(1, startingPoint);
  else startingSample_.add(startingPoint);
}

void OptimizationAlgorithmImplementation::appendStartingSample(const Sample & startingSample)
{
  const UnsignedInteger problemDimension = getProblem().getDimension();
  if (problemDimension > 0) // only perform the check if the problem has been set
  {
    if (problemDimension != startingSample.getDimension())
      throw InvalidArgumentException(HERE) << "Proposed starting sample has dimension " << startingSample.getDimension()
                                           << ", but the optimization problem has dimension " << problemDimension;
  }

  if (startingSample_.getSize() == 0) startingSample_ = startingSample;
  else startingSample_.add(startingSample);
}

/* Automatically select starting points */
void OptimizationAlgorithmImplementation::generateAdditionalStartingPoints(const UnsignedInteger nbStartingSampleToBeSelected, const LowDiscrepancySequence & generator)
{
  // Select points only if there remain points to be selected and the problem has been properly set.
  if ((nbStartingSampleToBeSelected > 0) && (getProblem().getDimension() > 0))
  {
    LOGINFO(OSS() << "Use multi-start with base optimization algorithm=" << solver_.getImplementation()->getClassName());
    const Point lowerBound(getProblem().getBounds().getLowerBound());
    const Point upperBound(getProblem().getBounds().getUpperBound());
    const UnsignedInteger dimension = lowerBound.getDimension();
    generator.initialize(dimension);
    Sample population(nbStartingSampleToBeSelected, dimension);
    for (UnsignedInteger i = 0; i < nbStartingSampleToBeSelected; ++i)
      {
        const Point u(generator.generate());
        for (UnsignedInteger j = 0; j < dimension; ++j)
          population(i, j) = lowerBound[j] + (upperBound[j] - lowerBound[j]) * u[j];
      } // i

    if (startingSample_.getSize() == 0) startingSample_ = population;
    else startingSample_.add(population);
  }
}

void OptimizationAlgorithmImplementation::generateAdditionalStartingPoints(Experiment & experiment)
{
  addToStartingSample(experiment.generate());
}

void OptimizationAlgorithmImplementation::run()
{
  if (startingSample_.getSize() == 0) throw InvalidArgumentException(HERE) << "No starting points are set.";
  const UnsignedInteger problemDimension = getProblem().getDimension();
  if (problemDimension == 0) throw InvalidArgumentException(HERE) << "No problem has been set.";
  if (problemDimension != startingSample_.getDimension())
    throw InvalidArgumentException(HERE) << "The starting points dimension (" << startingSample_.getDimension()
                                         << ") and the problem dimension (" << problemDimension << ") do not match.";

  // run the solver with each starting point
  resultCollection_.clear();
  Scalar bestValue = getProblem().isMinimization() ? SpecFunc::MaxScalar : SpecFunc::LowestScalar;
  const UnsignedInteger size = startingSample_.getSize();
  const UnsignedInteger initialEvaluationNumber = getProblem().getObjective().getEvaluationCallsNumber();
  UnsignedInteger evaluationNumber = 0;
  UnsignedInteger successNumber = 0;
  UnsignedInteger improvementNumber = 0;
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    // ensure we do not exceed the global budget if the maximum eval number is set
    const UnsignedInteger remainingEval = std::max(static_cast<SignedInteger>(getMaximumEvaluationNumber() - evaluationNumber), 0L);
    LOGDEBUG(OSS() << "Working with starting point[" << i << "]=" << startingSample_[i] << ", " << remainingEval << " remaining evaluations");

    const OptimizationResult result();
    try
    {
      result = runFromStartingPoint(startingSample_[i], remainingEval);
      ++successNumber;
    }
    catch (Exception & ex)
    {
      LOGDEBUG(OSS() << "StartingPoint " << i << " failed. Reason=" << ex);
      continue;
    }

    if (keepResults_) resultCollection_.add(result);
    Scalar currentValue = result.getOptimalValue()[0];
    if ((getProblem().isMinimization() && (currentValue < bestValue))
        || (!getProblem().isMinimization() && (currentValue > bestValue)))
    {
      bestValue = currentValue;
      setResult(result);
      LOGINFO(OSS() << "Best initial point so far=" << result.getOptimalPoint() << " value=" << result.getOptimalValue());
      ++improvementNumber;
    }

    evaluationNumber += getProblem().getObjective().getEvaluationCallsNumber() - initialEvaluationNumber;
    LOGDEBUG(OSS() << "Number of evaluations so far=" << evaluationNumber);
    if (evaluationNumber > getMaximumEvaluationNumber())
    {
      break;
    }

    // callbacks
    if (progressCallback_.first)
    {
      progressCallback_.first((100.0 * evaluationNumber) / getMaximumEvaluationNumber(), progressCallback_.second);
    }
    if (stopCallback_.first)
    {
      Bool stop = stopCallback_.first(stopCallback_.second);
      if (stop)
      {
        LOGWARN(OSS() << "Multi-start was stopped by user");
        break;
      }
    }
  }
  LOGINFO(OSS() << successNumber << " out of " << size << " local searches succeeded, " << improvementNumber << " improvements");

  if (successNumber == 0)
  {
    throw InternalException(HERE) << "None of the local searches succeeded.";
  }
}

/* Starting sample accessor */
void OptimizationAlgorithmImplementation::setStartingSample(const Sample & startingSample, const Bool resetProblem)
{
  if (resetProblem) setProblem(OptimizationProblem());

  const UnsignedInteger problemDimension = getProblem().getDimension();
  if (problemDimension > 0) // only perform the check if the problem has been set
  {
    if (problemDimension != startingSample.getDimension())
      throw InvalidArgumentException(HERE) << "Proposed starting sample has dimension " << startingSample.getDimension()
                                           << ", but the optimization problem has dimension " << problemDimension;
  }

  startingSample_ = startingSample;
}

/* Starting points accessor */
Sample OptimizationAlgorithmImplementation::getStartingSample() const
{
  return startingSample_;
}


/* Flag for results management accessors */
Bool OptimizationAlgorithmImplementation::getKeepResults() const
{
  return keepResults_;
}

void OptimizationAlgorithmImplementation::setKeepResults(const Bool keepResults)
{
  keepResults_ = keepResults;
}

OptimizationAlgorithmImplementation::OptimizationResultCollection OptimizationAlgorithmImplementation::getResultCollection() const
{
  return resultCollection_;
}

/* Method save() stores the object through the StorageManager */
void OptimizationAlgorithmImplementation::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute( "startingSample_", startingSample_);
  adv.saveAttribute( "problem_", problem_);
  adv.saveAttribute( "maximumIterationNumber_", maximumIterationNumber_);
  adv.saveAttribute( "maximumEvaluationNumber_", maximumEvaluationNumber_);
  adv.saveAttribute( "maximumAbsoluteError_", maximumAbsoluteError_);
  adv.saveAttribute( "maximumRelativeError_", maximumRelativeError_);
  adv.saveAttribute( "maximumResidualError_", maximumResidualError_);
  adv.saveAttribute( "maximumConstraintError_", maximumConstraintError_);
  adv.saveAttribute( "verbose_", verbose_);
  adv.saveAttribute("keepResults_", keepResults_);
  adv.saveAttribute("resultCollection_", resultCollection_);
}


/* Method load() reloads the object from the StorageManager */
void OptimizationAlgorithmImplementation::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute( "startingSample_", startingSample_);
  adv.loadAttribute( "problem_", problem_);
  adv.loadAttribute( "maximumIterationNumber_", maximumIterationNumber_);
  adv.loadAttribute( "maximumEvaluationNumber_", maximumEvaluationNumber_);
  adv.loadAttribute( "maximumAbsoluteError_", maximumAbsoluteError_);
  adv.loadAttribute( "maximumRelativeError_", maximumRelativeError_);
  adv.loadAttribute( "maximumResidualError_", maximumResidualError_);
  adv.loadAttribute( "maximumConstraintError_", maximumConstraintError_);
  adv.loadAttribute( "verbose_", verbose_);
  adv.loadAttribute("keepResults_", keepResults_);
  adv.loadAttribute("resultCollection_", resultCollection_);
}


void OptimizationAlgorithmImplementation::setProgressCallback(ProgressCallback callBack, void * state)
{
  progressCallback_ = std::pair<ProgressCallback, void *>(callBack, state);
}


void OptimizationAlgorithmImplementation::setStopCallback(StopCallback callBack, void * state)
{
  stopCallback_ = std::pair<StopCallback, void *>(callBack, state);
}


END_NAMESPACE_OPENTURNS
