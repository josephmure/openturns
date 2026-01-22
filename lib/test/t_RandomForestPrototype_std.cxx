//                                               -*- C++ -*-
/**
 *  @brief The test file of class Triangular for standard methods
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
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
#include "openturns/OT.hxx"
#include "openturns/OTtestcode.hxx"

using namespace OT;
using namespace OT::Test;


int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);

  try
  {
    // Test basic functionnalities
    SymbolicFunction model(Description({"x", "y"}), Description({"x * sin(y)"}));
    Uniform uniform;
    JointDistribution joint({uniform, uniform});
    Distribution distribution(joint);

    UnsignedInteger samplingSize = 25;
    MonteCarloExperiment experiment(distribution, samplingSize);
    Sample inputSample(experiment.generate());
    Sample outputSample(model(inputSample));

    RandomForestPrototype rf(inputSample, outputSample);
    fullprint << "Objet rf construit " << std::endl;
    rf.run();
    fullprint << "Entrainement termine " << std::endl;
    Sample prediction(rf.predict(inputSample));
    fullprint << "Prediction terminee " << std::endl;
    fullprint << "Prediction          =" << prediction << std::endl;
    fullprint << "Reference           =" << outputSample << std::endl;
    fullprint << "Residuals           =" << prediction - outputSample << std::endl;
    fullprint << "R2                  =" << 1.0 - (prediction - outputSample).asPoint().norm() / outputSample.computeStandardDeviation()[0] / outputSample.getSize() << std::endl;


  }
  catch (const TestFailed & ex)
  {
    std::cerr << ex << std::endl;
    return ExitCode::Error;
  }


  return ExitCode::Success;
}
