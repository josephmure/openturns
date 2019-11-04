//                                               -*- C++ -*-
/**
 *  @brief The test file of class CorrectedLeaveOneOut for standard methods
 *
 *  Copyright 2005-2019 Airbus-EDF-IMACS-ONERA-Phimeca
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
  setRandomGenerator();

  UnsignedInteger size = 100;
  Uniform xuniform(0.9, 1.1);
  Sample x( xuniform.getSample(size) );
  Uniform yuniform(1.9, 2.1);
  Sample y( yuniform.getSample(size) );
  SymbolicFunction f( Description(1, "x"), Description(1, "2.0*x") );
  Collection<Function> basis(1, f);
  Indices indices(1);
  indices.fill();
  FittingAlgorithm fittingAlgo = CorrectedLeaveOneOut();

  Scalar result = fittingAlgo.run(x, y, basis, indices);

  fullprint << "result = " << result << std::endl;


  try
  {
    size = 2;
    x = xuniform.getSample(size);
    y = yuniform.getSample(size);
    SymbolicFunction f2( Description(1, "x"), Description(1, "x^2") );
    basis.add(f2);
    SymbolicFunction f3( Description(1, "x"), Description(1, "x^3") );
    basis.add(f3);
    indices = Indices(basis.getSize());
    indices.fill();
    result = fittingAlgo.run( x, y, basis, indices );
    fullprint << "result=" << result << std::endl;
  }
  catch (InvalidArgumentException & ex)
  {
    fullprint << "caught exception: " << ex.what() << std::endl;
  }

  return ExitCode::Success;
}
