%feature("docstring") OT::RandomForestPrototype
"First attempt at a random forest class


Notes
-----
This class wraps classes from the C++ part of the Ranger package.

Examples
--------
Create the model:

>>> import openturns as ot
>>> rf = ot.RandomForestPrototype()
>>> inputDimension = 1
>>> model = ot.SymbolicFunction(['x'], ['x * sin(x)'])
>>> distribution = ot.JointDistribution([ot.Uniform()] * inputDimension)

Define the evaluation strategy of the  coefficients:

>>> samplingSize = 50
>>> experiment = ot.MonteCarloExperiment(distribution, samplingSize)
>>> inputSample = experiment.generate()
>>> outputSample = model(inputSample)
>>> rf.train(inputSample, outputSample)
>>> rf.predict(inputSample)

Test it:

>>> X = [[0.5]]
>>> print(rf.predict(X))
[0.239713]"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomForestPrototype::train
"Train the random forest.

Parameters
----------
inputSample : :class:`~openturns.Sample`
    Input sample.

outputSample : :class:`~openturns.Sample`
    Output sample."

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomForestPrototype::predict
"Predict with the random forest.

Parameters
----------
inputSample : :class:`~openturns.Sample`
    Input sample.

Returns
-------
outputSample : :class:`~openturns.Sample`
    Prediction."
