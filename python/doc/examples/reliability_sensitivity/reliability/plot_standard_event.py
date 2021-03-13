"""
Map a sample or an event to the standard space
==============================================

The FORM/SORM recipe has several ingredients,
but certainly one of the most fundamental is the transformation
from physical space (or X-space) to standard space (or U-space).
This is performed by the :class:`~openturns.StandardEvent` class.

The present example shows how the transformation affects
a multivariate normal distribution. It illustrates one of its main properties:
**rotational invariance**.

"""

# %%
import openturns as ot
from openturns.viewer import View
ot.Log.Show(ot.Log.NONE)

# %%
# Set a distribution in the physical space.

mu = [10.0] * 2
sigma = [3.0] * 2
R = ot.CorrelationMatrix([[1.0, 0.5], [0.5, 1.0]])
c = ot.NormalCopula(R)
X1 = ot.Normal(mu[0], sigma[0])
X1.setDescription([r"$x_1$"])
X2 = ot.Normal(mu[1], sigma[1])
X2.setDescription([r"$x_2$"])
distribution_X = ot.ComposedDistribution([X1, X2], c)

# %%
# Get the corresponding standard distribution.

distribution_U = distribution_X.getStandardDistribution()
distribution_U.setDescription([r"$u_1$", r"$u_2$"])

# %%
# Draw the PDF of the distribution in the physical space.

graph_X = distribution_X.drawPDF()
scientific_legend = ['{:.3e}'.format(float(x)) for x in graph_X.getLegends()]
graph_X.setLegends(scientific_legend)
legend_kw = {'bbox_to_anchor':(1,1), 'loc':"upper left"}
view = View(graph_X, square_axes=True, legend_kw=legend_kw)

# %%
# Draw the PDF of the corresponding standard distribution.
graph_U = distribution_U.drawPDF()
scientific_legend = ['{:.3e}'.format(float(x)) for x in graph_U.getLegends()]
graph_U.setLegends(scientific_legend)
view = View(graph_U, square_axes=True, legend_kw=legend_kw)

# %%
# Map a sample from the physical space to the standard space.

sample_X = distribution_X.getSample(10000)
transformation = distribution_X.getIsoProbabilisticTransformation()
sample_U = transformation(sample_X)

# %%
# Plot the original sample - the one in the physical space.
bb = graph_X.getBoundingBox()
cloud = ot.Cloud(sample_X)
cloud.setPointStyle("dot")
graph_X.add(cloud)
graph_X.setTitle("Sample in the physical space")
graph_X.setLegends([""])
low = bb.getLowerBound()[0]
up = bb.getUpperBound()[0]
view = View(graph_X, square_axes=True)._ax[0].axis([low, up, low, up])

# %%
# Plot the transformed sample - the one in the standard space.

# sphinx_gallery_thumbnail_number = 4
bb = graph_U.getBoundingBox()
cloud = ot.Cloud(sample_U)
cloud.setPointStyle("dot")
graph_U.add(cloud)
graph_U.setTitle("Sample in the standard space")
graph_U.setLegends([""])
low = bb.getLowerBound()[0]
up = bb.getUpperBound()[0]
view = View(graph_U, square_axes=True)._ax[0].axis([low, up, low, up])

# %%
# Define a limit state function.

formula = "2.5 - (x1 - x2) + (x1 + x2)^2"
limitStateFunction = ot.SymbolicFunction(["x1", "x2"], [formula])

# %%
# Set the event "the limit state function is nonpositive".
threshold = 0.0
random_vector = ot.RandomVector(distribution_X)
lsf_random_vector = ot.CompositeRandomVector(limitStateFunction, random_vector)
event = ot.ThresholdEvent(lsf_random_vector, ot.LessOrEqual(), threshold)

# %%
# Compute the corresponding event in the standard space.
myStandardEvent = ot.StandardEvent(event)
