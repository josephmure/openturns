//                                               -*- C++ -*-
/**
 *  @brief Implement Ranger wrapping for random forests
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

#include "openturns/RandomForestPrototype.hxx"

#ifdef OPENTURNS_HAVE_RANGER
#include <Forest.h>
#include <ForestRegression.h>
#include <Data.h>
#include <DataDouble.h>
#endif

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class RandomForestPrototype
 */

CLASSNAMEINIT(RandomForestPrototype)

/* Default constructor */
RandomForestPrototype::RandomForestPrototype()
  : MetaModelAlgorithm()
  , num_trees_(0)
  , child_node_ids_()
  , split_values_()
  , is_ordered_variable_()
{
  // Nothing to do
}

RandomForestPrototype::RandomForestPrototype(const Sample & inputSample,
    const Sample & outputSample)
  : MetaModelAlgorithm(inputSample, outputSample)
  , num_trees_(0)
  , child_node_ids_()
  , split_values_()
  , is_ordered_variable_()
{
  // Nothing to do
}


/* Virtual constructor */
RandomForestPrototype * RandomForestPrototype::clone() const
{
  return new RandomForestPrototype(*this);
}

void RandomForestPrototype::run()
{
#ifdef OPENTURNS_HAVE_RANGER

  auto data = std::make_unique<RandomForestPrototype::DataRanger>(inputSample_, outputSample_);
  data->setIsOrderedVariable({}); // no unordered variable

  is_ordered_variable_ = data->getIsOrderedVariable();

  // Create forest
  std::shared_ptr<ranger::ForestRegression> forest(
      new ranger::ForestRegression()
  );

  // Initialize forest
  std::vector<uint> zero_uint_vector = {0};
  std::vector<std::vector<long unsigned int>> empty_uint_sample = {};
  std::vector<std::vector<double>> empty_sample = {};
  std::vector<double> empty_point = {};
  std::vector<double> sample_fraction = {1.0}; // valeur par defaut en cas de remplacement
  forest->initR(
      std::move(data),
      /* mtry */ 0,  // 0 = auto-select
      /* num_trees */ 500,
      /* verbose_out */ &std::cout,
      /* seed */ 42,
      /* num_threads */ 4,
      /* importance_mode */ ranger::IMP_NONE,
      /* min_node_size */ zero_uint_vector, // {0} = auto-select
      /* min_bucket */ zero_uint_vector, // {0} = auto-select
      /* split_select_weights */ empty_sample, // {} = desactivation
      /* always_split_variable_names */ {},
      /* prediction_mode */ false,
      /* sample_with_replacement */ true,
      /* unordered_variable_names */ {},
      /* memory_saving_splitting */ false,
      /* splitrule */ ranger::DEFAULT_SPLITRULE,
      /* case_weights */ empty_point,
      /* manual_inbag */ empty_uint_sample,
      /* predict_all */ false,
      /* keep_inbag */ false,
      /* sample_fraction */ sample_fraction, 
      /* alpha */ 0.5, // sans importance
      /* minprop */ 0.1, // sans importance
      /* poisson_tau */ 0.1, //sans_importance
      /* holdout */ false,
      /* prediction_type */ ranger::RESPONSE,
      /* num_random_splits */ 1, // sans importance
      /* order_snps */ false, // sans importance
      /* max_depth */ 0, // no maximum depth
      /* regularization_factor */ {}, // no regularization
      /* regularization_usedepth */ false, // inutile sans regularisation
      /* node_stats */ false
  );

  // Train
  forest->run(/* verbose */ false, /* compute_oob_error */ true);

  child_node_ids_ = forest->getChildNodeIDs();
  split_var_ids_ = forest->getSplitVarIDs();
  split_values_ = forest->getSplitValues();
  num_trees_ = forest->getNumTrees();
        
#else
        throw NotYetImplementedException(HERE) 
            << "Random forest requires Ranger library";
#endif
}


Sample RandomForestPrototype::predict(Sample& inputSample)
{
#ifdef OPENTURNS_HAVE_RANGER

  Sample syntheticOutputSample(inputSample.getSize(), outputSample_.getDimension());
  syntheticOutputSample.setDescription(outputSample_.getDescription());
  inputSample.setDescription(inputSample_.getDescription());
  auto data = std::make_unique<RandomForestPrototype::DataRanger>(inputSample, syntheticOutputSample);

  // Create forest
  std::shared_ptr<ranger::ForestRegression> forest(
      new ranger::ForestRegression()
  );

  // Initialize forest
  std::vector<uint> zero_uint_vector = {0};
  std::vector<std::vector<long unsigned int>> empty_uint_sample = {};
  std::vector<std::vector<double>> empty_sample = {};
  std::vector<double> empty_point = {};
  std::vector<double> sample_fraction = {1.0}; // valeur par defaut en cas de remplacement

  forest->initR(
      std::move(data),
      /* mtry */ 0,  // 0 = auto-select
      /* num_trees */ 500,
      /* verbose_out */ &std::cout,
      /* seed */ 42,
      /* num_threads */ 1,
      /* importance_mode */ ranger::IMP_NONE,
      /* min_node_size */ zero_uint_vector, // {0} = auto-select
      /* min_bucket */ zero_uint_vector, // {0} = auto-select
      /* split_select_weights */ empty_sample, // {} = desactivation
      /* always_split_variable_names */ {},
      /* prediction_mode */ true,
      /* sample_with_replacement */ true,
      /* unordered_variable_names */ {},
      /* memory_saving_splitting */ false,
      /* splitrule */ ranger::DEFAULT_SPLITRULE,
      /* case_weights */ empty_point,
      /* manual_inbag */ empty_uint_sample,
      /* predict_all */ false,
      /* keep_inbag */ false,
      /* sample_fraction */ sample_fraction, 
      /* alpha */ 0.5, // sans importance
      /* minprop */ 0.1, // sans importance
      /* poisson_tau */ 0.1, //sans_importance
      /* holdout */ false,
      /* prediction_type */ ranger::RESPONSE,
      /* num_random_splits */ 1, // sans importance
      /* order_snps */ false, // sans importance
      /* max_depth */ 0, // no maximum depth
      /* regularization_factor */ {}, // no regularization
      /* regularization_usedepth */ false, // inutile sans regularisation
      /* node_stats */ false
  );

  forest->loadForest(
    num_trees_,
    child_node_ids_,
    split_var_ids_,
    split_values_,
    is_ordered_variable_
  );

  //LOGWARN(OSS() << "Dans predict : apres le loadForest ");

  // Predict
  forest->run(/* verbose */ true, /* compute_oob_error */ false);

  //LOGWARN(OSS() << "Dans predict : apres le run ");


  std::vector<std::vector<double>> predictions(forest->getPredictions()[0]);

  Sample result(predictions[0].size(), predictions.size());
  result.setDescription(outputSample_.getDescription());
  for (UnsignedInteger i = 0; i < result.getSize(); ++i)
  {
    for (UnsignedInteger j = 0; j < result.getDimension(); ++j)
    {
      result(i, j) = predictions[j][i];
    }
  }

  return result;
        
#else
        throw NotYetImplementedException(HERE) 
            << "Random forest requires Ranger library";
#endif
}

END_NAMESPACE_OPENTURNS
