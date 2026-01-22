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
#ifndef OPENTURNS_RANDOMFORESTPROTOTYPE_HXX
#define OPENTURNS_RANDOMFORESTPROTOTYPE_HXX

#include "openturns/MetaModelAlgorithm.hxx"
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

class OT_API RandomForestPrototype
  : public MetaModelAlgorithm
{
  CLASSNAME

public:

  /** Default constructor */
  RandomForestPrototype();

  /** Parameters constructor */
  RandomForestPrototype(const Sample & inputSample,
                      const Sample & outputSample);

  /** Virtual constructor */
  RandomForestPrototype* clone() const override;

  void run() override;
  Sample predict(Sample& inputSample);
    
private:

#ifdef OPENTURNS_HAVE_RANGER
  // convertToRangerData(const Sample& input, const Sample& output);
  // std::shared_ptr<ranger::ForestRegression> forest_ = 0;
  UnsignedInteger num_trees_;
  std::vector<std::vector<long unsigned int>> split_var_ids_;
  std::vector<std::vector<std::vector<long unsigned int>>> child_node_ids_;
  std::vector<std::vector<double>> split_values_;
  std::vector<bool> is_ordered_variable_;

  class DataRanger
    : public ranger::Data
    {
      public:
        DataRanger() = default;
        DataRanger(const Sample& input, const Sample& output)
          :ranger::Data()
          // ,num_rows(input.getSize())
          // ,num_rows_rounded(0)
          // ,num_cols(input.getSize() + 1)
          // ,snp_data(0)
          // ,num_cols_no_snp(input.getSize() + 1)
          // ,externalData(true)
          // ,index_data(0)
          // ,max_num_unique_values(0)
          // ,order_snps(false)
          // ,any_na(false) 
          ,x(input)
          ,y(output)
        {
          num_rows = input.getSize();
          num_cols = input.getDimension();
          num_cols_no_snp = num_cols;
          const Description inputDescription(input.getDescription());
          const Description outputDescription(output.getDescription());
          for (size_t i = 0; i < input.getDimension(); ++i) {
              variable_names.push_back("var_" + inputDescription[i]);
          }
          variable_names.push_back("var_" + outputDescription[0]);
        };
        
        DataRanger(const DataRanger&) = delete;
        DataRanger& operator=(const DataRanger&) = delete;

        virtual ~DataRanger() override = default;

        double get_x(size_t row, size_t col) const override {
          // Use permuted data for corrected impurity importance
          //LOGWARN(OSS() << "Dans get_x(" << row << ", "<< col << ")");
          //LOGWARN(OSS() << "valeur = " << x(row, col));
          size_t col_permuted = col;
          if (col >= num_cols) {
            col = getUnpermutedVarID(col);
            row = getPermutedSampleID(row);
          }

          if (col < num_cols_no_snp) {
            return x(row, col);
          } else {
            throw InvalidArgumentException(HERE) << "Outside the table";
          }
        }

        double get_y(size_t row, size_t col) const override {
          //LOGWARN(OSS() << "Dans get_y(" << row << ", "<< col << ")");
          return y(row, col);
        }

        void reserveMemory(size_t y_cols) override {
          //LOGWARN(OSS() << "Dans reserveMemory avec y_cols = " << y_cols);
          if (y_cols != 1) throw InvalidArgumentException(HERE) << "Only 1 output dimension possible";
        }

        void set_x(size_t col, size_t row, double value, bool& error) override {
          //LOGWARN(OSS() << "Dans set_x(" << row << ", "<< col << ") avec value = " << value);
          x(row, col) = value;
        }

        void set_y(size_t col, size_t row, double value, bool& error) override {
          //LOGWARN(OSS() << "Dans set_y(" << row << ", "<< col << ") avec value = " << value);
          y(row, col) = value;
        }

      private:
        Sample x;
        Sample y;
      };
#endif
};

END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_RANDOMFORESTPROTOTYPE_HXX */
