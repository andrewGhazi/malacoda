/*
    malacoda is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    malacoda is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with malacoda.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.17.0

#include <stan/model/model_header.hpp>

namespace model_crispr_model_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_crispr_model");
    reader.add_event(81, 81, "end", "model_crispr_model");
    return reader;
}

#include <meta_header.hpp>
 class model_crispr_model : public prob_grad {
private:
    int n_rna_samples;
    int n_dna_samples;
    int n_ref;
    int n_alt;
    vector<vector<int> > ref_dna_counts;
    vector<vector<int> > alt_dna_counts;
    vector<vector<int> > ref_rna_counts;
    vector<vector<int> > alt_rna_counts;
    vector<double> rna_depths;
    vector<double> dna_depths;
    double dna_m_a;
    double dna_m_b;
    double dna_p_a;
    double dna_p_b;
    vector<double> rna_m_a;
    vector<double> rna_m_b;
    vector<double> rna_p_a;
    vector<double> rna_p_b;
public:
    model_crispr_model(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_crispr_model(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_crispr_model_namespace::model_crispr_model";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "n_rna_samples", "int", context__.to_vec());
            n_rna_samples = int(0);
            vals_i__ = context__.vals_i("n_rna_samples");
            pos__ = 0;
            n_rna_samples = vals_i__[pos__++];
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "n_dna_samples", "int", context__.to_vec());
            n_dna_samples = int(0);
            vals_i__ = context__.vals_i("n_dna_samples");
            pos__ = 0;
            n_dna_samples = vals_i__[pos__++];
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "n_ref", "int", context__.to_vec());
            n_ref = int(0);
            vals_i__ = context__.vals_i("n_ref");
            pos__ = 0;
            n_ref = vals_i__[pos__++];
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "n_alt", "int", context__.to_vec());
            n_alt = int(0);
            vals_i__ = context__.vals_i("n_alt");
            pos__ = 0;
            n_alt = vals_i__[pos__++];
            current_statement_begin__ = 6;
            validate_non_negative_index("ref_dna_counts", "n_ref", n_ref);
            validate_non_negative_index("ref_dna_counts", "n_dna_samples", n_dna_samples);
            context__.validate_dims("data initialization", "ref_dna_counts", "int", context__.to_vec(n_ref,n_dna_samples));
            validate_non_negative_index("ref_dna_counts", "n_ref", n_ref);
            validate_non_negative_index("ref_dna_counts", "n_dna_samples", n_dna_samples);
            ref_dna_counts = std::vector<std::vector<int> >(n_ref,std::vector<int>(n_dna_samples,int(0)));
            vals_i__ = context__.vals_i("ref_dna_counts");
            pos__ = 0;
            size_t ref_dna_counts_limit_1__ = n_dna_samples;
            for (size_t i_1__ = 0; i_1__ < ref_dna_counts_limit_1__; ++i_1__) {
                size_t ref_dna_counts_limit_0__ = n_ref;
                for (size_t i_0__ = 0; i_0__ < ref_dna_counts_limit_0__; ++i_0__) {
                    ref_dna_counts[i_0__][i_1__] = vals_i__[pos__++];
                }
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("alt_dna_counts", "n_alt", n_alt);
            validate_non_negative_index("alt_dna_counts", "n_dna_samples", n_dna_samples);
            context__.validate_dims("data initialization", "alt_dna_counts", "int", context__.to_vec(n_alt,n_dna_samples));
            validate_non_negative_index("alt_dna_counts", "n_alt", n_alt);
            validate_non_negative_index("alt_dna_counts", "n_dna_samples", n_dna_samples);
            alt_dna_counts = std::vector<std::vector<int> >(n_alt,std::vector<int>(n_dna_samples,int(0)));
            vals_i__ = context__.vals_i("alt_dna_counts");
            pos__ = 0;
            size_t alt_dna_counts_limit_1__ = n_dna_samples;
            for (size_t i_1__ = 0; i_1__ < alt_dna_counts_limit_1__; ++i_1__) {
                size_t alt_dna_counts_limit_0__ = n_alt;
                for (size_t i_0__ = 0; i_0__ < alt_dna_counts_limit_0__; ++i_0__) {
                    alt_dna_counts[i_0__][i_1__] = vals_i__[pos__++];
                }
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("ref_rna_counts", "n_ref", n_ref);
            validate_non_negative_index("ref_rna_counts", "n_rna_samples", n_rna_samples);
            context__.validate_dims("data initialization", "ref_rna_counts", "int", context__.to_vec(n_ref,n_rna_samples));
            validate_non_negative_index("ref_rna_counts", "n_ref", n_ref);
            validate_non_negative_index("ref_rna_counts", "n_rna_samples", n_rna_samples);
            ref_rna_counts = std::vector<std::vector<int> >(n_ref,std::vector<int>(n_rna_samples,int(0)));
            vals_i__ = context__.vals_i("ref_rna_counts");
            pos__ = 0;
            size_t ref_rna_counts_limit_1__ = n_rna_samples;
            for (size_t i_1__ = 0; i_1__ < ref_rna_counts_limit_1__; ++i_1__) {
                size_t ref_rna_counts_limit_0__ = n_ref;
                for (size_t i_0__ = 0; i_0__ < ref_rna_counts_limit_0__; ++i_0__) {
                    ref_rna_counts[i_0__][i_1__] = vals_i__[pos__++];
                }
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("alt_rna_counts", "n_alt", n_alt);
            validate_non_negative_index("alt_rna_counts", "n_rna_samples", n_rna_samples);
            context__.validate_dims("data initialization", "alt_rna_counts", "int", context__.to_vec(n_alt,n_rna_samples));
            validate_non_negative_index("alt_rna_counts", "n_alt", n_alt);
            validate_non_negative_index("alt_rna_counts", "n_rna_samples", n_rna_samples);
            alt_rna_counts = std::vector<std::vector<int> >(n_alt,std::vector<int>(n_rna_samples,int(0)));
            vals_i__ = context__.vals_i("alt_rna_counts");
            pos__ = 0;
            size_t alt_rna_counts_limit_1__ = n_rna_samples;
            for (size_t i_1__ = 0; i_1__ < alt_rna_counts_limit_1__; ++i_1__) {
                size_t alt_rna_counts_limit_0__ = n_alt;
                for (size_t i_0__ = 0; i_0__ < alt_rna_counts_limit_0__; ++i_0__) {
                    alt_rna_counts[i_0__][i_1__] = vals_i__[pos__++];
                }
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("rna_depths", "n_rna_samples", n_rna_samples);
            context__.validate_dims("data initialization", "rna_depths", "double", context__.to_vec(n_rna_samples));
            validate_non_negative_index("rna_depths", "n_rna_samples", n_rna_samples);
            rna_depths = std::vector<double>(n_rna_samples,double(0));
            vals_r__ = context__.vals_r("rna_depths");
            pos__ = 0;
            size_t rna_depths_limit_0__ = n_rna_samples;
            for (size_t i_0__ = 0; i_0__ < rna_depths_limit_0__; ++i_0__) {
                rna_depths[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 12;
            validate_non_negative_index("dna_depths", "n_dna_samples", n_dna_samples);
            context__.validate_dims("data initialization", "dna_depths", "double", context__.to_vec(n_dna_samples));
            validate_non_negative_index("dna_depths", "n_dna_samples", n_dna_samples);
            dna_depths = std::vector<double>(n_dna_samples,double(0));
            vals_r__ = context__.vals_r("dna_depths");
            pos__ = 0;
            size_t dna_depths_limit_0__ = n_dna_samples;
            for (size_t i_0__ = 0; i_0__ < dna_depths_limit_0__; ++i_0__) {
                dna_depths[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 14;
            context__.validate_dims("data initialization", "dna_m_a", "double", context__.to_vec());
            dna_m_a = double(0);
            vals_r__ = context__.vals_r("dna_m_a");
            pos__ = 0;
            dna_m_a = vals_r__[pos__++];
            current_statement_begin__ = 15;
            context__.validate_dims("data initialization", "dna_m_b", "double", context__.to_vec());
            dna_m_b = double(0);
            vals_r__ = context__.vals_r("dna_m_b");
            pos__ = 0;
            dna_m_b = vals_r__[pos__++];
            current_statement_begin__ = 17;
            context__.validate_dims("data initialization", "dna_p_a", "double", context__.to_vec());
            dna_p_a = double(0);
            vals_r__ = context__.vals_r("dna_p_a");
            pos__ = 0;
            dna_p_a = vals_r__[pos__++];
            current_statement_begin__ = 18;
            context__.validate_dims("data initialization", "dna_p_b", "double", context__.to_vec());
            dna_p_b = double(0);
            vals_r__ = context__.vals_r("dna_p_b");
            pos__ = 0;
            dna_p_b = vals_r__[pos__++];
            current_statement_begin__ = 20;
            validate_non_negative_index("rna_m_a", "2", 2);
            context__.validate_dims("data initialization", "rna_m_a", "double", context__.to_vec(2));
            validate_non_negative_index("rna_m_a", "2", 2);
            rna_m_a = std::vector<double>(2,double(0));
            vals_r__ = context__.vals_r("rna_m_a");
            pos__ = 0;
            size_t rna_m_a_limit_0__ = 2;
            for (size_t i_0__ = 0; i_0__ < rna_m_a_limit_0__; ++i_0__) {
                rna_m_a[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 21;
            validate_non_negative_index("rna_m_b", "2", 2);
            context__.validate_dims("data initialization", "rna_m_b", "double", context__.to_vec(2));
            validate_non_negative_index("rna_m_b", "2", 2);
            rna_m_b = std::vector<double>(2,double(0));
            vals_r__ = context__.vals_r("rna_m_b");
            pos__ = 0;
            size_t rna_m_b_limit_0__ = 2;
            for (size_t i_0__ = 0; i_0__ < rna_m_b_limit_0__; ++i_0__) {
                rna_m_b[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 22;
            validate_non_negative_index("rna_p_a", "2", 2);
            context__.validate_dims("data initialization", "rna_p_a", "double", context__.to_vec(2));
            validate_non_negative_index("rna_p_a", "2", 2);
            rna_p_a = std::vector<double>(2,double(0));
            vals_r__ = context__.vals_r("rna_p_a");
            pos__ = 0;
            size_t rna_p_a_limit_0__ = 2;
            for (size_t i_0__ = 0; i_0__ < rna_p_a_limit_0__; ++i_0__) {
                rna_p_a[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 23;
            validate_non_negative_index("rna_p_b", "2", 2);
            context__.validate_dims("data initialization", "rna_p_b", "double", context__.to_vec(2));
            validate_non_negative_index("rna_p_b", "2", 2);
            rna_p_b = std::vector<double>(2,double(0));
            vals_r__ = context__.vals_r("rna_p_b");
            pos__ = 0;
            size_t rna_p_b_limit_0__ = 2;
            for (size_t i_0__ = 0; i_0__ < rna_p_b_limit_0__; ++i_0__) {
                rna_p_b[i_0__] = vals_r__[pos__++];
            }

            // validate, data variables
            current_statement_begin__ = 2;
            check_greater_or_equal(function__,"n_rna_samples",n_rna_samples,0);
            current_statement_begin__ = 3;
            check_greater_or_equal(function__,"n_dna_samples",n_dna_samples,0);
            current_statement_begin__ = 4;
            check_greater_or_equal(function__,"n_ref",n_ref,1);
            current_statement_begin__ = 5;
            check_greater_or_equal(function__,"n_alt",n_alt,1);
            current_statement_begin__ = 6;
            for (int k0__ = 0; k0__ < n_ref; ++k0__) {
                for (int k1__ = 0; k1__ < n_dna_samples; ++k1__) {
                    check_greater_or_equal(function__,"ref_dna_counts[k0__][k1__]",ref_dna_counts[k0__][k1__],0);
                }
            }
            current_statement_begin__ = 7;
            for (int k0__ = 0; k0__ < n_alt; ++k0__) {
                for (int k1__ = 0; k1__ < n_dna_samples; ++k1__) {
                    check_greater_or_equal(function__,"alt_dna_counts[k0__][k1__]",alt_dna_counts[k0__][k1__],0);
                }
            }
            current_statement_begin__ = 8;
            for (int k0__ = 0; k0__ < n_ref; ++k0__) {
                for (int k1__ = 0; k1__ < n_rna_samples; ++k1__) {
                    check_greater_or_equal(function__,"ref_rna_counts[k0__][k1__]",ref_rna_counts[k0__][k1__],0);
                }
            }
            current_statement_begin__ = 9;
            for (int k0__ = 0; k0__ < n_alt; ++k0__) {
                for (int k1__ = 0; k1__ < n_rna_samples; ++k1__) {
                    check_greater_or_equal(function__,"alt_rna_counts[k0__][k1__]",alt_rna_counts[k0__][k1__],0);
                }
            }
            current_statement_begin__ = 11;
            for (int k0__ = 0; k0__ < n_rna_samples; ++k0__) {
                check_greater_or_equal(function__,"rna_depths[k0__]",rna_depths[k0__],0);
            }
            current_statement_begin__ = 12;
            for (int k0__ = 0; k0__ < n_dna_samples; ++k0__) {
                check_greater_or_equal(function__,"dna_depths[k0__]",dna_depths[k0__],0);
            }
            current_statement_begin__ = 14;
            check_greater_or_equal(function__,"dna_m_a",dna_m_a,0);
            current_statement_begin__ = 15;
            check_greater_or_equal(function__,"dna_m_b",dna_m_b,0);
            current_statement_begin__ = 17;
            check_greater_or_equal(function__,"dna_p_a",dna_p_a,0);
            current_statement_begin__ = 18;
            check_greater_or_equal(function__,"dna_p_b",dna_p_b,0);
            current_statement_begin__ = 20;
            for (int k0__ = 0; k0__ < 2; ++k0__) {
                check_greater_or_equal(function__,"rna_m_a[k0__]",rna_m_a[k0__],0);
            }
            current_statement_begin__ = 21;
            for (int k0__ = 0; k0__ < 2; ++k0__) {
                check_greater_or_equal(function__,"rna_m_b[k0__]",rna_m_b[k0__],0);
            }
            current_statement_begin__ = 22;
            for (int k0__ = 0; k0__ < 2; ++k0__) {
                check_greater_or_equal(function__,"rna_p_a[k0__]",rna_p_a[k0__],0);
            }
            current_statement_begin__ = 23;
            for (int k0__ = 0; k0__ < 2; ++k0__) {
                check_greater_or_equal(function__,"rna_p_b[k0__]",rna_p_b[k0__],0);
            }
            // initialize data variables


            // validate transformed data

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 26;
            validate_non_negative_index("dna_m_ref", "n_ref", n_ref);
            num_params_r__ += n_ref;
            current_statement_begin__ = 27;
            validate_non_negative_index("dna_m_alt", "n_alt", n_alt);
            num_params_r__ += n_alt;
            current_statement_begin__ = 28;
            ++num_params_r__;
            current_statement_begin__ = 30;
            validate_non_negative_index("rna_m", "2", 2);
            num_params_r__ += 2;
            current_statement_begin__ = 31;
            validate_non_negative_index("rna_p", "2", 2);
            num_params_r__ += 2;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_crispr_model() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("dna_m_ref")))
            throw std::runtime_error("variable dna_m_ref missing");
        vals_r__ = context__.vals_r("dna_m_ref");
        pos__ = 0U;
        validate_non_negative_index("dna_m_ref", "n_ref", n_ref);
        context__.validate_dims("initialization", "dna_m_ref", "vector_d", context__.to_vec(n_ref));
        vector_d dna_m_ref(static_cast<Eigen::VectorXd::Index>(n_ref));
        for (int j1__ = 0U; j1__ < n_ref; ++j1__)
            dna_m_ref(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,dna_m_ref);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable dna_m_ref: ") + e.what());
        }

        if (!(context__.contains_r("dna_m_alt")))
            throw std::runtime_error("variable dna_m_alt missing");
        vals_r__ = context__.vals_r("dna_m_alt");
        pos__ = 0U;
        validate_non_negative_index("dna_m_alt", "n_alt", n_alt);
        context__.validate_dims("initialization", "dna_m_alt", "vector_d", context__.to_vec(n_alt));
        vector_d dna_m_alt(static_cast<Eigen::VectorXd::Index>(n_alt));
        for (int j1__ = 0U; j1__ < n_alt; ++j1__)
            dna_m_alt(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,dna_m_alt);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable dna_m_alt: ") + e.what());
        }

        if (!(context__.contains_r("dna_p")))
            throw std::runtime_error("variable dna_p missing");
        vals_r__ = context__.vals_r("dna_p");
        pos__ = 0U;
        context__.validate_dims("initialization", "dna_p", "double", context__.to_vec());
        double dna_p(0);
        dna_p = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,dna_p);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable dna_p: ") + e.what());
        }

        if (!(context__.contains_r("rna_m")))
            throw std::runtime_error("variable rna_m missing");
        vals_r__ = context__.vals_r("rna_m");
        pos__ = 0U;
        validate_non_negative_index("rna_m", "2", 2);
        context__.validate_dims("initialization", "rna_m", "vector_d", context__.to_vec(2));
        vector_d rna_m(static_cast<Eigen::VectorXd::Index>(2));
        for (int j1__ = 0U; j1__ < 2; ++j1__)
            rna_m(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,rna_m);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable rna_m: ") + e.what());
        }

        if (!(context__.contains_r("rna_p")))
            throw std::runtime_error("variable rna_p missing");
        vals_r__ = context__.vals_r("rna_p");
        pos__ = 0U;
        validate_non_negative_index("rna_p", "2", 2);
        context__.validate_dims("initialization", "rna_p", "vector_d", context__.to_vec(2));
        vector_d rna_p(static_cast<Eigen::VectorXd::Index>(2));
        for (int j1__ = 0U; j1__ < 2; ++j1__)
            rna_p(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,rna_p);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable rna_p: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<T__> in__(params_r__,params_i__);

            Eigen::Matrix<T__,Eigen::Dynamic,1>  dna_m_ref;
            (void) dna_m_ref;  // dummy to suppress unused var warning
            if (jacobian__)
                dna_m_ref = in__.vector_lb_constrain(0,n_ref,lp__);
            else
                dna_m_ref = in__.vector_lb_constrain(0,n_ref);

            Eigen::Matrix<T__,Eigen::Dynamic,1>  dna_m_alt;
            (void) dna_m_alt;  // dummy to suppress unused var warning
            if (jacobian__)
                dna_m_alt = in__.vector_lb_constrain(0,n_alt,lp__);
            else
                dna_m_alt = in__.vector_lb_constrain(0,n_alt);

            T__ dna_p;
            (void) dna_p;  // dummy to suppress unused var warning
            if (jacobian__)
                dna_p = in__.scalar_lb_constrain(0,lp__);
            else
                dna_p = in__.scalar_lb_constrain(0);

            Eigen::Matrix<T__,Eigen::Dynamic,1>  rna_m;
            (void) rna_m;  // dummy to suppress unused var warning
            if (jacobian__)
                rna_m = in__.vector_lb_constrain(0,2,lp__);
            else
                rna_m = in__.vector_lb_constrain(0,2);

            Eigen::Matrix<T__,Eigen::Dynamic,1>  rna_p;
            (void) rna_p;  // dummy to suppress unused var warning
            if (jacobian__)
                rna_p = in__.vector_lb_constrain(0,2,lp__);
            else
                rna_p = in__.vector_lb_constrain(0,2);


            // transformed parameters



            // validate transformed parameters

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            // model body

            current_statement_begin__ = 35;
            for (int bc = 1; bc <= n_ref; ++bc) {

                current_statement_begin__ = 36;
                lp_accum__.add(gamma_log<propto__>(get_base1(dna_m_ref,bc,"dna_m_ref",1), dna_m_a, dna_m_b));
            }
            current_statement_begin__ = 39;
            for (int bc = 1; bc <= n_alt; ++bc) {

                current_statement_begin__ = 40;
                lp_accum__.add(gamma_log<propto__>(get_base1(dna_m_alt,bc,"dna_m_alt",1), dna_m_a, dna_m_b));
            }
            current_statement_begin__ = 43;
            lp_accum__.add(gamma_log<propto__>(dna_p, dna_p_a, dna_p_b));
            current_statement_begin__ = 48;
            for (int s = 1; s <= n_dna_samples; ++s) {

                current_statement_begin__ = 49;
                for (int bc = 1; bc <= n_ref; ++bc) {

                    current_statement_begin__ = 50;
                    lp_accum__.add(neg_binomial_2_log<propto__>(get_base1(get_base1(ref_dna_counts,bc,"ref_dna_counts",1),s,"ref_dna_counts",2), (get_base1(dna_m_ref,bc,"dna_m_ref",1) * get_base1(dna_depths,s,"dna_depths",1)), dna_p));
                }
                current_statement_begin__ = 53;
                for (int bc = 1; bc <= n_alt; ++bc) {

                    current_statement_begin__ = 54;
                    lp_accum__.add(neg_binomial_2_log<propto__>(get_base1(get_base1(alt_dna_counts,bc,"alt_dna_counts",1),s,"alt_dna_counts",2), (get_base1(dna_m_alt,bc,"dna_m_alt",1) * get_base1(dna_depths,s,"dna_depths",1)), dna_p));
                }
            }
            current_statement_begin__ = 58;
            for (int allele = 1; allele <= 2; ++allele) {

                current_statement_begin__ = 59;
                lp_accum__.add(gamma_log<propto__>(get_base1(rna_m,allele,"rna_m",1), get_base1(rna_m_a,allele,"rna_m_a",1), get_base1(rna_m_b,allele,"rna_m_b",1)));
                current_statement_begin__ = 60;
                lp_accum__.add(gamma_log<propto__>(get_base1(rna_p,allele,"rna_p",1), get_base1(rna_p_a,allele,"rna_p_a",1), get_base1(rna_p_b,allele,"rna_p_b",1)));
            }
            current_statement_begin__ = 63;
            for (int s = 1; s <= n_rna_samples; ++s) {

                current_statement_begin__ = 64;
                for (int bc = 1; bc <= n_ref; ++bc) {

                    current_statement_begin__ = 65;
                    lp_accum__.add(neg_binomial_2_log<propto__>(get_base1(get_base1(ref_rna_counts,bc,"ref_rna_counts",1),s,"ref_rna_counts",2), ((get_base1(rna_m,1,"rna_m",1) * get_base1(rna_depths,s,"rna_depths",1)) * get_base1(dna_m_ref,bc,"dna_m_ref",1)), get_base1(rna_p,1,"rna_p",1)));
                }
                current_statement_begin__ = 68;
                for (int bc = 1; bc <= n_alt; ++bc) {

                    current_statement_begin__ = 69;
                    lp_accum__.add(neg_binomial_2_log<propto__>(get_base1(get_base1(alt_rna_counts,bc,"alt_rna_counts",1),s,"alt_rna_counts",2), ((get_base1(rna_m,2,"rna_m",1) * get_base1(rna_depths,s,"rna_depths",1)) * get_base1(dna_m_alt,bc,"dna_m_alt",1)), get_base1(rna_p,2,"rna_p",1)));
                }
            }

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("dna_m_ref");
        names__.push_back("dna_m_alt");
        names__.push_back("dna_p");
        names__.push_back("rna_m");
        names__.push_back("rna_p");
        names__.push_back("ref_act");
        names__.push_back("alt_act");
        names__.push_back("transcription_shift");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(n_ref);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_alt);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(2);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(2);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        vars__.resize(0);
        stan::io::reader<double> in__(params_r__,params_i__);
        static const char* function__ = "model_crispr_model_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector_d dna_m_ref = in__.vector_lb_constrain(0,n_ref);
        vector_d dna_m_alt = in__.vector_lb_constrain(0,n_alt);
        double dna_p = in__.scalar_lb_constrain(0);
        vector_d rna_m = in__.vector_lb_constrain(0,2);
        vector_d rna_p = in__.vector_lb_constrain(0,2);
            for (int k_0__ = 0; k_0__ < n_ref; ++k_0__) {
            vars__.push_back(dna_m_ref[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < n_alt; ++k_0__) {
            vars__.push_back(dna_m_alt[k_0__]);
            }
        vars__.push_back(dna_p);
            for (int k_0__ = 0; k_0__ < 2; ++k_0__) {
            vars__.push_back(rna_m[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < 2; ++k_0__) {
            vars__.push_back(rna_p[k_0__]);
            }

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {



            // validate transformed parameters

            // write transformed parameters

            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 75;
            double ref_act(0.0);
            (void) ref_act;  // dummy to suppress unused var warning

            stan::math::initialize(ref_act, std::numeric_limits<double>::quiet_NaN());
            stan::math::fill(ref_act,DUMMY_VAR__);
            current_statement_begin__ = 76;
            double alt_act(0.0);
            (void) alt_act;  // dummy to suppress unused var warning

            stan::math::initialize(alt_act, std::numeric_limits<double>::quiet_NaN());
            stan::math::fill(alt_act,DUMMY_VAR__);
            current_statement_begin__ = 77;
            double transcription_shift(0.0);
            (void) transcription_shift;  // dummy to suppress unused var warning

            stan::math::initialize(transcription_shift, std::numeric_limits<double>::quiet_NaN());
            stan::math::fill(transcription_shift,DUMMY_VAR__);


            current_statement_begin__ = 78;
            stan::math::assign(ref_act, log(get_base1(rna_m,1,"rna_m",1)));
            current_statement_begin__ = 79;
            stan::math::assign(alt_act, log(get_base1(rna_m,2,"rna_m",1)));
            current_statement_begin__ = 80;
            stan::math::assign(transcription_shift, (alt_act - ref_act));

            // validate generated quantities
            current_statement_begin__ = 75;
            current_statement_begin__ = 76;
            current_statement_begin__ = 77;

            // write generated quantities
        vars__.push_back(ref_act);
        vars__.push_back(alt_act);
        vars__.push_back(transcription_shift);

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_crispr_model";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= n_ref; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "dna_m_ref" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_alt; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "dna_m_alt" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "dna_p";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rna_m" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rna_p" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "ref_act";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alt_act";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "transcription_shift";
        param_names__.push_back(param_name_stream__.str());
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= n_ref; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "dna_m_ref" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= n_alt; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "dna_m_alt" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "dna_p";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rna_m" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rna_p" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "ref_act";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alt_act";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "transcription_shift";
        param_names__.push_back(param_name_stream__.str());
    }

}; // model

}

typedef model_crispr_model_namespace::model_crispr_model stan_model;


#endif
