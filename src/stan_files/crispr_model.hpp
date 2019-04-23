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
// Code generated by Stan version 2.18.1

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

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_crispr_model");
    reader.add_event(69, 67, "end", "model_crispr_model");
    return reader;
}

#include <meta_header.hpp>
 class model_crispr_model : public prob_grad {
private:
    int n_gRNA;
    vector<int> t0_counts;
    vector<int> case_counts;
    vector<int> ctrl_counts;
    double t0_depths;
    double case_depths;
    double ctrl_depths;
    double t0_m_a;
    double t0_m_b;
    double t0_p_a;
    double t0_p_b;
    vector<double> output_m_a;
    vector<double> output_m_b;
    vector<double> output_p_a;
    vector<double> output_p_b;
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
        typedef double local_scalar_t__;

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
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "n_gRNA", "int", context__.to_vec());
            n_gRNA = int(0);
            vals_i__ = context__.vals_i("n_gRNA");
            pos__ = 0;
            n_gRNA = vals_i__[pos__++];
            current_statement_begin__ = 4;
            validate_non_negative_index("t0_counts", "n_gRNA", n_gRNA);
            context__.validate_dims("data initialization", "t0_counts", "int", context__.to_vec(n_gRNA));
            validate_non_negative_index("t0_counts", "n_gRNA", n_gRNA);
            t0_counts = std::vector<int>(n_gRNA,int(0));
            vals_i__ = context__.vals_i("t0_counts");
            pos__ = 0;
            size_t t0_counts_limit_0__ = n_gRNA;
            for (size_t i_0__ = 0; i_0__ < t0_counts_limit_0__; ++i_0__) {
                t0_counts[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 5;
            validate_non_negative_index("case_counts", "n_gRNA", n_gRNA);
            context__.validate_dims("data initialization", "case_counts", "int", context__.to_vec(n_gRNA));
            validate_non_negative_index("case_counts", "n_gRNA", n_gRNA);
            case_counts = std::vector<int>(n_gRNA,int(0));
            vals_i__ = context__.vals_i("case_counts");
            pos__ = 0;
            size_t case_counts_limit_0__ = n_gRNA;
            for (size_t i_0__ = 0; i_0__ < case_counts_limit_0__; ++i_0__) {
                case_counts[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 6;
            validate_non_negative_index("ctrl_counts", "n_gRNA", n_gRNA);
            context__.validate_dims("data initialization", "ctrl_counts", "int", context__.to_vec(n_gRNA));
            validate_non_negative_index("ctrl_counts", "n_gRNA", n_gRNA);
            ctrl_counts = std::vector<int>(n_gRNA,int(0));
            vals_i__ = context__.vals_i("ctrl_counts");
            pos__ = 0;
            size_t ctrl_counts_limit_0__ = n_gRNA;
            for (size_t i_0__ = 0; i_0__ < ctrl_counts_limit_0__; ++i_0__) {
                ctrl_counts[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 8;
            context__.validate_dims("data initialization", "t0_depths", "double", context__.to_vec());
            t0_depths = double(0);
            vals_r__ = context__.vals_r("t0_depths");
            pos__ = 0;
            t0_depths = vals_r__[pos__++];
            current_statement_begin__ = 9;
            context__.validate_dims("data initialization", "case_depths", "double", context__.to_vec());
            case_depths = double(0);
            vals_r__ = context__.vals_r("case_depths");
            pos__ = 0;
            case_depths = vals_r__[pos__++];
            current_statement_begin__ = 10;
            context__.validate_dims("data initialization", "ctrl_depths", "double", context__.to_vec());
            ctrl_depths = double(0);
            vals_r__ = context__.vals_r("ctrl_depths");
            pos__ = 0;
            ctrl_depths = vals_r__[pos__++];
            current_statement_begin__ = 12;
            context__.validate_dims("data initialization", "t0_m_a", "double", context__.to_vec());
            t0_m_a = double(0);
            vals_r__ = context__.vals_r("t0_m_a");
            pos__ = 0;
            t0_m_a = vals_r__[pos__++];
            current_statement_begin__ = 13;
            context__.validate_dims("data initialization", "t0_m_b", "double", context__.to_vec());
            t0_m_b = double(0);
            vals_r__ = context__.vals_r("t0_m_b");
            pos__ = 0;
            t0_m_b = vals_r__[pos__++];
            current_statement_begin__ = 15;
            context__.validate_dims("data initialization", "t0_p_a", "double", context__.to_vec());
            t0_p_a = double(0);
            vals_r__ = context__.vals_r("t0_p_a");
            pos__ = 0;
            t0_p_a = vals_r__[pos__++];
            current_statement_begin__ = 16;
            context__.validate_dims("data initialization", "t0_p_b", "double", context__.to_vec());
            t0_p_b = double(0);
            vals_r__ = context__.vals_r("t0_p_b");
            pos__ = 0;
            t0_p_b = vals_r__[pos__++];
            current_statement_begin__ = 18;
            validate_non_negative_index("output_m_a", "2", 2);
            context__.validate_dims("data initialization", "output_m_a", "double", context__.to_vec(2));
            validate_non_negative_index("output_m_a", "2", 2);
            output_m_a = std::vector<double>(2,double(0));
            vals_r__ = context__.vals_r("output_m_a");
            pos__ = 0;
            size_t output_m_a_limit_0__ = 2;
            for (size_t i_0__ = 0; i_0__ < output_m_a_limit_0__; ++i_0__) {
                output_m_a[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 19;
            validate_non_negative_index("output_m_b", "2", 2);
            context__.validate_dims("data initialization", "output_m_b", "double", context__.to_vec(2));
            validate_non_negative_index("output_m_b", "2", 2);
            output_m_b = std::vector<double>(2,double(0));
            vals_r__ = context__.vals_r("output_m_b");
            pos__ = 0;
            size_t output_m_b_limit_0__ = 2;
            for (size_t i_0__ = 0; i_0__ < output_m_b_limit_0__; ++i_0__) {
                output_m_b[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 20;
            validate_non_negative_index("output_p_a", "2", 2);
            context__.validate_dims("data initialization", "output_p_a", "double", context__.to_vec(2));
            validate_non_negative_index("output_p_a", "2", 2);
            output_p_a = std::vector<double>(2,double(0));
            vals_r__ = context__.vals_r("output_p_a");
            pos__ = 0;
            size_t output_p_a_limit_0__ = 2;
            for (size_t i_0__ = 0; i_0__ < output_p_a_limit_0__; ++i_0__) {
                output_p_a[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 21;
            validate_non_negative_index("output_p_b", "2", 2);
            context__.validate_dims("data initialization", "output_p_b", "double", context__.to_vec(2));
            validate_non_negative_index("output_p_b", "2", 2);
            output_p_b = std::vector<double>(2,double(0));
            vals_r__ = context__.vals_r("output_p_b");
            pos__ = 0;
            size_t output_p_b_limit_0__ = 2;
            for (size_t i_0__ = 0; i_0__ < output_p_b_limit_0__; ++i_0__) {
                output_p_b[i_0__] = vals_r__[pos__++];
            }

            // validate, data variables
            current_statement_begin__ = 2;
            check_greater_or_equal(function__,"n_gRNA",n_gRNA,1);
            current_statement_begin__ = 4;
            for (int k0__ = 0; k0__ < n_gRNA; ++k0__) {
                check_greater_or_equal(function__,"t0_counts[k0__]",t0_counts[k0__],0);
            }
            current_statement_begin__ = 5;
            for (int k0__ = 0; k0__ < n_gRNA; ++k0__) {
                check_greater_or_equal(function__,"case_counts[k0__]",case_counts[k0__],0);
            }
            current_statement_begin__ = 6;
            for (int k0__ = 0; k0__ < n_gRNA; ++k0__) {
                check_greater_or_equal(function__,"ctrl_counts[k0__]",ctrl_counts[k0__],0);
            }
            current_statement_begin__ = 8;
            check_greater_or_equal(function__,"t0_depths",t0_depths,0);
            current_statement_begin__ = 9;
            check_greater_or_equal(function__,"case_depths",case_depths,0);
            current_statement_begin__ = 10;
            check_greater_or_equal(function__,"ctrl_depths",ctrl_depths,0);
            current_statement_begin__ = 12;
            check_greater_or_equal(function__,"t0_m_a",t0_m_a,0);
            current_statement_begin__ = 13;
            check_greater_or_equal(function__,"t0_m_b",t0_m_b,0);
            current_statement_begin__ = 15;
            check_greater_or_equal(function__,"t0_p_a",t0_p_a,0);
            current_statement_begin__ = 16;
            check_greater_or_equal(function__,"t0_p_b",t0_p_b,0);
            current_statement_begin__ = 18;
            for (int k0__ = 0; k0__ < 2; ++k0__) {
                check_greater_or_equal(function__,"output_m_a[k0__]",output_m_a[k0__],0);
            }
            current_statement_begin__ = 19;
            for (int k0__ = 0; k0__ < 2; ++k0__) {
                check_greater_or_equal(function__,"output_m_b[k0__]",output_m_b[k0__],0);
            }
            current_statement_begin__ = 20;
            for (int k0__ = 0; k0__ < 2; ++k0__) {
                check_greater_or_equal(function__,"output_p_a[k0__]",output_p_a[k0__],0);
            }
            current_statement_begin__ = 21;
            for (int k0__ = 0; k0__ < 2; ++k0__) {
                check_greater_or_equal(function__,"output_p_b[k0__]",output_p_b[k0__],0);
            }
            // initialize data variables


            // validate transformed data

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 24;
            validate_non_negative_index("t0_m", "n_gRNA", n_gRNA);
            num_params_r__ += n_gRNA;
            current_statement_begin__ = 25;
            ++num_params_r__;
            current_statement_begin__ = 27;
            validate_non_negative_index("output_m", "2", 2);
            num_params_r__ += 2;
            current_statement_begin__ = 28;
            validate_non_negative_index("output_p", "2", 2);
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

        if (!(context__.contains_r("t0_m")))
            throw std::runtime_error("variable t0_m missing");
        vals_r__ = context__.vals_r("t0_m");
        pos__ = 0U;
        validate_non_negative_index("t0_m", "n_gRNA", n_gRNA);
        context__.validate_dims("initialization", "t0_m", "vector_d", context__.to_vec(n_gRNA));
        vector_d t0_m(static_cast<Eigen::VectorXd::Index>(n_gRNA));
        for (int j1__ = 0U; j1__ < n_gRNA; ++j1__)
            t0_m(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,t0_m);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable t0_m: ") + e.what());
        }

        if (!(context__.contains_r("t0_p")))
            throw std::runtime_error("variable t0_p missing");
        vals_r__ = context__.vals_r("t0_p");
        pos__ = 0U;
        context__.validate_dims("initialization", "t0_p", "double", context__.to_vec());
        double t0_p(0);
        t0_p = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,t0_p);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable t0_p: ") + e.what());
        }

        if (!(context__.contains_r("output_m")))
            throw std::runtime_error("variable output_m missing");
        vals_r__ = context__.vals_r("output_m");
        pos__ = 0U;
        validate_non_negative_index("output_m", "2", 2);
        context__.validate_dims("initialization", "output_m", "vector_d", context__.to_vec(2));
        vector_d output_m(static_cast<Eigen::VectorXd::Index>(2));
        for (int j1__ = 0U; j1__ < 2; ++j1__)
            output_m(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,output_m);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable output_m: ") + e.what());
        }

        if (!(context__.contains_r("output_p")))
            throw std::runtime_error("variable output_p missing");
        vals_r__ = context__.vals_r("output_p");
        pos__ = 0U;
        validate_non_negative_index("output_p", "2", 2);
        context__.validate_dims("initialization", "output_p", "vector_d", context__.to_vec(2));
        vector_d output_p(static_cast<Eigen::VectorXd::Index>(2));
        for (int j1__ = 0U; j1__ < 2; ++j1__)
            output_p(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,output_p);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable output_p: ") + e.what());
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

        typedef T__ local_scalar_t__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);

            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  t0_m;
            (void) t0_m;  // dummy to suppress unused var warning
            if (jacobian__)
                t0_m = in__.vector_lb_constrain(0,n_gRNA,lp__);
            else
                t0_m = in__.vector_lb_constrain(0,n_gRNA);

            local_scalar_t__ t0_p;
            (void) t0_p;  // dummy to suppress unused var warning
            if (jacobian__)
                t0_p = in__.scalar_lb_constrain(0,lp__);
            else
                t0_p = in__.scalar_lb_constrain(0);

            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  output_m;
            (void) output_m;  // dummy to suppress unused var warning
            if (jacobian__)
                output_m = in__.vector_lb_constrain(0,2,lp__);
            else
                output_m = in__.vector_lb_constrain(0,2);

            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  output_p;
            (void) output_p;  // dummy to suppress unused var warning
            if (jacobian__)
                output_p = in__.vector_lb_constrain(0,2,lp__);
            else
                output_p = in__.vector_lb_constrain(0,2);


            // transformed parameters



            // validate transformed parameters

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            // model body

            current_statement_begin__ = 33;
            for (int g = 1; g <= n_gRNA; ++g) {

                current_statement_begin__ = 34;
                lp_accum__.add(gamma_log<propto__>(get_base1(t0_m,g,"t0_m",1), t0_m_a, t0_m_b));
            }
            current_statement_begin__ = 37;
            lp_accum__.add(gamma_log<propto__>(t0_p, t0_p_a, t0_p_b));
            current_statement_begin__ = 42;
            for (int g = 1; g <= n_gRNA; ++g) {

                current_statement_begin__ = 43;
                lp_accum__.add(neg_binomial_2_log<propto__>(get_base1(t0_counts,g,"t0_counts",1), (get_base1(t0_m,g,"t0_m",1) * t0_depths), t0_p));
            }
            current_statement_begin__ = 47;
            for (int group = 1; group <= 2; ++group) {

                current_statement_begin__ = 48;
                lp_accum__.add(gamma_log<propto__>(get_base1(output_m,group,"output_m",1), get_base1(output_m_a,group,"output_m_a",1), get_base1(output_m_b,group,"output_m_b",1)));
                current_statement_begin__ = 49;
                lp_accum__.add(gamma_log<propto__>(get_base1(output_p,group,"output_p",1), get_base1(output_p_a,group,"output_p_a",1), get_base1(output_p_b,group,"output_p_b",1)));
            }
            current_statement_begin__ = 52;
            for (int g = 1; g <= n_gRNA; ++g) {

                current_statement_begin__ = 53;
                lp_accum__.add(neg_binomial_2_log<propto__>(get_base1(ctrl_counts,g,"ctrl_counts",1), ((get_base1(output_m,1,"output_m",1) * ctrl_depths) * get_base1(t0_m,g,"t0_m",1)), get_base1(output_p,1,"output_p",1)));
            }
            current_statement_begin__ = 56;
            for (int g = 1; g <= n_gRNA; ++g) {

                current_statement_begin__ = 57;
                lp_accum__.add(neg_binomial_2_log<propto__>(get_base1(case_counts,g,"case_counts",1), ((get_base1(output_m,2,"output_m",1) * case_depths) * get_base1(t0_m,g,"t0_m",1)), get_base1(output_p,2,"output_p",1)));
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
        names__.push_back("t0_m");
        names__.push_back("t0_p");
        names__.push_back("output_m");
        names__.push_back("output_p");
        names__.push_back("ctrl_act");
        names__.push_back("case_act");
        names__.push_back("condition_shift");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(n_gRNA);
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
        typedef double local_scalar_t__;

        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "model_crispr_model_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector_d t0_m = in__.vector_lb_constrain(0,n_gRNA);
        double t0_p = in__.scalar_lb_constrain(0);
        vector_d output_m = in__.vector_lb_constrain(0,2);
        vector_d output_p = in__.vector_lb_constrain(0,2);
            for (int k_0__ = 0; k_0__ < n_gRNA; ++k_0__) {
            vars__.push_back(t0_m[k_0__]);
            }
        vars__.push_back(t0_p);
            for (int k_0__ = 0; k_0__ < 2; ++k_0__) {
            vars__.push_back(output_m[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < 2; ++k_0__) {
            vars__.push_back(output_p[k_0__]);
            }

        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {



            // validate transformed parameters

            // write transformed parameters
            if (include_tparams__) {
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 61;
            local_scalar_t__ ctrl_act;
            (void) ctrl_act;  // dummy to suppress unused var warning

            stan::math::initialize(ctrl_act, DUMMY_VAR__);
            stan::math::fill(ctrl_act,DUMMY_VAR__);
            current_statement_begin__ = 62;
            local_scalar_t__ case_act;
            (void) case_act;  // dummy to suppress unused var warning

            stan::math::initialize(case_act, DUMMY_VAR__);
            stan::math::fill(case_act,DUMMY_VAR__);
            current_statement_begin__ = 63;
            local_scalar_t__ condition_shift;
            (void) condition_shift;  // dummy to suppress unused var warning

            stan::math::initialize(condition_shift, DUMMY_VAR__);
            stan::math::fill(condition_shift,DUMMY_VAR__);


            current_statement_begin__ = 64;
            stan::math::assign(ctrl_act, stan::math::log(get_base1(output_m,1,"output_m",1)));
            current_statement_begin__ = 65;
            stan::math::assign(case_act, stan::math::log(get_base1(output_m,2,"output_m",1)));
            current_statement_begin__ = 66;
            stan::math::assign(condition_shift, (case_act - ctrl_act));

            // validate generated quantities
            current_statement_begin__ = 61;
            current_statement_begin__ = 62;
            current_statement_begin__ = 63;

            // write generated quantities
        vars__.push_back(ctrl_act);
        vars__.push_back(case_act);
        vars__.push_back(condition_shift);

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
        for (int k_0__ = 1; k_0__ <= n_gRNA; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "t0_m" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "t0_p";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "output_m" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "output_p" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
        }


        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "ctrl_act";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "case_act";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "condition_shift";
        param_names__.push_back(param_name_stream__.str());
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= n_gRNA; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "t0_m" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "t0_p";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "output_m" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "output_p" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
        }


        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "ctrl_act";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "case_act";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "condition_shift";
        param_names__.push_back(param_name_stream__.str());
    }

}; // model

}

typedef model_crispr_model_namespace::model_crispr_model stan_model;


#endif
