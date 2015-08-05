// Copyright (c) 2011 Morten Bakkedal
// This code is published under the MIT License.

using System;
using System.Diagnostics;

using FuncLib.Optimization.Ipopt.Cureos;
using System.Collections.Generic;
using System.Globalization;

namespace FuncLib.Optimization.Ipopt
{
    public enum IpoptOptionType
    {
        String,
        Integer,
        Number
    }
    public enum IpoptOptionCategory
    {
        IpoptOutput,
        IpoptTermination,
        NLPScaling,
        NLP,
        Initialization,
        BarrierParameter,
        MultiplierUpdates,
        LineSearch,
        WarmStart,
        RestorationPhase,
        LinearSolver,
        HessianPerturbation,
        QuasiNewton,
        DerivativeTest,
        MA27LinearSolver,
        MA57LinearSolver,
        MA77LinearSolver,
        MA86LinearSolver,
        MA97LinearSolver,
        MUMPSLinearSolver,
        PardisoLinearSolver,
        WSMPLinearSolver
    }

    #region OPTION ENUMS
    public enum IpoptOutputOptions
    {
        print_level,
        print_user_options,
        print_options_documentation,
        print_frequency_iter,
        print_frequency_time,
        output_file,
        file_print_level,
        option_file_name,
        print_info_string,
        inf_pr_output,
        print_timing_statistics
    }
    public enum IpoptTerminationOptions
    {
        tol,
        max_iter,
        max_cpu_time,
        dual_inf_tol,
        constr_viol_tol,
        compl_inf_tol,
        acceptable_tol,
        acceptable_iter,
        acceptable_constr_viol_tol,
        acceptable_dual_inf_tol,
        acceptable_compl_inf_tol,
        acceptable_obj_change_tol,
        diverging_iterates_tol
    }
    public enum IpoptNLPScalingOptions
    {
        obj_scaling_factor,
        nlp_scaling_method,
        nlp_scaling_max_gradient,
        nlp_scaling_min_value
    }
    public enum IpoptNLPOptions
    {
        bound_relax_factor,
        honor_original_bounds,
        check_derivatives_for_naninf,
        nlp_lower_bound_inf,
        nlp_upper_bound_inf,
        fixed_variable_treatment,
        jac_c_constant,
        jac_d_constant,
        hessian_constant
    }
    public enum IpoptInitializationOptions
    {
        bound_frac,
        bound_push,
        slack_bound_frac,
        slack_bound_push,
        bound_mult_init_val,
        constr_mult_init_max,
        bound_mult_init_method
    }
    public enum IpoptBarrierParameterOptions
    {
        mehrotra_algorithm,
        mu_strategy,
        mu_oracle,
        quality_function_max_section_steps,
        fixed_mu_oracle,
        adaptive_mu_globalization,
        mu_init,
        mu_max_fact,
        mu_max,
        mu_min,
        mu_target,
        barrier_tol_factor,
        mu_linear_decrease_factor,
        mu_superlinear_decrease_power
    }
    public enum IpoptMultiplierUpdatesOptions
    {
        alpha_for_y,
        alpha_for_y_tol,
        recalc_y,
        recalc_y_feas_tol
    }
    public enum IpoptLineSearchOptions
    {
        max_soc,
        watchdog_shortened_iter_trigger,
        watchdog_trial_iter_max,
        accept_every_trial_step,
        corrector_type
    }
    public enum IpoptWarmStartOptions
    {
        warm_start_init_point,
        warm_start_bound_frac,
        warm_start_slack_bound_frac,
        warm_start_slack_bound_push,
        warm_start_mult_bound_push,
        warm_start_mult_init_max
    }
    public enum IpoptRestorationPhaseOptions
    {
        expect_infeasible_problem,
        expect_infeasible_problem_ctol,
        expect_infeasible_problem_ytol,
        start_with_resto,
        soft_resto_pderror_reduction_factor,
        required_infeasibility_reduction,
        bound_mult_reset_threshold,
        constr_mult_reset_threshold,
        evaluate_orig_obj_at_resto_trial
    }
    public enum IpoptLinearSolverOptions
    {
        linear_solver,
        linear_system_scaling,
        linear_scaling_on_demand,
        max_refinement_steps,
        min_refinement_steps,
    }
    public enum IpoptHessianPerturbationOptions
    {
        max_hessian_perturbation,
        min_hessian_perturbation,
        first_hessian_perturbation,
        perturb_inc_fact_first,
        perturb_inc_fact,
        perturb_dec_fact,
        jacobian_regularization_value
    }
    public enum IpoptQuasiNewtonOptions
    {
        hessian_approximation,
        limited_memory_update_type,
        limited_memory_max_history,
        limited_memory_max_skipping,
        limited_memory_initialization,
        limited_memory_init_val,
        limited_memory_init_val_max,
        limited_memory_init_val_min,
        limited_memory_special_for_resto
    }
    public enum IpoptDerivativeTestOptions
    {
        derivative_test,
        derivative_test_perturbation,
        derivative_test_tol,
        derivative_test_print_all,
        derivative_test_first_index,
        point_perturbation_radius
    }
    public enum IpoptMA27LinearSolverOptions
    {
        ma27_pivtol,
        ma27_pivtolmax,
        ma27_liw_init_factor,
        ma27_la_init_factor,
        ma27_meminc_factor
    }
    public enum IpoptMA57LinearSolverOptions
    {
        ma57_pivtol,
        ma57_pivtolmax,
        ma57_pre_alloc,
        ma57_pivot_order,
        ma57_automatic_scaling,
        ma57_block_size,
        ma57_node_amalgamation,
        ma57_small_pivot_flag
    }
    public enum IpoptMA77LinearSolverOptions
    {
        ma77_print_level,
        ma77_buffer_lpage,
        ma77_buffer_npage,
        ma77_file_size,
        ma77_maxstore,
        ma77_nemin,
        ma77_order,
        ma77_small,
        ma77_static,
        ma77_u,
        ma77_umax
    }
    public enum IpoptMA86LinearSolverOptions
    {
        ma86_print_level,
        ma86_nemin,
        ma86_order,
        ma86_scaling,
        ma86_small,
        ma86_static,
        ma86_u,
        ma86_umax
    }
    public enum IpoptMA97LinearSolverOptions
    {
        ma97_print_level,
        ma97_nemin,
        ma97_order,
        ma97_scaling,
        ma97_scaling1,
        ma97_scaling2,
        ma97_scaling3,
        ma97_small,
        ma97_solve_blas3,
        ma97_switch1,
        ma97_switch2,
        ma97_switch3,
        ma97_u,
        ma97_umax
    }
    public enum IpoptMUMPSLinearSolverOptions
    {
        mumps_pivtol,
        mumps_pivtolmax,
        mumps_mem_percent,
        mumps_permuting_scaling,
        mumps_pivot_order,
        mumps_scaling
    }
    public enum IpoptPardisoLinearSolverOptions
    {
        pardiso_msglvl,
        pardiso_matching_strategy,
        pardiso_out_of_core_power
    }
    public enum IpoptWSMPLinearSolverOptions
    {
        wsmp_ordering_option,
        wsmp_pivtol,
        wsmp_pivtolmax,
        wsmp_scaling
    }
    #endregion

    public static class IpoptOptionManager
    {
        #region GENERIC OPTIONS
        public static IpoptOption GenericUnsafeOption(string name, string value)
        {
            return new IpoptStringOption(name, "generic unsafe string option", "", value, value);
        }
        public static IpoptOption GenericUnsafeOption(string name, int value)
        {
            return new IpoptIntegerOption(name, "generic unsafe integer option", "", value, value);
        }
        public static IpoptOption GenericUnsafeOption(string name, double value)
        {
            return new IpoptNumberOption(name, "generic unsafe number option", "", value, value);
        }
        #endregion

        #region OUTPUT OPTIONS
        /// <summary>
        /// print_level : output verbosity level.
        /// Sets the default verbosity level for console output.
        /// The larger this value the more detailed is the output.
        /// The valid range for this integer option is 0 ≤ print_level ≤ 12 and its default value is 5.
        /// </summary>
        public static IpoptOption print_level { get { return _print_level.Clone(); } }
        private static IpoptOption _print_level = new IpoptIntegerOption(
            "print_level",
            "Output verbosity level. Sets the default verbosity level for console output. " +
            "The larger this value the more detailed is the output.",
            "[0,12]", 5);

        /// <summary>
        /// print_user_options : print all options set by the user.
        /// If selected, the algorithm will print the list of all options set by the user including their values
        /// and whether they have been used. The default value for this string option is ”no”.
        /// </summary>
        /// <list type="bullet">
        /// Possible values: 
        /// <item><term>no</term><description>don’t print options</description></item>
        /// <item><term>yes</term><description>print options</description></item>
        /// </list>
        public static IpoptOption print_user_options { get { return _print_user_options.Clone(); } }
        private static IpoptOption _print_user_options = new IpoptStringOption(
            "print_user_options",
            "Print all options set by the user. " +
            "If selected, the algorithm will print the list of all options set by the user including their values and whether they have been used." +
            "In some cases this information might be incorrect, due to the internal program flow.",
            "[no,yes]", "no");

        public static IpoptOption print_options_documentation { get { return _print_options_documentation.Clone(); } }
        private static IpoptOption _print_options_documentation = new IpoptStringOption(
            "print_options_documentation",
            "Switch to print all algorithmic options. " +
            "If selected, the algorithm will print the list of all available algorithmic options with some documentation before solving the optimization problem.",
            "[no,yes]", "no");
        //print_options_documentation,
        //print_frequency_iter,
        //print_frequency_time,
        //output_file,
        //file_print_level,
        //option_file_name,
        //print_info_string,
        //inf_pr_output,
        //print_timing_statistics
        #endregion

        #region TERMINATION OPTIONS
        public static IpoptOption tol { get { return _tol.Clone(); } }
        private static IpoptOption _tol = new IpoptNumberOption(
           "tol",
           "Desired convergence tolerance (relative). " +
           "Determines the convergence tolerance for the algorithm. " +
           "The algorithm terminates successfully, if the (scaled) NLP error becomes smaller than this value, " +
           "and if the (absolute) criteria according to \"dual inf tol\", \"primal inf tol\", and \"compl inf tol\" are met",
           "]0,+inf[", 1.0e-8);

        public static IpoptOption max_iter { get { return _max_iter.Clone(); } }
        private static IpoptOption _max_iter = new IpoptIntegerOption(
           "max_iter",
           "Maximum number of iterations. " +
           "The algorithm terminates with an error message if the number of iterations exceeded this number",
           "[0,+inf[", 3000);

        public static IpoptOption max_cpu_time { get { return _max_cpu_time.Clone(); } }
        private static IpoptOption _max_cpu_time = new IpoptNumberOption(
           "max_cpu_time",
           "Maximum number of CPU seconds. " +
           "A limit on CPU seconds that Ipopt can use to solve one problem. " +
           "If during the convergence check this limit is exceeded, Ipopt will terminate with a corresponding error message.",
           "]0,+inf[", 1.0e6);
        //,
        //dual_inf_tol,
        //constr_viol_tol,
        //compl_inf_tol,
        //acceptable_tol,
        //acceptable_iter,
        //acceptable_constr_viol_tol,
        //acceptable_dual_inf_tol,
        //acceptable_compl_inf_tol,
        //acceptable_obj_change_tol,
        //diverging_iterates_tol
        #endregion

        #region NLP SCALING OPTIONS
        public static IpoptOption obj_scaling_factor { get { return _obj_scaling_factor.Clone(); } }
        private static IpoptOption _obj_scaling_factor = new IpoptNumberOption(
            "obj_scaling_factor",
            "Scaling factor for the objective function. " +
            "This option sets a scaling factor for the objective function. " +
            "The scaling is seen internally by Ipopt but the unscaled objective is reported in the console output. " +
            "If additional scaling parameters are computed (e.g. user-scaling or gradient-based), both factors are multiplied. " +
            "If this value is chosen to be negative, Ipopt will maximize the objective function instead of minimizing it.",
            "]-inf,+inf[", 1);
        //nlp_scaling_method,
        //nlp_scaling_max_gradient,
        //nlp_scaling_min_value

        #endregion

        #region QUASI-NEWTON OPTIONS
        /// <summary>
        /// hessian_approximation : indicates what Hessian information is to be used.
        /// This determines which kind of information for the Hessian of the Lagrangian function is used by the algorithm. 
        /// The default value for this string option is ”exact”.
        /// </summary>
        /// <list type="bullet">
        /// Possible values: 
        /// <item><term>exact</term><description>Use second derivatives provided by the NLP.</description></item>
        /// <item><term>limited-memory:</term><description>Perform a limited-memory quasi-Newton approximation</description></item>
        /// </list>
        public static IpoptOption hessian_approximation { get { return _hessian_approximation.Clone(); } }
        private static IpoptOption _hessian_approximation = new IpoptStringOption(
            "hessian_approximation", "",
            "[exact,limited-memory]", "exact");
        //limited_memory_update_type,
        //limited_memory_max_history,
        //limited_memory_max_skipping,
        //limited_memory_initialization,
        //limited_memory_init_val,
        //limited_memory_init_val_max,
        //limited_memory_init_val_min,
        //limited_memory_special_for_resto
        #endregion
    }

    #region OPTION CLASSES

    [Serializable]
    public class IpoptOptionCollection : List<IpoptOption>
    {
        internal IpoptOptionCollection()
        {
        }

        public bool Add(IpoptOption option)
        {
            if (!Contains(option.Name))
            {
                base.Add(option);
                return true;
            }
            else { return false; }
        }
        public bool Add(IpoptOption option, object value)
        {
            if (Add(option))
            {
                option.Value = value;
                return true;
            }
            else { return false; }
        }
        public void Remove(string name)
        {
            RemoveAll(delegate(IpoptOption option) { return option.Name == name; });
        }
        public bool Contains(string name)
        {
            return Find(delegate(IpoptOption option) { return option.Name == name; }) != null;
        }
        public IpoptOption this[string name]
        {
            get
            {
                return Find(delegate(IpoptOption option) { return option.Name == name; });
            }
        }
    }

    [Serializable]
    [DebuggerDisplay("{ToString(),nq}")]
    public abstract class IpoptOption
    {
        protected IpoptOptionType _type;
        private string _name;
        private string _description;
        private string _domain;

        public IpoptOptionType Type
        {
            get { return _type; }
        }
        public string Name
        {
            get { return _name; }
        }
        public string Description
        {
            get { return _description; }
        }
        public string Domain
        {
            get { return _domain; }
        }
        public abstract object DefaultValue { get; }
        public abstract object Value { get; set; }

        public IpoptOption(string name, string description,  string domain)
        {
            this._name = name;
            this._description = description;
            this._domain = domain;
        }

        internal abstract void Prepare(IpoptProblem ipopt);
        internal abstract IpoptOption Clone();
        public abstract bool IsValueValid();
        public override string ToString()
        {
            return Name + " = " + Value;
        }
    }

    [Serializable]
    public class IpoptStringOption : IpoptOption
    {
        // "[yes,no]" with no spaces
        private List<string> _valid_set;
        private string _default_value;
        private string _value;

        public override object Value
        {
            get { return _value; }
            set { this._value = (string)value; }
        }
        public override object DefaultValue
        {
            get { return _default_value; }
        }

        public IpoptStringOption(string name, string description, string domain, string default_value, string value = "")
            : base(name, description, domain)
        {
            this._default_value = default_value;
            if (value != "") { this.Value = value; }
            _valid_set = new List<string>(domain.Substring(1, domain.Length - 2).Split(','));
        }

        internal override void Prepare(IpoptProblem ipopt)
        {
            ipopt.AddOption(Name, _value);
        }
        internal override IpoptOption Clone()
        {
            return new IpoptStringOption(Name, Description, Domain, _default_value, _value);
        }
        public override bool IsValueValid()
        {
            if (Domain == "") return true;
            return _valid_set.Exists(str => str == Value);
        }

    }

    [Serializable]
    public class IpoptIntegerOption : IpoptOption
    {
        // [0,150] or [-inf,10[ or ]-inf,10]
        private char _lower_bracket;
        private char _upper_bracket;
        private string _valid_lower_bound;
        private string _valid_upper_bound;
        private int _default_value;
        private int _value;

        public override object Value
        {
            get { return _value; }
            set { this._value = (int)value; }
        }
        public override object DefaultValue
        {
            get { return _default_value; }
        }

        public IpoptIntegerOption(string name, string description, string domain, int default_value)
            : base(name, description, domain)
        {
            this._default_value = default_value;

            _lower_bracket = domain[0];
            _upper_bracket = domain[domain.Length - 1];

            string[] bounds = domain.Substring(1, domain.Length - 2).Split(',');
            _valid_lower_bound = bounds[0];
            _valid_upper_bound = bounds[1];
        }
        public IpoptIntegerOption(string name, string description, string domain, int default_value, int value)
            : this(name, description, domain, default_value)
        {
            this.Value = value;
        }

        internal override void Prepare(IpoptProblem ipopt)
        {
            ipopt.AddOption(Name, _value);
        }
        internal override IpoptOption Clone()
        {
            return new IpoptIntegerOption(Name, Description, Domain, _default_value, _value);
        }
        public override bool IsValueValid()
        {
            if (Domain == "") return true;

            bool lower = false;
            bool upper = false;
            int lower_bound = Int32.Parse(_valid_lower_bound, NumberStyles.Float, CultureInfo.InvariantCulture);
            int upper_bound = Int32.Parse(_valid_upper_bound, NumberStyles.Float, CultureInfo.InvariantCulture);

            if (_valid_lower_bound != "-inf")
            {
                if (_lower_bracket == '[')
                {
                    lower = (lower_bound <= _value);
                }
                if (_lower_bracket == ']')
                {
                    lower = (lower_bound < _value);
                }
            }
            else { lower = true; }

            if (_valid_lower_bound != "+inf")
            {
                if (_upper_bracket == ']')
                {
                    upper = (_value <= upper_bound);
                }
                if (_lower_bracket == '[')
                {
                    upper = (_value < upper_bound);
                }
            }
            else { upper = true; }

            return (lower && upper);
        }
    }

    [Serializable]
    public class IpoptNumberOption : IpoptOption
    {
        // [0.27923,150.09] or [-inf,10.98[ or ]-inf,10.210912]
        private char _lower_bracket;
        private char _upper_bracket;
        private string _valid_lower_bound;
        private string _valid_upper_bound;
        private double _default_value;
        private double _value;

        public override object Value
        {
            get { return _value; }
            set { this._value = (double)value; }
        }
        public override object DefaultValue
        {
            get { return _default_value; }
        }

        public IpoptNumberOption(string name, string description, string domain, double default_value)
            : base(name, description, domain)
        {
            this._default_value = default_value;

            _lower_bracket = domain[0];
            _upper_bracket = domain[domain.Length - 1];

            string[] bounds = domain.Substring(1, domain.Length - 2).Split(',');
            _valid_lower_bound = bounds[0];
            _valid_upper_bound = bounds[1];
        }
        public IpoptNumberOption(string name, string description, string domain, double default_value, double value)
            : this(name, description, domain, default_value)
        {
            this._value = value;
        }

        internal override void Prepare(IpoptProblem ipopt)
        {
            ipopt.AddOption(Name, _value);
        }
        internal override IpoptOption Clone()
        {
            return new IpoptNumberOption(Name, Description, Domain, _default_value, _value);
        }
        public override bool IsValueValid()
        {
            if (Domain == "") return true;

            bool lower = false;
            bool upper = false;
            double lower_bound = Double.Parse(_valid_lower_bound, NumberStyles.Float, CultureInfo.InvariantCulture);
            double upper_bound = Double.Parse(_valid_upper_bound, NumberStyles.Float, CultureInfo.InvariantCulture);

            if (_valid_lower_bound != "-inf")
            {
                if (_lower_bracket == '[')
                {
                    lower = (lower_bound <= _value);
                }
                if (_lower_bracket == ']')
                {
                    lower = (lower_bound < _value);
                }
            }
            else { lower = true; }

            if (_valid_lower_bound != "+inf")
            {
                if (_upper_bracket == ']')
                {
                    upper = (_value <= upper_bound);
                }
                if (_lower_bracket == '[')
                {
                    upper = (_value < upper_bound);
                }
            }
            else { upper = true; }

            return (lower && upper);
        }
    }

    #endregion

    
}
