import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.integrate import quad
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter, MultipleLocator

# Read the data
## US and states
df = pd.read_csv('data/mobility_t_values_states.csv')
## county
# df = pd.read_csv('data/mobility_t_values_key_county.csv')
## routes
# df = pd.read_csv('data/average_weekday_boardings_rapidride.csv')


# Define V-Shape
def Q_V_Shape(t, a_V, b_V, t0, t2, t4, Q_t4):
    # Apply logic element-wise
    b_V = (Q_t4-a_V*(t2-t0))/(t4-t2)
    Q_D = a_V * (t-t0)  # Disruption phase equation
    Q_R = a_V * (t2-t0) - b_V * (t2-t0) + b_V * (t-t0)   # Recovery phase equation
    Q_normal = a_V * t2 - b_V * t2 + b_V * t4
    return np.where((t0 <= t) & (t <= t2), Q_D, np.where((t2 <= t) & (t <= t4), Q_R, Q_R))

def pi_V_Shape(t, a_V, b_V, t0, t2, t4):
    # Apply logic element-wise
    return np.where((t0 <= t) & (t <= t2), a_V, np.where((t2 <= t) & (t <= t4), b_V, b_V))


# Define the objective function for V-shape fitting
def objective_V_Shape(params, t, Q_observed):
    a_V, b_V = params
    Q_V_estimated = Q_V_Shape(t, a_V, b_V, t0, t2, t4, Q_t4)
    return np.sum((Q_observed - Q_V_estimated) ** 2)

def Q_piecewise_exp(t, A, k, B, A_prime, k_prime, t0, t2, t4, Q_t4):
    # # Calculate A' based on the boundary condition at t4
    # denominator = np.exp(-k_prime * (t4 - t0)) - np.exp(-k_prime * (t2 - t0))
    # if denominator == 0:
    #     raise ValueError("Denominator in expression for A' results in division by zero.")
    #
    # A_prime = (Q_t4 - A * np.exp(-k * (t2 - t0)) - B) / denominator

    Q_D = A * np.exp(-k * (t-t0)) + B
    Q_R = A_prime * np.exp(-k_prime * (t-t0)) + A * np.exp(-k * (t2-t0)) + B - A_prime * np.exp(-k_prime * (t2-t0))
    # Apply the correct piece based on the interval
    return np.where((t0 <= t) & (t <= t2), Q_D, np.where((t2 <= t) & (t <= t4), Q_R, Q_R))

def pi_piecewise_exp(t, A, k, B, A_prime, k_prime, t0, t2, t4):
    pi_D = -A * k * np.exp(-k * (t-t0))
    pi_R = -A_prime * k_prime * np.exp(-k_prime * (t-t0))
    # Apply the correct piece based on the interval
    return np.where((t0 <= t) & (t <= t2), pi_D, np.where((t2 <= t) & (t <= t4), pi_R, pi_R))

# Objective function for the piecewise exponential fitting
def objective_piecewise_exp(params, t, Q_observed):
    A, k, B, A_prime, k_prime = params
    Q_piecewise_estimated = Q_piecewise_exp(t, A, k, B, A_prime, k_prime, t0, t2, t4, Q_t4)
    return np.sum((Q_observed - Q_piecewise_estimated) ** 2)


# Define the DQQ
def Q_DQQ(t, alpha, beta, t0, t2, t4, Q_t4):
    beta = -6*(Q_t4 + (alpha / 6) * (t2-t0)**3)/(t4-t2)**3
    Q_D = (alpha / 6) * (t-t0)**2 * (2*(t-t0) - 3*(t2-t0))
    Q_R = -(alpha / 6) * (t2-t0)**3 - beta * (t - t2)**2 * ((t4 - t2)/2 - (t - t2)/3)
    # Q_normal = -(alpha / 6) * (t2-t0)**3 - beta * (t4 - t2)**2 * ((t4 - t2)/2 - (t4 - t2)/3)
    Q_normal = Q_t4
    return np.where((t0 <= t) & (t <= t2), Q_D, np.where((t2 <= t) & (t <= t4), Q_R, Q_normal))

def pi_DQQ(t, alpha, beta, t0, t2, t4):
    pi_D = alpha * (t - t0) * (t - t2)
    pi_R = beta * (t - t2) * (t - t4)
    return np.where((t0 <= t) & (t <= t2), pi_D, np.where((t2 <= t) & (t <= t4), pi_R, 0))


# Objective function for fitting
def objective_DQQ(params, t, Q_observed, pi_observed):
    theta =0.5
    alpha, beta = params
    Q_DQQ_estimated = Q_DQQ(t, alpha, beta, t0, t2, t4, Q_t4)
    pi_DQQ_estimated = pi_DQQ(t, alpha, beta, t0, t2, t4)
    return theta*np.sum((Q_observed - Q_DQQ_estimated) ** 2) + (1-theta)*np.sum((pi_observed - pi_DQQ_estimated) ** 2)

# Function to calculate R²
def calculate_r2(observed, predicted):
    SSR = np.sum((observed - predicted) ** 2)
    SST = np.sum((observed - np.mean(observed)) ** 2)
    return 1 - (SSR / SST)

# Function to calculate Root Mean Square Error (RMSE)
def calculate_rmse(observed, predicted):
    return np.sqrt(np.mean((observed - predicted) ** 2))

# Function to calculate Akaike Information Criterion (AIC)
def calculate_aic(observed, predicted, num_params):
    rss = np.sum((observed - predicted) ** 2)
    n = len(observed)
    return n * np.log(rss / n) + 2 * num_params

# Define a formatter function
def two_decimal_formatter(x, pos):
    return f'{x:.2f}'

def three_decimal_formatter(x, pos):
    return f'{x:.3f}'

# Initialize empty lists to store results
states = []
mobility_types = []
# -------------------
a_V_Shape = []
b_V_Shape = []
r2_pi_V_Shape = []
r2_pi_V_Shape_disruption = []
r2_pi_V_Shape_recovery = []
r2_Q_V_Shape = []
r2_Q_V_Shape_disruption = []
r2_Q_V_Shape_recovery = []
# -------------------
A_list = []
k_list = []
B_list = []
A_prime_list = []
k_prime_list = []
r2_pi_exp = []
r2_pi_exp_disruption = []
r2_pi_exp_recovery = []
r2_Q_exp = []
r2_Q_exp_disruption = []
r2_Q_exp_recovery = []
# -------------------
alpha_DQQ = []
beta_DQQ = []
r2_pi_DQQ = []
r2_pi_DQQ_disruption = []
r2_pi_DQQ_recovery = []
r2_Q_DQQ = []
r2_Q_DQQ_disruption = []
r2_Q_DQQ_recovery = []
Q_observed_t2_list = []
Q_observed_t4_list = []
Q_predicted_DQQ_t2_list = []
Q_predicted_DQQ_t4_list = []
# -------------------
rmse_Q_V_Shape = []
rmse_pi_V_Shape = []
rmse_Q_piecewise_exp = []
rmse_pi_piecewise_exp = []
rmse_Q_DQQ = []
rmse_pi_DQQ = []
# -------------------
aic_Q_V_Shape = []
aic_pi_V_Shape = []
aic_Q_piecewise_exp = []
aic_pi_piecewise_exp = []
aic_Q_DQQ = []
aic_pi_DQQ = []



# Filter for Washington state and residential, workplace mobility types
# df = df[(df['state'] == 'Washington') & (df['mobility_type'].isin(['residential_percent_change_from_baseline', 'workplaces_percent_change_from_baseline']))]
df = df[(df['state'] == 'US') & (df['mobility_type'].isin(['transit_stations_percent_change_from_baseline']))]
# df = df[(df['state'] == 'Washington') & (df['mobility_type'].isin(['residential_percent_change_from_baseline', 'workplaces_percent_change_from_baseline']))]
# US
# California, Arizona, and Washington
# 'San Francisco County', 'Maricopa County', 'King County'
# A Line Northbound
# A Line Southbound
# B Line Eastbound
# B Line Westbound
# C Line Northbound
# C Line Southbound
# D Line Northbound
# D Line Southbound
# E Line Northbound
# E Line Southbound
# F Line Eastbound
# F Line Westbound

# Loop over unique combinations of state and mobility_type
for (state, mobility_type), group in df.groupby(['state', 'mobility_type']):
    t = np.array(np.arange(len(group['t'])))
    t_lst = list(np.array(np.arange(len(group['t']))))

    t0 = group['t0'].values[0]
    t2 = group['t2'].values[0]
    t4 = group['t4'].values[0]

    Q_observed = group['ridership'].values
    Q_t4 = Q_observed[int(t4)]
    # Calculate the gradient of the ridership values
    pi_observed = np.gradient(Q_observed)
    pi_t4 = pi_observed[int(t4)]

    # Filter time and observations between t0 and t4
    mask = (t >= t0) & (t <= t4+40)
    t_fit = t[mask]
    t_fit_lst = list(t_fit)
    model_Q_observed = group['ridership'].values[mask]
    model_pi_observed = np.gradient(model_Q_observed)

    # Fitting V-Shape model
    init_params_v_shape = [1, 1]  # Initial guess for A (adjust as needed)
    res_v_shape = minimize(objective_V_Shape, init_params_v_shape, args=(t_fit, model_Q_observed))

    # Initial guesses for parameters
    init_params_piecewise_exp = [1, 0.1, 0, 1, 0.1]  # Adjust as necessary
    # Perform optimization for the piecewise exponential model
    res_piecewise_exp = minimize(objective_piecewise_exp, init_params_piecewise_exp, args=(t_fit, model_Q_observed))


    # Fitting DQQ model
    # Initial guesses for alpha and beta
    init_params_DQQ = [1, 1]
    # Perform optimization to fit the custom hazard model
    res_DQQ = minimize(objective_DQQ, init_params_DQQ, args=(t_fit, model_Q_observed, model_pi_observed))  # Note: upper bound for s is determined by top 5% of s

    # Extract the fitted parameters
    a_V, b_V = res_v_shape.x
    # Extract the fitted parameters
    A, k, B, A_prime, k_prime = res_piecewise_exp.x
    # A_log_logistic, B_log_logistic = res_log_logistic.x
    # A_weibull, B_weibull = res_weibull.x
    alpha, beta = res_DQQ.x

    # Calibrated model
    Q_predicted_V_Shape = Q_V_Shape(t_fit, a_V, b_V, t0, t2, t4, Q_t4)
    # r_predicted_log_logistic = Log_Logistic(t, A_log_logistic, B_log_logistic)
    # r_predicted_Weibull = weibull(t, A_weibull, B_weibull)
    Q_predicted_DQQ = Q_DQQ(t_fit, alpha, beta, t0, t2, t4, Q_t4)
    pi_predicted_V_Shape = pi_V_Shape(t_fit, a_V, b_V, t0, t2, t4)
    pi_predicted_DQQ = pi_DQQ(t_fit, alpha, beta, t0, t2, t4)
    # Predicted models
    Q_predicted_piecewise_exp = Q_piecewise_exp(t_fit, A, k, B, A_prime, k_prime, t0, t2, t4, Q_t4)
    pi_predicted_piecewise_exp = pi_piecewise_exp(t_fit, A, k, B, A_prime, k_prime, t0, t2, t4)

    # Assuming you know the indices for t0, t2 within t_fit or recalculate them if t0, t2 are actual time points
    index_t2 = np.searchsorted(t_fit, t2)  # Find index of t2 within t_fit

    # Calculate R² for V-Shape
    R2_Q_v_shape = calculate_r2(model_Q_observed, Q_predicted_V_Shape)
    # Calculate R² for intervals within the filtered range
    R2_Q_v_shape_disruption = calculate_r2(model_Q_observed[:index_t2], Q_predicted_V_Shape[:index_t2])
    R2_Q_v_shape_recovery = calculate_r2(model_Q_observed[index_t2:], Q_predicted_V_Shape[index_t2:])
    # Calculate R² for V-Shape
    # R2_pi_v_shape = calculate_r2(pi_observed[5:], pi_predicted_V_Shape[5:])
    R2_pi_v_shape = calculate_r2(model_pi_observed, pi_predicted_V_Shape)
    # Calculate R² for V-Shape in intervals
    # R2_pi_v_shape_disruption = calculate_r2(pi_observed[5:int(t2)], pi_predicted_V_Shape[5:int(t2)])
    R2_pi_v_shape_disruption = calculate_r2(model_pi_observed[:index_t2], pi_predicted_V_Shape[:index_t2])
    R2_pi_v_shape_recovery = calculate_r2(model_pi_observed[index_t2:], pi_predicted_V_Shape[index_t2:])


    # Calculate R² for the piecewise exponential model
    R2_Q_piecewise_exp = calculate_r2(model_Q_observed, Q_predicted_piecewise_exp)
    # Calculate R² for exp in intervals
    R2_Q_piecewise_exp_disruption = calculate_r2(model_Q_observed[:index_t2], Q_predicted_piecewise_exp[:index_t2])
    R2_Q_piecewise_exp_recovery = calculate_r2(model_Q_observed[index_t2:], Q_predicted_piecewise_exp[index_t2:])
    # Calculate R² for exp
    R2_pi_piecewise_exp = calculate_r2(model_pi_observed, pi_predicted_piecewise_exp)
    # Calculate R² for exp in intervals
    R2_pi_piecewise_exp_disruption = calculate_r2(model_pi_observed[:index_t2], pi_predicted_piecewise_exp[:index_t2])
    R2_pi_piecewise_exp_recovery = calculate_r2(model_pi_observed[index_t2:], pi_predicted_piecewise_exp[index_t2:])

    # Calculate R² for DQQ
    R2_Q_DQQ = calculate_r2(model_Q_observed, Q_predicted_DQQ)
    # Calculate R² for DQQ in intervals
    R2_Q_DQQ_disruption = calculate_r2(model_Q_observed[:index_t2], Q_predicted_DQQ[:index_t2])
    R2_Q_DQQ_recovery = calculate_r2(model_Q_observed[index_t2:], Q_predicted_DQQ[index_t2:])
    # Calculate R² for DQQ
    # R2_pi_DQQ = calculate_r2(pi_observed[5:], pi_predicted_DQQ[5:])
    R2_pi_DQQ = calculate_r2(model_pi_observed, pi_predicted_DQQ)
    # Calculate R² for DQQ in intervals
    # R2_pi_DQQ_disruption = calculate_r2(pi_observed[5:int(t2)], pi_predicted_DQQ[5:int(t2)])
    R2_pi_DQQ_disruption = calculate_r2(model_pi_observed[:index_t2], pi_predicted_DQQ[:index_t2])
    R2_pi_DQQ_recovery = calculate_r2(model_pi_observed[index_t2:], pi_predicted_DQQ[index_t2:])


    # Extract values at t2 and t4
    Q_observed_t2 = Q_observed[int(t2)]
    Q_observed_t4 = Q_observed[int(t4)]
    Q_predicted_DQQ_t2 = Q_DQQ(t2, alpha, beta, t0, t2, t4, Q_t4)
    Q_predicted_DQQ_t4 = Q_DQQ(t4, alpha, beta, t0, t2, t4, Q_t4)

    # RMSE calculations
    RMSE_Q_V_Shape = calculate_rmse(model_Q_observed, Q_predicted_V_Shape)
    RMSE_pi_V_Shape = calculate_rmse(model_pi_observed, pi_predicted_V_Shape)
    RMSE_Q_piecewise_exp = calculate_rmse(model_Q_observed, Q_predicted_piecewise_exp)
    RMSE_pi_piecewise_exp = calculate_rmse(model_pi_observed, pi_predicted_piecewise_exp)
    RMSE_Q_DQQ = calculate_rmse(model_Q_observed, Q_predicted_DQQ)
    RMSE_pi_DQQ = calculate_rmse(model_pi_observed, pi_predicted_DQQ)

    # AIC calculations
    AIC_Q_V_Shape = calculate_aic(model_Q_observed, Q_predicted_V_Shape, 2)  # Assuming 2 parameters for V-Shape
    AIC_pi_V_Shape = calculate_aic(model_pi_observed, pi_predicted_V_Shape, 2)  # Assuming 2 parameters for V-Shape
    AIC_Q_piecewise_exp = calculate_aic(model_Q_observed, Q_predicted_piecewise_exp, 5)  # Assuming 5 parameters for Piecewise Exponential
    AIC_pi_piecewise_exp = calculate_aic(model_pi_observed, pi_predicted_piecewise_exp,5)  # Assuming 5 parameters for Piecewise Exponential
    AIC_Q_DQQ = calculate_aic(model_Q_observed, Q_predicted_DQQ, 2)  # Assuming 2 parameters for DQQ
    AIC_pi_DQQ = calculate_aic(model_pi_observed, pi_predicted_DQQ, 2)  # Assuming 2 parameters for DQQ

    # Create a single figure with 2 subplots
    # fig_1, axs = plt.subplots(2, 1, figsize=(16, 24), sharex=True)
    fig, axs = plt.subplots(2, 1, figsize=(32, 24))
    plt.subplots_adjust(hspace=0.3)  # Adjust vertical spacing between subplots
    plt.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.2)  # Adjust the margins

    key_times = [t0, t2, t4]  # example key times
    # for ktime in key_times:
    #     ymax_value_Q = Q_DQQ(ktime, alpha, beta, t0, t2, t4, Q_t4)  # Get the ymax value from your function
    #     ymax_value_pi = pi_DQQ(ktime, alpha, beta, t0, t2, t4)
    #     # Use vlines instead of axvline to specify actual data coordinates for ymin and ymax
    #     axs[0].vlines(x=ktime, ymin=-60, ymax=ymax_value_Q, colors='r', linestyles='--', linewidth=4, alpha=0.5)
    #     axs[1].vlines(x=ktime, ymin=-60, ymax=ymax_value_pi, colors='r', linestyles='--', linewidth=4, alpha=0.5)

    for ax in axs:
        for ktime in key_times:
            ax.axvline(ktime, color='r', linestyle='--',linewidth=4, alpha=0.5)
        # # Set x-ticks manually if needed
        # ax.set_xticks([t0, t2, t4])
        # # Custom labels
        # labels = ['$t_0$', '$t_2$', '$t_4$']
        # ax.set_xticklabels(labels, fontsize=24)  # Apply custom labels

        ax.legend(loc='lower right', fontsize=20)
        ax.grid(True)

        # Setting x-axis tick intervals to 25
        ax.xaxis.set_major_locator(MultipleLocator(25))

    # Create predictions for plotting
    axs[0].set_xlim(0, 250)
    axs[0].set_ylim(-60, 15)
    axs[1].set_xlim(0, 250)
    axs[1].set_ylim(-6, 5)

    # Subplot 1: S_curve and Observed Data
    axs[0].plot(t, Q_observed, 'o-', label='Observed Data', markersize=15, linewidth=15, alpha=0.3)
    axs[0].plot(t_fit, Q_predicted_V_Shape, 'x-', label=f'V-Shape, $R^2$: {R2_Q_v_shape:.3f}')
    axs[0].plot(t_fit, Q_predicted_piecewise_exp, 'x-', label=f'Exp, $R^2$: {R2_Q_piecewise_exp:.3f}')
    # axs[0].plot(t, F_predicted_log_logistic, 'x-', label='Log Logistic')
    # axs[0].plot(t, F_predicted_Weibull, 'x-', label='Weibull')
    axs[0].plot(t_fit, Q_predicted_DQQ, 'x-', label=f'DQQ, $R^2$: {R2_Q_DQQ:.3f}', linewidth=6)
    # axs[0].set_xlabel('Time (Days)\n\n(a)', fontsize=30, labelpad=20)
    axs[0].set_xlabel('Time (Days)', fontsize=30, labelpad=20)
    # Apply the formatter to the y-axis
    axs[0].yaxis.set_major_formatter(FuncFormatter(two_decimal_formatter))
    axs[0].set_ylabel('Q(t)', fontsize=30, labelpad=20)
    axs[0].tick_params(axis='x', labelsize=27.5)  # Set font size for x-axis labels
    axs[0].tick_params(axis='y', labelsize=27.5)  # Set font size for y-axis labels
    # axs[0].set_title(f'Observed vs. Modeled for {state}, {mobility_type}', fontsize=27.5)
    axs[0].legend(loc='lower right', fontsize=27.5)
    axs[0].grid(True)

    # Subplot 2: Hazard Model and Observed Data
    axs[1].plot(t, pi_observed, 'o-', label='Observed Data', markersize=15, linewidth=15, alpha=0.3)
    axs[1].plot(t_fit, pi_predicted_V_Shape, 'x-', label=f'V-Shape, $R^2$: {R2_pi_v_shape:.3f}')
    axs[1].plot(t_fit, pi_predicted_piecewise_exp, 'x-', label=f'Exp, $R^2$: {R2_pi_piecewise_exp:.3f}')
    # axs[1].plot(t, r_predicted_log_logistic, 'x-', label='Log Logistic')
    # axs[1].plot(t, r_predicted_Weibull, 'x-', label='Weibull')
    axs[1].plot(t_fit, pi_predicted_DQQ, 'x-', label=f'DQQ, $R^2$: {R2_pi_DQQ:.3f}', linewidth=6)
    # axs[1].set_xlabel('Time (Days)\n\n(b)\n\n(2) Residential in Washington', fontsize=30, labelpad=20)
    axs[1].set_xlabel('Time (Days)', fontsize=30, labelpad=20)
    axs[1].set_ylabel(f'$\pi(t)$', fontsize=30, labelpad=20)
    axs[1].tick_params(axis='x', labelsize=27.5)  # Set font size for x-axis labels
    # Apply the formatter to the y-axis
    axs[1].yaxis.set_major_formatter(FuncFormatter(three_decimal_formatter))
    axs[1].tick_params(axis='y', labelsize=27.5)  # Set font size for y-axis labels
    # axs[1].set_title(f'Observed vs. Modeled for {state}, {mobility_type}')
    axs[1].legend(loc='lower right', fontsize=27.5)
    axs[1].grid(True)


    # Save the plot
    filename = f"fig/DQQ_{state}_{mobility_type}.png"
    plt.savefig(filename, format='png')

    # plt.tight_layout()
    plt.show()

    # Create a subplot to zoom in extended days
    fig, axs = plt.subplots(2, 1, figsize=(16, 24))
    plt.subplots_adjust(hspace=0.3)  # Adjust vertical spacing between subplots
    plt.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.2)  # Adjust the margins
    # Customize subplots
    for ax in axs:
        key_times = [t0, t2, t4]
        for ktime in key_times:
            ax.axvline(ktime, color='r', linestyle='--', linewidth=4, alpha=0.5)

        # ax.legend(loc='lower right', fontsize=20)
        ax.grid(True)

        # Setting x-axis tick intervals to 25
        ax.xaxis.set_major_locator(MultipleLocator(5))

    # Create predictions for plotting
    axs[0].set_xlim(t4-5, t_fit[-1]+5)
    axs[0].set_ylim(Q_t4-10, Q_t4+10)
    axs[1].set_xlim(t4-5, t_fit[-1]+5)
    axs[1].set_ylim(-1.5, 1.5)


    # Subplot 1: S_curve and Observed Data
    axs[0].plot(t, Q_observed, 'o-', label='Observed Data', markersize=15, linewidth=15, alpha=0.3)
    axs[0].plot(t_fit, Q_predicted_V_Shape, 'x-', label=f'V-Shape, $R^2$: {R2_Q_v_shape:.3f}', linewidth=4)
    axs[0].plot(t_fit, Q_predicted_piecewise_exp, 'x-', label=f'Exp, $R^2$: {R2_Q_piecewise_exp:.3f}', linewidth=4)
    # axs[0].plot(t, F_predicted_log_logistic, 'x-', label='Log Logistic')
    # axs[0].plot(t, F_predicted_Weibull, 'x-', label='Weibull')
    axs[0].plot(t_fit, Q_predicted_DQQ, 'x-', label=f'DQQ, $R^2$: {R2_Q_DQQ:.3f}', linewidth=8)
    # axs[0].set_xlabel('Time (Days)\n\n(a)', fontsize=30, labelpad=20)
    axs[0].set_xlabel('Time (Days)', fontsize=30, labelpad=20)
    # Apply the formatter to the y-axis
    axs[0].yaxis.set_major_formatter(FuncFormatter(two_decimal_formatter))
    axs[0].set_ylabel('Q(t)', fontsize=30, labelpad=20)
    axs[0].tick_params(axis='x', labelsize=27.5)  # Set font size for x-axis labels
    axs[0].tick_params(axis='y', labelsize=27.5)  # Set font size for y-axis labels
    # axs[0].set_title(f'Observed vs. Modeled for {state}, {mobility_type}', fontsize=27.5)
    # axs[0].legend(loc='upper right', fontsize=27.5)
    axs[0].grid(True)

    # Subplot 2: Hazard Model and Observed Data
    axs[1].plot(t, pi_observed, 'o-', label='Observed Data', markersize=15, linewidth=15, alpha=0.3)
    axs[1].plot(t_fit, pi_predicted_V_Shape, 'x-', label=f'V-Shape, $R^2$: {R2_pi_v_shape:.3f}', linewidth=4)
    axs[1].plot(t_fit, pi_predicted_piecewise_exp, 'x-', label=f'Exp, $R^2$: {R2_pi_piecewise_exp:.3f}', linewidth=4)
    # axs[1].plot(t, r_predicted_log_logistic, 'x-', label='Log Logistic')
    # axs[1].plot(t, r_predicted_Weibull, 'x-', label='Weibull')
    axs[1].plot(t_fit, pi_predicted_DQQ, 'x-', label=f'DQQ, $R^2$: {R2_pi_DQQ:.3f}', linewidth=8)
    # axs[1].set_xlabel('Time (Days)\n\n(b)\n\n(2) Residential in Washington', fontsize=30, labelpad=20)
    axs[1].set_xlabel('Time (Days)', fontsize=30, labelpad=20)
    axs[1].set_ylabel(f'$\pi(t)$', fontsize=30, labelpad=20)
    axs[1].tick_params(axis='x', labelsize=27.5)  # Set font size for x-axis labels
    # Apply the formatter to the y-axis
    axs[1].yaxis.set_major_formatter(FuncFormatter(three_decimal_formatter))
    axs[1].tick_params(axis='y', labelsize=27.5)  # Set font size for y-axis labels
    # axs[1].set_title(f'Observed vs. Modeled for {state}, {mobility_type}')
    # axs[1].legend(loc='upper right', fontsize=27.5)
    axs[1].grid(True)

    # Save the plot
    filename = f"fig/DQQ_{state}_{mobility_type}_zoomin.png"
    plt.savefig(filename, format='png')

    # plt.tight_layout()
    plt.show()


    states.append(state)
    mobility_types.append(mobility_type)
    # -------------------
    a_V_Shape.append(a_V)
    b_V_Shape.append(b_V)
    r2_pi_V_Shape.append(R2_pi_v_shape)
    r2_pi_V_Shape_disruption.append(R2_pi_v_shape_disruption)
    r2_pi_V_Shape_recovery.append(R2_pi_v_shape_recovery)
    r2_Q_V_Shape.append(R2_Q_v_shape)
    r2_Q_V_Shape_disruption.append(R2_Q_v_shape_disruption)
    r2_Q_V_Shape_recovery.append(R2_Q_v_shape_recovery)
    # -------------------
    A_list.append(A)
    k_list.append(k)
    B_list.append(B)
    A_prime_list.append(A_prime)
    k_prime_list.append(k_prime)
    r2_pi_exp.append(R2_pi_piecewise_exp)
    r2_pi_exp_disruption.append(R2_pi_piecewise_exp_disruption)
    r2_pi_exp_recovery.append(R2_pi_piecewise_exp_recovery)
    r2_Q_exp.append(R2_Q_piecewise_exp)
    r2_Q_exp_disruption.append(R2_Q_piecewise_exp_disruption)
    r2_Q_exp_recovery.append(R2_Q_piecewise_exp_recovery)
    # -------------------
    alpha_DQQ.append(alpha)
    beta_DQQ.append(beta)
    r2_pi_DQQ.append(R2_pi_DQQ)
    r2_pi_DQQ_disruption.append(R2_pi_DQQ_disruption)
    r2_pi_DQQ_recovery.append(R2_pi_DQQ_recovery)
    r2_Q_DQQ.append(R2_Q_DQQ)
    r2_Q_DQQ_disruption.append(R2_Q_DQQ_disruption)
    r2_Q_DQQ_recovery.append(R2_Q_DQQ_recovery)
    # Store values at t2 and t4
    Q_observed_t2_list.append(Q_observed_t2)
    Q_observed_t4_list.append(Q_observed_t4)
    Q_predicted_DQQ_t2_list.append(Q_predicted_DQQ_t2)
    Q_predicted_DQQ_t4_list.append(Q_predicted_DQQ_t4)
    # -------------------
    rmse_Q_V_Shape.append(RMSE_Q_V_Shape)
    rmse_pi_V_Shape.append(RMSE_pi_V_Shape)
    rmse_Q_piecewise_exp.append(RMSE_Q_piecewise_exp)
    rmse_pi_piecewise_exp.append(RMSE_pi_piecewise_exp)
    rmse_Q_DQQ.append(RMSE_Q_DQQ)
    rmse_pi_DQQ.append(RMSE_pi_DQQ)
    # -------------------
    aic_Q_V_Shape.append(AIC_Q_V_Shape)
    aic_pi_V_Shape.append(AIC_pi_V_Shape)
    aic_Q_piecewise_exp.append(AIC_Q_piecewise_exp)
    aic_pi_piecewise_exp.append(AIC_pi_piecewise_exp)
    aic_pi_DQQ.append(AIC_pi_DQQ)
    aic_Q_DQQ.append(AIC_Q_DQQ)

# Compile results into a dataframe
results_df = pd.DataFrame({
    'State': states,
    'Mobility Type': mobility_types,
    'a_V_Shape': a_V_Shape,
    'b_V_Shape': b_V_Shape,
    'R2_pi(t)_v_shape': r2_pi_V_Shape,
    'R2_pi(t)_v_shape_disruption': r2_pi_V_Shape_disruption,
    'R2_pi(t)_v_shape_recovery': r2_pi_V_Shape_recovery,
    'R2_Q(t)_v_shape': r2_Q_V_Shape,
    'R2_Q(t)_v_shape_disruption': r2_Q_V_Shape_disruption,
    'R2_Q(t)_v_shape_recovery': r2_Q_V_Shape_recovery,
    'A_exp': A_list,
    'k_exp': k_list,
    'B_exp': B_list,
    'A_prime_exp': A_prime_list,
    'k_prime_exp': k_prime_list,
    'R2_pi(t)_exp': r2_pi_exp,
    'R2_pi(t)_exp_disruption': r2_pi_exp_disruption,
    'R2_pi(t)_exp_recovery': r2_pi_exp_recovery,
    'R2_Q(t)_exp': r2_Q_exp,
    'R2_Q(t)_exp_disruption': r2_Q_exp_disruption,
    'R2_Q(t)_exp_recovery': r2_Q_exp_recovery,
    'alpha_DQQ': alpha_DQQ,
    'beta_DQQ': beta_DQQ,
    'R2_pi(t)_DQQ': r2_pi_DQQ,
    'R2_pi(t)_DQQ_disruption': r2_pi_DQQ_disruption,
    'R2_pi(t)_DQQ_recovery': r2_pi_DQQ_recovery,
    'R2_Q(t)_DQQ': r2_Q_DQQ,
    'R2_Q(t)_DQQ_disruption': r2_Q_DQQ_disruption,
    'R2_Q(t)_DQQ_recovery': r2_Q_DQQ_recovery,
    'Q_observed_t2': Q_observed_t2_list,
    'Q_observed_t4': Q_observed_t4_list,
    'Q_predicted_DQQ_t2': Q_predicted_DQQ_t2_list,
    'Q_predicted_DQQ_t4': Q_predicted_DQQ_t4_list,
    'RMSE_Q(t)_V_Shape': rmse_Q_V_Shape,
    'RMSE_pi(t)_V_Shape': rmse_pi_V_Shape,
    'RMSE_Q(t)_piecewise_exp': rmse_Q_piecewise_exp,
    'RMSE_pi(t)_piecewise_exp': rmse_pi_piecewise_exp,
    'RMSE_Q(t)_DQQ': rmse_Q_DQQ,
    'RMSE_pi(t)_DQQ': rmse_pi_DQQ,
    'AIC_Q(t)_V_Shape': aic_Q_V_Shape,
    'AIC_pi(t)_V_Shape': aic_pi_V_Shape,
    'AIC_Q(t)_piecewise_exp': aic_Q_piecewise_exp,
    'AIC_pi(t)_piecewise_exp': aic_pi_piecewise_exp,
    'AIC_Q(t)_DQQ': aic_Q_DQQ,
    'AIC_pi(t)_DQQ': aic_pi_DQQ
})

results_df.to_csv("Results/Model_parameters_R2.csv")
