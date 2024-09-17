import pandas as pd
import numpy as np
import matplotlib

matplotlib.use('TkAgg')
matplotlib.rcParams.update({'font.size': 18})
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib.ticker import FuncFormatter
import matplotlib.dates as mdates

fig_path = 'fig'
out_put_path = 'output/output_feb14'


def find_ts(df_for_t, state):
    t0_time_window = 100
    t3_time_window = 120  # First x days chunk to find t1
    t4_time_window = 60  # Time window to get std
    t4_std_threshhold = 1.2

    df_for_t['delta'] = df_for_t[f'{state}_rolling'] - df_for_t[f'{state}_rolling'].shift(1)
    # df_for_t['delta_2'] = df_for_t[f'delta'] - df_for_t[f'delta'].shift(1)
    df_for_t.loc[:, 'delta_2'] = df_for_t['delta'] - df_for_t['delta'].shift(1)
    # df_for_t['delta_3'] = df_for_t[f'delta_2'] - df_for_t[f'delta_2'].shift(1)
    df_for_t.loc[:, 'delta_3'] = df_for_t['delta_2'] - df_for_t['delta_2'].shift(1)

    # Get t0
    df_for_t['sign_change'] = df_for_t[f'{state}_rolling'] * df_for_t[f'{state}_rolling'].shift(1)
    if df_for_t['sign_change'].head(t0_time_window).min() < 0:

        two_day_after_t0 = df_for_t['sign_change'].head(t0_time_window).idxmin()
        t0 = df_for_t[df_for_t['index_num'] == df_for_t.loc[two_day_after_t0]['index_num'] - 2].index[0]
        rate_t0 = df_for_t.loc[t0][f'{state}_rolling']

    else:
        two_day_after_t0 = df_for_t['delta_2'].head(t0_time_window).idxmin()
        t0 = df_for_t[df_for_t['index_num'] == df_for_t.loc[two_day_after_t0]['index_num'] - 2].index[0]
        rate_t0 = df_for_t.loc[t0][f'{state}_rolling']

    # Get t1
    t1 = df_for_t['delta'].rolling(window=7).mean().head(t3_time_window).idxmin()
    rate_t1 = df_for_t.loc[t1][f'{state}_rolling']

    # Get t2
    t2 = df_for_t[f'{state}_rolling'].head(t3_time_window).idxmin()
    rate_t2 = df_for_t[f'{state}_rolling'].head(t3_time_window).min()

    # Get t3

    t3 = df_for_t[(df_for_t['index_num'] > df_for_t.loc[t2]['index_num']) & (df_for_t['index_num'] < 168)][
        'delta'].rolling(window=7).mean().head(t3_time_window).idxmax()
    rate_t3 = df_for_t.loc[t3][f'{state}_rolling']

    # Get t4

    t4 = None

    i = 0
    while (t4 is None):

        if (i >= len(df_for_t) - t4_time_window):
            # t4_time_window -= 10
            t4_std_threshhold += 0.05
            print(t4_std_threshhold)
            if t4_std_threshhold > 20:
                print('LOL')
            i = 0

        date_i = df_for_t.iloc[i].name

        if (date_i.month < 8 and date_i.month > t3.month) or (date_i.month == t3.month and date_i.day > t3.day):
        # if (8 > int(date_i[:2]) > int(t3[:2])) or (int(date_i[:2]) == int(t3[:2]) and int(date_i[-2:]) > int(t3[-2:])):

            std_i = df_for_t[f'{state}_rolling'][i:i + t4_time_window].std()
            # std_list.append(std_i)
            if std_i <= t4_std_threshhold:
                t4 = df_for_t.iloc[i].name
                rate_t4 = df_for_t.loc[t4][f'{state}_rolling']

        i += 1

    T = {
        't0': t0,
        'rate_t0': rate_t0,
        't1': t1,
        'rate_t1': rate_t1,
        't2': t2,
        'rate_t2': rate_t2,
        't3': t3,
        'rate_t3': rate_t3,
        't4': t4,
        'rate_t4': rate_t4,

    }
    return T


classes = ['retail_and_recreation_percent_change_from_baseline',
           'grocery_and_pharmacy_percent_change_from_baseline',
           'parks_percent_change_from_baseline',
           'transit_stations_percent_change_from_baseline',
           'workplaces_percent_change_from_baseline',
           'residential_percent_change_from_baseline']

# Create a mapping from category to its initial
category_initials = {
    'retail_and_recreation_percent_change_from_baseline': 'Retail and Recreation',
    'grocery_and_pharmacy_percent_change_from_baseline': 'Grocery and Pharmacy',
    'parks_percent_change_from_baseline': 'Parks',
    'transit_stations_percent_change_from_baseline': 'Transit Stations',
    'workplaces_percent_change_from_baseline': 'Workplaces',
    'residential_percent_change_from_baseline': 'Residential'
}

markers = [
    's', 'v', '1', '.', '+', 'p'
]

df = pd.read_csv('./data/Region_Mobility_Report_CSVs/US_Region_Mobility_Report.csv',
                 usecols=['sub_region_1', 'date'] + classes)

# df = pd.read_csv('./data/Region_Mobility_Report_CSVs/2020_US_Region_Mobility_Report.csv')

df.columns = ['state', 'date'] + classes

# df['date'] = df.apply(lambda x: x['date'][-10:], axis=1)
df['date'] = pd.to_datetime(df['date'])

df_full = df.dropna(axis=0, how='any')

classes_sequence = ['parks_percent_change_from_baseline',
                    'residential_percent_change_from_baseline',
                    'grocery_and_pharmacy_percent_change_from_baseline',
                    'retail_and_recreation_percent_change_from_baseline',
                    'transit_stations_percent_change_from_baseline',
                    'workplaces_percent_change_from_baseline']

# Define the desired legend labels and their corresponding handles in the correct order
legend_labels = ['Parks', 'Residential', 'Grocery and Pharmacy', 'Retail and Recreation', 'Transit Stations', 'Workplaces']


# for state in state_date.columns:
# for state in ['US', 'Arizona', 'California']:
for state in ['US']:

    plt.figure(figsize=(32, 12))
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.25)  # Adjust the margins

    for i in range(len(classes_sequence)):
        category = classes_sequence[i]
        marker = markers[i]

        # Fetch the initial letter for the category
        initial_letter = category_initials[category]

        state_date = pd.pivot_table(df_full, index='date', columns='state', values=category, aggfunc=np.mean)

        US_avg = []

        for date_pt, date_df in df_full.groupby('date'):
            US_avg.append(date_df[category].mean())

        state_date['US'] = US_avg
        state_date['US_rolling'] = state_date['US'].rolling(window=7).mean()

        state_date[f'{state}_rolling'] = state_date[f'{state}'].rolling(window=7).mean()

        state_date['index_num'] = np.array(range(len(state_date)))

        print()

        if category[0:7] == "transit":

            plt.plot(state_date[f'{state}_rolling'], label=initial_letter, color='red', alpha=0.5, linewidth=6)

            plt.plot(state_date[f'{state}'], color='red', alpha=0.2)

            ts = find_ts(state_date[[f'{state}', f'{state}_rolling', 'index_num']], state)

            # plt.scatter(state_date.loc[ts['t0']]['index_num'], ts['rate_t0'], color='r', s=10)
            # plt.scatter(state_date.loc[ts['t1']]['index_num'], ts['rate_t1'], color='r', s=10)
            # plt.scatter(state_date.loc[ts['t2']]['index_num'], ts['rate_t2'], color='r', s=10)
            # plt.scatter(state_date.loc[ts['t3']]['index_num'], ts['rate_t3'], color='r', s=10)
            # plt.scatter(state_date.loc[ts['t4']]['index_num'], ts['rate_t4'], color='r', s=10)
            #
            # y_lim = [min(-60, state_date[f'{state}'].min() - 10), max(30, state_date[f'{state}'].max() + 80)]
            # text_y = y_lim[0] + 3
            # plt.ylim(y_lim)
            # plt.annotate(xy=(state_date.loc[ts['t0']]['index_num'] - 10, text_y),
            #              text=f'$t_0$:{ts["t0"]}', rotation=90,
            #              bbox=dict(boxstyle='round,pad=0.5', fc='black', lw=1, alpha=0.12), fontsize=27.5, fontname='Arial')
            #
            # plt.annotate(xy=(state_date.loc[ts['t1']]['index_num'] + 2.2, text_y),
            #              text=f'$t_1$:{ts["t1"]}', rotation=90,
            #              bbox=dict(boxstyle='round,pad=0.5', fc='black', lw=1, alpha=0.12), fontsize=27.5, fontname='Arial')
            #
            # plt.annotate(xy=(state_date.loc[ts['t2']]['index_num'] + 2.2, text_y),
            #              text=f'$t_2$:{ts["t2"]}', rotation=90,
            #              bbox=dict(boxstyle='round,pad=0.5', fc='black', lw=1, alpha=0.12), fontsize=27.5, fontname='Arial')
            #
            # plt.annotate(xy=(state_date.loc[ts['t3']]['index_num'] + 2.2, text_y),
            #              text=f'$t_3$:{ts["t3"]}', rotation=90,
            #              bbox=dict(boxstyle='round,pad=0.5', fc='black', lw=1, alpha=0.12), fontsize=27.5, fontname='Arial')
            #
            # plt.annotate(xy=(state_date.loc[ts['t4']]['index_num'] + 2.2, text_y),
            #              text=f'$t_4$:{ts["t4"]}', rotation=90,
            #              bbox=dict(boxstyle='round,pad=0.5', fc='black', lw=1, alpha=0.12), fontsize=27.5, fontname='Arial')
            #
            # # dashline for annotations
            # dashline_x = [state_date.loc[ts['t0']]['index_num'],
            #               state_date.loc[ts['t1']]['index_num'],
            #               state_date.loc[ts['t2']]['index_num'],
            #               state_date.loc[ts['t3']]['index_num'],
            #               state_date.loc[ts['t4']]['index_num']]
            #
            # dashline_y = [ts['rate_t0'], ts['rate_t1'], ts['rate_t2'], ts['rate_t3'], ts['rate_t4']]
            # x_major_locator = MultipleLocator(20)
            # ax = plt.gca()
            #
            # ax.xaxis.set_major_locator(x_major_locator)
            # ax.set_xlabel(state_date.index.tolist())
            # ax.vlines(dashline_x, np.ones(len(dashline_x)) * y_lim[0], ymax=dashline_y,
            #           linestyles='dashed', colors='red', alpha=0.8)

        else:
            plt.plot(state_date[f'{state}_rolling'], label=initial_letter, marker=marker, markersize=10, markevery=10)
            plt.plot(state_date[f'{state}'], color='grey', alpha=0.2)
            x_major_locator = MultipleLocator(20)
            ax = plt.gca()

            ax.xaxis.set_major_locator(x_major_locator)
            ax.set_xlabel(state_date.index.tolist())

    ax = plt.gca()  # Get the current axes instance
    ax.tick_params(axis='y', labelsize=27.5)  # Set font size for y-axis labels
    # Set font type for y-axis tick labels to Arial
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontname('Arial')  # Use label1 instead of label

    # Set the x-axis major formatter to display mm-dd format
    # plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))

    plt.xticks(rotation=60, fontsize=27.5, fontname='Arial')
    plt.legend(loc='upper right', fontsize=27.5, prop={'family': 'Arial', 'size': 27.5})
    plt.xlabel('\n\n\nDate', fontsize=30, fontname='Arial')
    plt.ylabel('Mobility Change Rate (%)', fontsize=30, fontname='Arial')
    plt.savefig(f'{fig_path}/{state}_multi_category.png', dpi=200)
    print(f'Plotting {state}')
    plt.close()
