import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt


def get_range(df, c, comparison):
    """Get genomic position ranges where a denoised signal crosses a threshold."""
    if comparison == '<':
        datac = df.index[df['r'] < c].tolist()
    else:
        datac = df.index[df['r'] > c].tolist()
    if not datac:
        return []

    result = [[datac[0]]]
    for i in range(1, len(datac)):
        if datac[i - 1] + 1 == datac[i]:
            result[-1].append(datac[i])
        else:
            result.append([datac[i]])

    fin_range = [[df.loc[i[0], 'loc'], df.loc[i[-1], 'loc']] for i in result]
    return fin_range


def get_crange(df, c=-0.4):
    return get_range(df, c, '<')


def get_hrange(df, c=0.6):
    return get_range(df, c, '>')


def wavelet_transform(df, value_col, a, hc=0.6, lc=-0.4):
    """Apply wavelet denoising to one PMR signal column from pmr_data."""
    # PyWavelets may receive a read-only view from pandas in some environments.
    rec = np.array(df[value_col].astype(float).to_numpy(), dtype=float, copy=True)
    signal_len = len(rec)

    if a == '30':
        wavefunc = 'db4'
        level = 5
        threshold = 0.8
    else:
        wavefunc = 'db2'
        level = 8
        threshold = 0.8

    max_level = pywt.dwt_max_level(signal_len, pywt.Wavelet(wavefunc).dec_len)
    if signal_len == 0 or max_level < 1:
        y_denoised = rec.copy()
        y_re = pd.DataFrame({
            'r': y_denoised,
            'loc': df['pos'].to_numpy(),
        })
        y_re.index = range(signal_len)
        if a == '30':
            selected_positions = y_re[y_re['r'] > hc]['loc'].tolist()
        else:
            selected_positions = y_re[y_re['r'] < lc]['loc'].tolist()
        return rec, y_denoised, y_re, selected_positions

    level = min(level, max_level)

    coeffs = pywt.wavedec(rec, wavefunc, mode='periodic', level=level)
    for i in range(1, len(coeffs)):
        coeff = coeffs[i]
        max_coeff = np.max(np.abs(coeff)) if len(coeff) else 0
        thr = threshold * max_coeff
        coeffs[i] = coeff if thr == 0 else pywt.threshold(coeff, thr)

    y_denoised = pywt.waverec(coeffs, wavefunc, mode='periodic')[:signal_len]

    y_re = pd.DataFrame({
        'r': y_denoised,
        'loc': df['pos'].to_numpy(),
    })
    y_re.index = range(signal_len)

    if a == '30':
        selected_positions = y_re[y_re['r'] > hc]['loc'].tolist()
    else:
        selected_positions = y_re[y_re['r'] < lc]['loc'].tolist()

    return rec, y_denoised, y_re, selected_positions


def plot_and_save(get_crange, get_hrange, wavelet_transform, result_path, pmr_data, widthf=6 / 2.54, hc=0.6, lc=-0.4):
    """Plot denoised PMR signals and integrate hotspot/coldspot annotations."""
    params = {
        'font.size': 6,
        'axes.spines.right': False,
        'axes.spines.top': False,
        'axes.titlesize': 10,
        'axes.labelsize': 8,
        'legend.loc': 'upper right',
        'legend.fontsize': 8,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'pdf.fonttype': 42,
    }
    mpl.rcParams.update(params)

    df = pmr_data.copy().sort_values('pos').reset_index(drop=True)
    x = df['pos'].to_numpy()

    denoised_30 = None
    denoised_all = None
    hotspots = []
    coldspots = []

    configs = {
        '30': {'col': '30_pmr', 'threshold': hc},
        'all': {'col': 'all_pmr', 'threshold': lc},
    }

    for a in ['30', 'all']:
        rec, y_denoised, y_re, selected_positions = wavelet_transform(
            df=df,
            value_col=configs[a]['col'],
            a=a,
            hc=hc,
            lc=lc,
        )

        if a == '30':
            denoised_30 = y_denoised
            hotspots = selected_positions
        else:
            denoised_all = y_denoised
            coldspots = selected_positions

        plt.figure().set_size_inches(widthf, widthf * 0.618)
        plt.ylim(-1, 1)
        plt.tick_params(which='major', width=0.5, length=3)
        plt.scatter(
            x,
            rec,
            label='original signal',
            s=0.3,
            c='none',
            alpha=0.1,
            marker='o',
            edgecolors='grey',
            linewidths=0.3,
        )
        plt.plot(x, y_denoised, label='Denoised signal', linewidth=0.5)
        plt.xlabel('Genome Position (MB)')
        plt.ylabel('R value')

        if a == '30':
            hotspot_ranges = get_hrange(y_re, hc)
            if hotspot_ranges == []:
                print('no hot region')
            for n, i in enumerate(hotspot_ranges):
                plt.fill_between(
                    [i[0], i[1]],
                    0.95,
                    1,
                    facecolor='#FF0000',
                    label='hot spots' if n == 0 else None,
                )
            plt.axhline(y=hc, c='r', linewidth=0.5)
        else:
            coldspot_ranges = get_crange(y_re, lc)
            if coldspot_ranges == []:
                print('no cold region')
            for n, i in enumerate(coldspot_ranges):
                plt.fill_between(
                    [i[0], i[1]],
                    -1,
                    -0.95,
                    facecolor='#0000FF',
                    label='cold spots' if n == 0 else None,
                )
            plt.axhline(y=lc, c='r', linewidth=0.5)

        plt.savefig(f'{result_path}.fa.pmr.{a}.2.wt.l.pdf', dpi=300, bbox_inches='tight')
        plt.close()

    integrate_outputs(df, denoised_30, denoised_all, hotspots, coldspots, result_path)
    return denoised_30, denoised_all, hotspots, coldspots


def integrate_outputs(pmr_data, denoised_30, denoised_all, hotspots, coldspots, output_prefix):
    """Integrate raw PMR, denoised PMR and hot/cold labels into one TSV."""
    pmr_data = pmr_data.copy().sort_values('pos').reset_index(drop=True)

    if len(denoised_30) != len(pmr_data) or len(denoised_all) != len(pmr_data):
        raise ValueError('Denoised signal length does not match PMR table length.')

    hotspot_set = set(pmr_data.loc[pmr_data['pos'].isin(hotspots), 'pos'].tolist())
    coldspot_set = set(pmr_data.loc[pmr_data['pos'].isin(coldspots), 'pos'].tolist())

    pmr_data['30_pmr_wt'] = np.asarray(denoised_30, dtype=float)
    pmr_data['all_pmr_wt'] = np.asarray(denoised_all, dtype=float)
    pmr_data['is_hot'] = pmr_data['pos'].isin(hotspot_set).astype(int)
    pmr_data['is_cold'] = pmr_data['pos'].isin(coldspot_set).astype(int)

    order = ['pos', 'us_pmr', '30_pmr', '30_pmr_wt', 'is_hot', 'all_pmr', 'all_pmr_wt', 'is_cold']
    pmr_data = pmr_data[order]

    output_file = f'{output_prefix}.integrated.tsv'
    pmr_data.to_csv(output_file, sep='\t', index=False)
    print(f'Integrated output saved to {output_file}')
