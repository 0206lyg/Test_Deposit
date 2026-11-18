#!/usr/bin/env python3
"""Compare converged half-channel solutions at one imposed flow and mesh level.

python tools/compare_results.py GLY75_CASE CROSS_CASE EPTT_CASE --output comparison
Add --refined GLY75_FINE CROSS_FINE EPTT_FINE for a two-grid comparison.
--deposit-band LO HI supplies a measured step footprint in µm from the outer corner.
"""
import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from corner_vortex import vortex

FLUIDS = ('gly75', 'cross', 'peo')
WALLS = ('upstream', 'step', 'downstream')
LABELS = dict(gly75='Gly75', cross='PEO Cross', peo='PEO ePTT')
COLORS = dict(gly75='#2864a0', cross='#bc690c', peo='#147761')
STYLES = dict(gly75='-', cross='--', peo='-.')
WALL_COLS = ('iteration', 's_m', 'x_wall_m', 'y_wall_m', 'x_sample_m', 'y_sample_m', 'ds_m',
             'u_toward_m_s', 'tau_total_Pa', 'tau_viscous_Pa', 'tau_mode1_Pa', 'tau_mode2_Pa')
FLOW_COLS = ('iteration', 'x_m', 'y_m', 'dx_m', 'dy_m', 'Ux_m_s', 'Uy_m_s', 'qx_left_m2_s')


def read_table(path, required, iteration):
    data = np.atleast_1d(np.genfromtxt(path, names=True, delimiter=',', dtype=None, encoding='utf-8'))
    names = data.dtype.names or ()
    if not len(data) or not set(required).issubset(names):
        raise ValueError(f'{path}: missing records or columns {set(required)-set(names)}')
    if any(not np.all(np.isfinite(data[k])) for k in required if k != 'wall'):
        raise ValueError(f'{path}: non-finite numerical field')
    if not np.all(data['iteration'] == iteration):
        raise ValueError(f'{path}: iteration does not match the stationary marker ({iteration})')
    return data


def write_csv(path, rows, fieldnames=None):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_case(path, fluid):
    meta = json.loads((path/'metadata.json').read_text())
    out = path/'postProcessing/deposition'
    marker = json.loads((out/'stationary.json').read_text())
    if meta['fluid_key'] != fluid or marker.get('status') != 'stationary':
        raise ValueError(f'{path}: need a converged {fluid} case')
    q = meta['Q_mL_h'] * 1e-6 / (3600 * 2 * meta['geometry_m']['height'])
    if q <= 0 or not math.isclose(q, marker['qHalfPerDepth_m2_s'], rel_tol=1e-7):
        raise ValueError(f'{path}: inlet normalization differs from Q/(2h)')
    iteration = marker['iteration']
    profiles = read_table(out/'wall_profile.csv', ('wall',) + WALL_COLS, iteration)
    walls = {}
    for key in WALLS:
        w = np.sort(profiles[profiles['wall'] == key], order='s_m')
        if len(w) < 2 or np.any(w['ds_m'] <= 0) or np.any(np.diff(w['s_m']) <= 0):
            raise ValueError(f'{path}: missing or invalid {key} wall faces')
        left, right = w['s_m']-w['ds_m']/2, w['s_m']+w['ds_m']/2
        if not np.allclose(left[1:], right[:-1], rtol=1e-8, atol=1e-12):
            raise ValueError(f'{path}: gaps/overlaps in {key} line quadrature')
        if not np.allclose(w['tau_total_Pa'], w['tau_viscous_Pa'] + w['tau_mode1_Pa'] + w['tau_mode2_Pa'],
                           rtol=1e-7, atol=1e-9):
            raise ValueError(f'{path}: viscous and two modal tractions do not sum to total on {key}')
        if not np.allclose(np.hypot(w['x_wall_m']-w['x_sample_m'], w['y_wall_m']-w['y_sample_m']),
                           2e-6, atol=1e-12, rtol=1e-7):
            raise ValueError(f'{path}: measurement offset differs from 2 µm')
        walls[key] = w
    cells = read_table(out/'flow_cells.csv', FLOW_COLS, iteration)
    dictionaries = {}
    for name in ('fvSchemes', 'fvSolution'):
        file = path/'system'/name
        if file.is_file():
            dictionaries[name] = hashlib.sha256(file.read_bytes()).hexdigest()
    scheme = None
    if fluid == 'peo' and 'fvSchemes' in dictionaries:
        schemes = (path/'system/fvSchemes').read_text()
        matches = re.findall(r'div\(phi,theta_mode[12]\)\s+Gauss\s+([^;]+);', schemes)
        scheme = ', '.join(dict.fromkeys(matches))
    return dict(advection=scheme, dictionaries=dictionaries, path=path, fluid=fluid, meta=meta, walls=walls, cells=cells, q=q,
                level=meta['mesh_level'], iteration=iteration)


def cases(paths):
    group = [load_case(path, key) for path, key in zip(paths, FLUIDS)]
    first = group[0]
    for c in group[1:]:
        if (c['meta']['geometry_m'] != first['meta']['geometry_m'] or c['level'] != first['level']
                or not math.isclose(first['q'], c['q'], rel_tol=1e-10)):
            raise ValueError('All three cases must use the same geometry, flow rate and mesh')
        for key in WALLS:
            for coord in ('s_m', 'ds_m'):
                if not np.allclose(first['walls'][key][coord], c['walls'][key][coord], rtol=1e-8, atol=1e-12):
                    raise ValueError(f'Model grids do not match on {key}')
    return group


def stats(c, key, near=False):
    w = c['walls'][key]
    s, ds, u, tau = (w[k] for k in ('s_m', 'ds_m', 'u_toward_m_s', 'tau_total_Pa'))
    left, right = s - .5*ds, s + .5*ds
    lo, hi = left[0], right[-1]
    if near:
        wu = c['meta']['geometry_m']['upstream_width']
        lo, hi = max(lo, -wu), min(hi, wu)
        active = (right > lo) & (left < hi)
        w, s, u, tau = w[active], s[active], u[active], tau[active]
        left, right = np.maximum(left[active], lo), np.minimum(right[active], hi)
        ds = right - left
    positive = u > 0
    starts = np.flatnonzero(positive & ~np.r_[False, positive[:-1]])
    ends = np.flatnonzero(positive & ~np.r_[positive[1:], False])
    runs = []
    for i, j in zip(starts, ends):
        segment = slice(i, j+1)
        peak = i + np.argmax(u[segment])
        f = float(np.sum(u[segment] * ds[segment]))
        runs.append(dict(fluid=c['fluid'], wall=key, start_s_m=left[i], end_s_m=right[j],
                         width_m=ds[segment].sum(), peak_s_m=s[peak], peak_u_m_s=u[peak],
                         Qplus_m2_s=f, Qplus_over_Qhalf=f/c['q'],
                         tau_at_peak_Pa=tau[peak], mean_abs_tau_Pa=np.average(abs(tau[segment]), weights=ds[segment]),
                         viscous_at_peak_Pa=w['tau_viscous_Pa'][peak], mode1_at_peak_Pa=w['tau_mode1_Pa'][peak],
                         mode2_at_peak_Pa=w['tau_mode2_Pa'][peak]))
    full = float(np.sum(np.maximum(u, 0) * ds))
    return dict(fluid=c['fluid'], wall=key, region='within_one_Wu' if near else 'full', s_start_m=lo, s_end_m=hi,
                Qhalf_m2_s=c['q'], Qplus_m2_s=full, Qplus_over_Qhalf=full/c['q'],
                mean_abs_tau_Pa=np.average(abs(tau), weights=ds),
                mean_signed_tau_Pa=np.average(tau, weights=ds)), runs


def ratio(n, d):
    return n/d if d != 0 else None


def profile_comparison(group):
    rows = []
    for key in WALLS:
        profiles = [c['walls'][key] for c in group]
        for gly, cross, peo in zip(*profiles):
            row = dict(wall=key, s_m=gly['s_m'], ds_m=gly['ds_m'],
                       x_wall_m=gly['x_wall_m'], y_wall_m=gly['y_wall_m'],
                       x_sample_m=gly['x_sample_m'], y_sample_m=gly['y_sample_m'])
            for name, f in zip(FLUIDS, (gly, cross, peo)):
                for k in ('u_toward_m_s', 'tau_total_Pa'):
                    row[f'{name}_{k}'] = f[k]
            row.update(peo_viscous_Pa=peo['tau_viscous_Pa'], peo_mode1_Pa=peo['tau_mode1_Pa'],
                       peo_mode2_Pa=peo['tau_mode2_Pa'],
                       abs_cross_over_gly=ratio(abs(cross['tau_total_Pa']), abs(gly['tau_total_Pa'])),
                       abs_eptt_over_gly=ratio(abs(peo['tau_total_Pa']), abs(gly['tau_total_Pa'])))
            rows.append(row)
    return rows


def wall_figure(group, band, band_label, out):
    fig, axes = plt.subplots(2, 2, figsize=(9.3, 6.2), sharex='col', sharey='row',
                             gridspec_kw={'width_ratios': (1.55, 1)}, constrained_layout=True)
    g = group[0]['meta']['geometry_m']
    span = .5e6*(g['upstream_width']-g['downstream_width'])
    zoom_start = max(0., span-10.)
    markers = dict(gly75='o', cross='s', peo='^')
    for c in group:
        w, key = c['walls']['step'], c['fluid']
        for i, (column, scale) in enumerate((('u_toward_m_s', 1e3), ('tau_total_Pa', 1))):
            s = w['s_m']*1e6
            for j, ax in enumerate(axes[i]):
                ax.plot(s, w[column]*scale, STYLES[key], color=COLORS[key], lw=1.8, label=LABELS[key],
                        marker=markers[key] if j else None, markersize=4.3,
                        markerfacecolor='white', markeredgewidth=1.1,
                        markevery=np.flatnonzero(s >= zoom_start).tolist() if j else None)
    for i in range(2):
        # The shading identifies the magnified coordinate range on the left.
        axes[i, 0].axvspan(zoom_start, span, color='.45', alpha=.08, zorder=-2)
        for ax in axes[i]:
            ax.axhline(0, lw=.8, color='.42', zorder=0)
            if band is not None:
                ax.axvspan(*band, color='#8053a2', alpha=.17, label=band_label)
            ax.grid(axis='y', alpha=.22)
    axes[0, 0].set_ylabel('Velocity toward step (mm/s)')
    axes[1, 0].set_ylabel('Signed wall traction (Pa)')
    axes[1, 0].set_xlim(0, span)
    axes[1, 1].set_xlim(zoom_start, span)
    full_ticks = [x for x in np.arange(0., span, 10.) if x < span-.05*span]
    axes[1, 0].set_xticks(full_ticks + [span])
    axes[1, 1].set_xticks(np.linspace(zoom_start, span, 5))
    for ax in axes[1]:
        ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:g}'))
    if band is None and math.isclose(group[0]['meta']['Q_mL_h'], .3):
        # Point to the geometrical landmark on the axis. The outline wraps
        # around the corner; it cannot define a numerical capture footprint.
        axes[1, 0].annotate('Fig. 5(c), 135 min\nDeposit near lip (qualitative)',
                            xy=(1, 0), xycoords='axes fraction', xytext=(.82, .11),
                            textcoords='axes fraction', ha='right', va='bottom', fontsize=8, color='#77355b',
                            bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=2),
                            arrowprops=dict(arrowstyle='->', color='#77355b'))
    axes[0, 0].legend(ncol=2, fontsize=9, loc='upper left')
    axes[0, 0].set_title('Full step', loc='left', fontsize=11)
    axes[0, 1].set_title('Lip detail · symbols at wall-face centres', loc='left', fontsize=11)
    fig.supxlabel('Distance along step from outer corner toward lip, s (µm)', fontsize=11)
    scheme = next((c['advection'] for c in group if c['advection']), None)
    title = f"Step face   |   {group[0]['meta']['Q_mL_h']:g} mL/h   |   {group[0]['level']} mesh"
    fig.suptitle(title + (f"   |   ePTT {scheme}" if scheme else ""), fontsize=14)
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'step_wall_comparison.{ext}', dpi=220)
    plt.close(fig)


def grid_comparison(coarse, fine):
    if (coarse['q'] != fine['q'] or coarse['meta']['geometry_m'] != fine['meta']['geometry_m']
            or coarse['meta']['fluid'] != fine['meta']['fluid'] or coarse['level'] == fine['level']):
        raise ValueError('Mesh comparison requires different grids with matched flow, geometry and material')
    rows = []
    for key in WALLS:
        w, f = coarse['walls'][key], fine['walls'][key]
        if len(f) <= len(w):
            raise ValueError(f'The --refined profile on {key} must have more wall faces')
        # Split at both grids' face edges. This preserves the native midpoint
        # quadratures without smoothing a singularity or inventing endpoint data.
        def edges(a):
            return np.r_[a['s_m'] - .5*a['ds_m'], a['s_m'][-1] + .5*a['ds_m'][-1]]
        ec, ef = edges(w), edges(f)
        if not np.allclose(ec[[0, -1]], ef[[0, -1]], rtol=1e-8, atol=1e-12):
            raise ValueError(f'Different wall extents for {key}')
        merged = np.unique(np.r_[ec, ef])
        lo, hi = max(ec[0], ef[0]), min(ec[-1], ef[-1])
        regions = [('full', lo, hi), ('trimmed_2um', lo+2e-6, hi-2e-6)]
        if key != 'step':
            wu = coarse['meta']['geometry_m']['upstream_width']
            regions.append(('within_one_Wu', max(lo, -wu), min(hi, wu)))
        for region, low, high in regions:
            if high <= low:
                continue
            e = np.r_[low, merged[(merged > low) & (merged < high)], high]
            segment, pos = np.diff(e), .5*(e[:-1]+e[1:])
            cw = w[np.clip(np.searchsorted(ec, pos, side='right')-1, 0, len(w)-1)]
            fw = f[np.clip(np.searchsorted(ef, pos, side='right')-1, 0, len(f)-1)]
            qs = [np.sum(segment*np.maximum(v['u_toward_m_s'], 0)) for v in (cw, fw)]
            means = [np.average(abs(v['tau_total_Pa']), weights=segment) for v in (cw, fw)]
            def l1(k):
                return ratio(np.sum(segment*abs(fw[k]-cw[k])), np.sum(segment*abs(fw[k])))
            rows.append(dict(fluid=coarse['fluid'], wall=key, region=region,
                             base_grid=coarse['level'], refined_grid=fine['level'],
                             Qplus_base_m2_s=qs[0], Qplus_refined_m2_s=qs[1],
                             delta_Qplus_over_Qhalf=(qs[1]-qs[0])/coarse['q'],
                             relative_Qplus_change=ratio(qs[1]-qs[0], qs[1]),
                             mean_abs_tau_base_Pa=means[0], mean_abs_tau_refined_Pa=means[1],
                             relative_mean_abs_tau_change=ratio(means[1]-means[0], means[1]),
                             relative_L1_velocity=l1('u_toward_m_s'), relative_L1_shear=l1('tau_total_Pa')))
    return rows



def developed_cross_stress(q, h, fluid):
    """Planar pressure-driven Cross flow, q/h² = integral_0^Tw t*gammadot(t)dt/Tw².

    Changing from stress to shear rate and integrating by parts avoids an
    inverse Cross-law solve at each quadrature point.
    """
    k, a = fluid['k_s'], fluid['a']
    eta0, etai = fluid['eta_zero_Pa_s'], fluid['eta_inf_Pa_s']
    def stress(rate):
        return rate * (etai + (eta0-etai)/(1 + (k*rate)**a))
    nodes, weights = np.polynomial.legendre.leggauss(64)
    r, wt = .5*(nodes+1), .5*weights
    target = q/h**2
    # q/h² lies between the Newtonian limiting relationships with eta0 and etai.
    lo, hi = 0., 3*target*(eta0/etai)
    for _ in range(70):
        g = .5*(lo+hi)
        flux = .5*g*(1. - np.dot(wt, (stress(r*g)/stress(g))**2))
        if flux < target:
            lo = g
        else:
            hi = g
    return stress(.5*(lo+hi))


def straight_benchmarks(group):
    gly, cross, _ = group
    g, q = gly['meta']['geometry_m'], gly['q']
    rows = []
    for key, pos in [('upstream', -500e-6), ('downstream', 500e-6)]:
        h = g[f'{key}_width']/2
        viscous = 3*gly['meta']['fluid']['eta_Pa_s']*q/h**2
        shear_thinning = developed_cross_stress(q, h, cross['meta']['fluid'])
        wall = gly['walls'][key]
        j = np.searchsorted(wall['s_m'], pos)
        if not 0 < j < len(wall):
            raise ValueError(f'{key}: location x={pos:g} outside wall profile')
        observed = [abs(np.interp(pos, c['walls'][key]['s_m'], c['walls'][key]['tau_total_Pa'])) for c in group]
        measured_ratio = ratio(observed[1], observed[0])
        expected_ratio = shear_thinning/viscous
        rows.append(dict(wall=key, x_m=pos, bracket_left_m=wall['s_m'][j-1], bracket_right_m=wall['s_m'][j],
                         abs_Gly75_Pa=observed[0], abs_Cross_Pa=observed[1], abs_ePTT_Pa=observed[2],
                         developed_Gly75_Pa=viscous, developed_Cross_Pa=shear_thinning,
                         observed_Cross_over_Gly75=measured_ratio, developed_Cross_over_Gly75=expected_ratio,
                         relative_ratio_error=(measured_ratio/expected_ratio-1 if measured_ratio is not None else None),
                         observed_ePTT_over_Gly75=ratio(observed[2], observed[0])))
    return rows


def vortex_figure(group, out):
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.9), sharex=True, sharey=True, constrained_layout=True)
    summaries, contours = [], []
    for c, ax in zip(group, axes):
        g = c['meta']['geometry_m']
        result, line, x, y, psi, speed = vortex(c['cells'], c['walls'], g, c['q'])
        summaries.append(dict(fluid=c['fluid'], **result))
        contours.extend(dict(fluid=c['fluid'], x_m=a, y_m=b) for a, b in line)
        x0 = -max(.55*g['upstream_width'], 1.2*(result['Lv_m'] or 0))
        # Common color range reproduces the slow-flow magnification in Fig. 8(c).
        speed_map = ax.pcolormesh(x*1e6, y*1e6, 1e3*speed.T, shading='flat',
                                  rasterized=True, vmin=0, vmax=2, cmap='cividis')
        stream = psi.T/c['q']
        levels = np.array([.0001, .0005, .002, .01, .03, .08, .16, .3, .55, .8])
        levels = levels[(levels > stream.min()) & (levels < stream.max())]
        if len(levels):
            ax.contour(x*1e6, y*1e6, stream, levels=levels, colors='white', linewidths=.6, alpha=.85)
        if len(line):
            closed_levels = stream.min()*np.array([.8, .45, .15])
            ax.contour(x*1e6, y*1e6, stream, levels=closed_levels, colors='white',
                       linewidths=.7, linestyles='dashed', alpha=.9)
            ax.plot(line[:, 0]*1e6, line[:, 1]*1e6, '-', lw=1.8, color='#ef502e')
        ax.set_title(LABELS[c['fluid']])
        ax.set(xlabel='x (µm)', aspect='equal')
        ax.plot([0, 0], [.5e6*g['downstream_width'], .5e6*g['upstream_width']], color='black', lw=2)
        ax.set_xlim(min(x0*1e6, axes[0].get_xlim()[0]) if ax is not axes[0] else x0*1e6, 1)
        ax.set_ylim(0, .5e6*g['upstream_width'])
    axes[0].set_ylabel('Distance from symmetry plane (µm)')
    fig.colorbar(speed_map, ax=axes, location='bottom', aspect=45, shrink=.67,
                 label='Upstream fluid speed (mm/s)', extend='max')
    fig.savefig(out/'upstream_vortex.png', dpi=200, bbox_inches='tight', pad_inches=.12)
    plt.close(fig)
    write_csv(out/'vortex_metrics.csv', summaries)
    write_csv(out/'separatrix.csv', contours, ['fluid', 'x_m', 'y_m'])


def compare(group, out, band, band_label):
    out.mkdir(parents=True, exist_ok=True)
    wall_figure(group, band, band_label, out)
    vortex_figure(group, out)
    summaries, intervals = [], []
    for c in group:
        for key in WALLS:
            summary, runs = stats(c, key)
            summaries.append(summary)
            if key != 'step':
                summaries.append(stats(c, key, near=True)[0])
            intervals.extend(runs)
    write_csv(out/'wall_integrals.csv', summaries)
    write_csv(out/'approach_intervals.csv', intervals,
              ['fluid', 'wall', 'start_s_m', 'end_s_m', 'width_m', 'peak_s_m', 'peak_u_m_s',
               'Qplus_m2_s', 'Qplus_over_Qhalf', 'tau_at_peak_Pa', 'mean_abs_tau_Pa',
               'viscous_at_peak_Pa', 'mode1_at_peak_Pa', 'mode2_at_peak_Pa'])
    write_csv(out/'wall_profiles.csv', profile_comparison(group))
    write_csv(out/'straight_wall_check.csv', straight_benchmarks(group))
    record = dict(cases=[dict(fluid=c['fluid'], source=str(c['path'].resolve()),
                              iteration=int(c['iteration']), mesh=c['level'],
                              mesh_definition=c['meta'].get('mesh'), material=c['meta']['fluid'],
                              numerics=c['meta'].get('numerics'), dictionaries_sha256=c['dictionaries'],
                              theta_advection=c['advection']) for c in group],
                  experimental_Q_mL_h=group[0]['meta']['Q_mL_h'], step_band_um=band,
                  band_source=band_label if band is not None else None,
                  figure5_note=('Fig. 5(c), Q=0.3 mL/h and 135 min: an outlined deposit adjacent to the throat lip; '
                                'the projected image does not quantify its first attachment coordinate.'),
                  conventions=[
                      's=x on straight walls; on step s=Wu/2-y from outer corner to lip. The normal points out of the fluid.',
                      'Positive u_toward points outward at a 2 µm inward offset. The FV normal velocity interpolates oriented phi/area between parallel faces of the containing cell.',
                      'Signed traction is t dot (T dot n), t=+ex on straight walls and -ey on the step. Ratios use magnitudes.',
                      'Positive intervals start/end at the edges of contiguous positive midpoint samples. Sum max(u,0)*ds per unit depth.',
                      'Straight-wall near-contraction integrals use |x| <= Wu = 180 µm. Whole-wall integrals include the inlet and outlet.',
                      'Straight-duct control checks use interpolated wall stress at x=±500 µm and the exact developed 2D flow relation.',
                      'Mesh changes use (refined-base)/refined; normalized L1 uses the refined profile norm. Zero-denominator ratios are empty.',
                      'All faces included. Ratios near zero Gly75 stress and point values at sharp corners are mesh sensitive.',
                      'The vortex criterion requires a wall-connected negative streamfunction covering at least two x columns and y rows.',
                      'The streamfunction integrates oriented x-face phi/depth from the no-slip wall. The dividing curve follows its interior zero.',
                      'Lv uses the interpolated upstream viscous wall-shear reversal. Velocity means and percentiles weight fluid area.',
                      'An unresolved vortex has empty lengths and speeds. No physical zero is assigned.',
                  ])
    (out/'comparison.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('gly75', type=Path)
    p.add_argument('cross', type=Path)
    p.add_argument('peo', type=Path)
    p.add_argument('--output', type=Path, default=Path('comparison'))
    p.add_argument('--refined', type=Path, nargs=3, metavar=('GLY', 'CROSS', 'EPTT'))
    p.add_argument('--deposit-band', type=float, nargs=2, metavar=('LO_UM', 'HI_UM'),
                   help='User-measured footprint, step distance s=0 at outer corner (µm)')
    args = p.parse_args()
    try:
        group = cases([args.gly75, args.cross, args.peo])
        g = group[0]['meta']['geometry_m']
        band, title = args.deposit_band, 'Measured deposit footprint'
        if band is not None and not 0 <= band[0] < band[1] <= .5e6*(g['upstream_width']-g['downstream_width']) + 1e-8:
            raise ValueError('Deposit band falls outside the step')
        if args.refined:
            refined = cases(args.refined)
            grid = [r for a, b in zip(group, refined) for r in grid_comparison(a, b)]
        compare(group, args.output, band, title)
        if args.refined:
            compare(refined, args.output/'refined', band, title)
            write_csv(args.output/'mesh_comparison.csv', grid)
    except (ValueError, OSError, KeyError) as err:
        p.exit(2, f'Comparison failed: {err}\n')
    print(f"Comparison for Q={group[0]['meta']['Q_mL_h']:g} mL/h: {args.output.resolve()}")


if __name__ == '__main__':
    main()
