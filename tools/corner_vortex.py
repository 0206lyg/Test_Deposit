"""Resolve the primary inlet-corner vortex on a rectilinear half-channel mesh."""
from collections import deque

import numpy as np


def corner_component(inside, first_step_row):
    """Find the largest recirculating region connected to both physical walls."""
    remaining = inside.copy()
    candidates = []
    for i in np.flatnonzero(inside[:, -1]):
        if not remaining[i, -1]:
            continue
        pending = deque([(i, inside.shape[1] - 1)])
        component = []
        while pending:
            a, b = pending.pop()
            if not (0 <= a < inside.shape[0] and 0 <= b < inside.shape[1]) or not remaining[a, b]:
                continue
            remaining[a, b] = False
            component.append((a, b))
            pending.extend(((a - 1, b), (a + 1, b), (a, b - 1), (a, b + 1)))
        if any(a == inside.shape[0] - 1 and b >= first_step_row for a, b in component):
            candidates.append(component)
    return candidates


def reverse_to_forward(s, stress, reverse_end):
    """Interpolate a wall-flow reversal from solvent/viscous shear (outward n)."""
    index = np.flatnonzero((stress[:-1] > 0) & (stress[1:] <= 0))
    if index.size:
        i = index[np.argmin(np.abs(s[index] - reverse_end))]
        return float(s[i] - stress[i] * (s[i + 1] - s[i]) / (stress[i + 1] - stress[i]))
    return None


def vortex(flow, walls, g, qhalf):
    """Integrate phi from the upper no-slip wall; psi<0 inside the primary vortex.

    A left x-face spans one y interval. Its discrete volume flux per unit depth
    is the exact change in streamfunction across that interval. Coordinates
    are clustered only at roundoff scale; the section sum provides an
    independent check against the imposed inlet flow.
    """
    upstream = flow[flow['x_m'] < 0]
    # Adjacent block rows can give the same x with different floating roundoff.
    def axis(values):
        order = np.argsort(values)
        groups = np.r_[0, np.cumsum(np.diff(values[order]) > 1e-12)]
        index = np.empty_like(order)
        index[order] = groups
        coordinates = np.bincount(groups, weights=values[order]) / np.bincount(groups)
        return coordinates, index
    xleft, ix = axis(upstream['x_m'] - upstream['dx_m'] / 2)
    yt, iy = axis(upstream['y_m'] + upstream['dy_m'] / 2)
    nx, ny = len(xleft), len(yt)
    if nx * ny != len(upstream) or np.unique(ix * ny + iy).size != len(upstream):
        raise ValueError('Upstream cells must form a complete, single-layer rectangular grid')
    # From solver: signed conservative flux across the left face of each cell.
    def grid(name):
        v = np.empty((nx, ny))
        v[ix, iy] = upstream[name]
        return v
    qx = grid('qx_left_m2_s')
    speed = np.hypot(grid('Ux_m_s'), grid('Uy_m_s'))
    x = np.r_[xleft, 0.]
    y = np.r_[0., yt]
    if not (np.isclose(y[-1], g['upstream_width']/2, rtol=1e-7, atol=1e-12)
            and np.allclose(np.diff(x)[:, None], grid('dx_m'), rtol=1e-5, atol=1e-12)
            and np.allclose(np.diff(y)[None, :], grid('dy_m'), rtol=1e-5, atol=1e-12)):
        raise ValueError('Unexpected upstream mesh coordinates or cell sizes')

    # Close x=0 using downstream left-face flux below the lip. The vertical
    # solid step carries exactly zero face flux from the lip to the outer wall.
    throat = flow[(flow['x_m'] > 0) & np.isclose(flow['x_m'] - flow['dx_m']/2, 0., atol=1e-12, rtol=0)]
    throat = np.sort(throat, order='y_m')
    core_rows = (yt <= g['downstream_width']/2 + 1e-12)
    if len(throat) != core_rows.sum() or not np.allclose(throat['y_m'] + throat['dy_m']/2, yt[core_rows],
                                                            rtol=1e-6, atol=1e-12):
        raise ValueError('Downstream entrance rows do not match upstream mesh')
    at_step = np.zeros(ny)
    at_step[core_rows] = throat['qx_left_m2_s']
    qx = np.vstack((qx, at_step))
    psi = np.zeros((nx + 1, ny + 1))
    psi[:, :-1] = np.cumsum(qx[:, ::-1], axis=1)[:, ::-1]
    psi_cell = .25 * (psi[:-1, :-1] + psi[:-1, 1:] + psi[1:, :-1] + psi[1:, 1:])
    area = np.diff(x)[:, None] * np.diff(y)[None, :]
    # A flux tolerance at roundoff scale only removes algebraic zero.
    candidates = corner_component(psi_cell < -qhalf * 1e-12, int(core_rows.sum()))
    stats = dict(status='unresolved', Lv_m=None, Lv_over_Wu=None, speed_mean_m_s=None,
                 speed_p95_m_s=None, speed_max_m_s=None, area_m2=None, resolved_cells=0,
                 max_section_flux_error=float(np.max(np.abs(psi[:, 0] - qhalf)) / qhalf))
    boundary = []
    chosen = np.zeros((nx, ny), dtype=bool)
    if candidates:
        component = max(candidates, key=lambda cells: sum(area[i, j] for i, j in cells))
        a, b = np.array(component).T
        if np.unique(a).size >= 2 and np.unique(b).size >= 2:
            chosen[a, b] = True
    if chosen.any():
        # The dividing streamline is the lower psi=0 boundary of the primary
        # region connected to the corner. Locate its transverse zero on each
        # resolved x column, leaving physical walls out of the interior curve.
        for i in np.flatnonzero(np.any(chosen, axis=1)):
            pc = .5 * (psi[i, :] + psi[i+1, :])
            crossing = np.flatnonzero((pc[:-1] >= 0) & (pc[1:] < 0))
            if not crossing.size:
                continue
            bottom = np.flatnonzero(chosen[i])[0]
            j = crossing[np.argmin(abs(crossing - bottom))]
            zero_y = y[j] - pc[j] * (y[j+1] - y[j]) / (pc[j+1] - pc[j])
            boundary.append((.5*(x[i]+x[i+1]), zero_y))
        if len(boundary) >= 2:
            up, step = walls['upstream'], walls['step']
            # On the upper upstream wall t=+x. Reverse the array to search
            # from the recirculating part back toward the upstream inlet.
            sep = reverse_to_forward(up['s_m'][::-1], up['tau_viscous_Pa'][::-1], boundary[0][0])
            if sep is not None:
                boundary.insert(0, (sep, y[-1]))
            attach_s = reverse_to_forward(step['s_m'], step['tau_viscous_Pa'], y[-1]-boundary[-1][1])
            if attach_s is not None:
                boundary.append((0., y[-1] - attach_s))
            v, weights = speed[chosen], area[chosen]
            order = np.argsort(v)
            p95 = v[order][np.searchsorted(np.cumsum(weights[order]), .95*weights.sum())]
            stats.update(status='resolved' if sep is not None else 'resolved_without_wall_shear_bracket',
                         Lv_m=None if sep is None else -sep,
                         Lv_over_Wu=None if sep is None else -sep/g['upstream_width'],
                         speed_mean_m_s=float(np.average(v, weights=weights)),
                         speed_p95_m_s=float(p95), speed_max_m_s=float(v.max()),
                         area_m2=float(weights.sum()), resolved_cells=int(chosen.sum()))
    return stats, np.asarray(boundary).reshape(-1, 2), x, y, psi, speed
