#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT RH / Xi-QED Shared Law — Master Evidence Audit
==================================================

GitHub-ready consolidated audit for the frozen Xi-QED shared-law candidate.

This script deliberately DISABLES formula discovery:
  * no coefficient menu search
  * no symbolic grammar search
  * no additive mixed-mass external Q-scale search
  * no replacement of the frozen main law

Frozen law tested:

    K2    = - d^2/dgamma^2 log Xi(1/2 + i gamma) |_{gamma=0}
    A0    = 2*pi/K2
    gap   = alpha(0)^(-1) - A0

    Q1    = (4/pi) * m_pi                 with electron loop
    Q2    = sqrt(2*pi) * red(m_mu,m_pi)   with electron + muon loops

where red(a,b) = ab/(a+b).

Main evidence blocks:
  1. frozen-law baseline score
  2. quadrature stability
  3. random-delta controls
  4. pseudo-mass controls
  5. local mass-ratio grid controls
  6. scheme robustness
  7. threshold-structure mechanism decomposition
  8. threshold stationarity/sensitivity diagnostic
  9. mass-ratio origin diagnostic
 10. final evidence table and global flag

Outputs:
  <prefix>_master_protocol.json
  <prefix>_master_summary.json
  <prefix>_master_meta.json
  <prefix>_master_evidence_table.csv
  plus CSV files for individual evidence blocks.

Author: Buddhika Weerasooriya / ROT-RH audit branch
License suggestion: MIT, Apache-2.0, or CC-BY-4.0 for data/results.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import statistics
import sys
import time
from dataclasses import dataclass, asdict
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

import mpmath as mp


# =============================================================================
# Utility / output
# =============================================================================

LINE = "=" * 120
SUBLINE = "-" * 120


def banner(title: str) -> None:
    print(LINE)
    print(title)
    print(LINE)


def section(title: str) -> None:
    print(LINE)
    print(title)
    print(LINE)


def fmt(x, digits: int = 16) -> str:
    try:
        if isinstance(x, mp.mpf):
            return mp.nstr(x, digits)
        if isinstance(x, float):
            return f"{x:.{digits}g}"
        return str(x)
    except Exception:
        return str(x)


def ffloat(x) -> float:
    return float(x)


def write_csv(path: str, rows: List[Dict[str, object]], fieldnames: Sequence[str] | None = None) -> None:
    if fieldnames is None:
        keys: List[str] = []
        for row in rows:
            for k in row.keys():
                if k not in keys:
                    keys.append(k)
        fieldnames = keys
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: str, obj: object) -> None:
    def conv(o):
        if isinstance(o, mp.mpf):
            return mp.nstr(o, 50)
        if isinstance(o, (float, int, str, bool)) or o is None:
            return o
        if isinstance(o, dict):
            return {str(k): conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [conv(v) for v in o]
        return str(o)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(conv(obj), fh, indent=2, sort_keys=True)


def sha256_json(obj: object) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def empirical_p(count_leq: int, n: int) -> float:
    # add-one correction to avoid zero p-values
    return (count_leq + 1.0) / (n + 1.0)


# =============================================================================
# Xi curvature invariant
# =============================================================================


def xi_completed(s: mp.mpc) -> mp.mpc:
    # Riemann xi: xi(s) = 1/2 s(s-1) pi^{-s/2} Gamma(s/2) zeta(s)
    return mp.mpf("0.5") * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def log_xi_gamma(gamma: mp.mpf) -> mp.mpf:
    # Along critical line: Xi(gamma) = xi(1/2 + i gamma) is real and positive near 0.
    val = xi_completed(mp.mpf("0.5") + 1j * gamma)
    return mp.log(mp.re(val))


def compute_xi_invariants(dps: int) -> Dict[str, mp.mpf]:
    mp.mp.dps = dps
    # Use numerical differentiation at gamma=0. K2 = -L''(0).
    K2 = -mp.diff(log_xi_gamma, mp.mpf("0"), 2)
    L4 = mp.diff(log_xi_gamma, mp.mpf("0"), 4)
    A0 = 2 * mp.pi / K2
    R4 = L4 / (K2 * K2)
    return {"K2": K2, "A0": A0, "L4": L4, "R4": R4}


# =============================================================================
# QED vacuum polarization kernels
# =============================================================================


def gauss_legendre_nodes(n: int, dps: int) -> Tuple[List[mp.mpf], List[mp.mpf]]:
    # mpmath returns nodes/weights on [-1,1]
    mp.mp.dps = dps
    nodes, weights = mp.gauss_quadrature(n, "legendre")
    xs = [(x + 1) / 2 for x in nodes]
    ws = [w / 2 for w in weights]
    return xs, ws


def pi_spacelike_one_loop_alpha_inv_shift(Q: mp.mpf, m: mp.mpf, alpha_inv0: mp.mpf, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf]) -> mp.mpf:
    """Exact on-shell spacelike one-loop contribution to alpha^{-1} shift.

    alpha(Q)^(-1) = alpha(0)^(-1) - Delta(Q)
    Delta_l(Q) = (2/pi) int_0^1 dx x(1-x) log(1 + Q^2 x(1-x)/m_l^2)

    This quantity is independent of alpha_inv0 at one-loop in inverse-coupling units.
    alpha_inv0 is included for protocol clarity and future extension.
    """
    z = (Q / m) ** 2
    total = mp.mpf("0")
    for x, w in zip(xs, ws):
        xx = x * (1 - x)
        total += w * xx * mp.log(1 + z * xx)
    return (2 / mp.pi) * total


def pi_asymp_minus_5over3(Q: mp.mpf, m: mp.mpf) -> mp.mpf:
    return (mp.mpf("1") / (3 * mp.pi)) * (mp.log((Q / m) ** 2) - mp.mpf("5") / 3)


def pi_log_only(Q: mp.mpf, m: mp.mpf) -> mp.mpf:
    return (mp.mpf("1") / (3 * mp.pi)) * mp.log((Q / m) ** 2)


def pi_small_q_series_diag(Q: mp.mpf, m: mp.mpf) -> mp.mpf:
    # Leading small-Q term for spacelike kernel; diagnostic only.
    return (mp.mpf("1") / (15 * mp.pi)) * (Q / m) ** 2


def pi_kernel_component(kind: str, Q: mp.mpf, m: mp.mpf, alpha_inv0: mp.mpf, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf]) -> mp.mpf:
    if kind == "exact_spacelike":
        return pi_spacelike_one_loop_alpha_inv_shift(Q, m, alpha_inv0, xs, ws)
    if kind == "asymptotic_minus_5over3":
        return pi_asymp_minus_5over3(Q, m)
    if kind == "asymptotic_log_only":
        return pi_log_only(Q, m)
    if kind == "small_q_series":
        return pi_small_q_series_diag(Q, m)
    raise ValueError(f"unknown scheme: {kind}")


def delta_for_loops(Q: mp.mpf, masses: Dict[str, mp.mpf], loops: Sequence[str], scheme: str, alpha_inv0: mp.mpf, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf]) -> mp.mpf:
    return mp.fsum(pi_kernel_component(scheme, Q, masses[name], alpha_inv0, xs, ws) for name in loops)


def reduced(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    return (a * b) / (a + b)


# =============================================================================
# Scoring
# =============================================================================


@dataclass
class LawScore:
    Q1: float
    Q2: float
    delta1: float
    delta2: float
    rel1: float
    rel2: float
    joint_max: float
    joint_rms: float
    joint_mean: float
    agreement_rel: float


def relerr(pred: mp.mpf, target: mp.mpf) -> mp.mpf:
    return abs(pred - target) / abs(target)


def score_two_channel(delta_gap: mp.mpf, Q1: mp.mpf, Q2: mp.mpf, masses: Dict[str, mp.mpf], alpha_inv0: mp.mpf, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf], scheme: str = "exact_spacelike", q1_loops: Sequence[str] = ("e",), q2_loops: Sequence[str] = ("e", "mu")) -> LawScore:
    d1 = delta_for_loops(Q1, masses, q1_loops, scheme, alpha_inv0, xs, ws)
    d2 = delta_for_loops(Q2, masses, q2_loops, scheme, alpha_inv0, xs, ws)
    r1 = relerr(d1, delta_gap)
    r2 = relerr(d2, delta_gap)
    return LawScore(
        Q1=ffloat(Q1),
        Q2=ffloat(Q2),
        delta1=ffloat(d1),
        delta2=ffloat(d2),
        rel1=ffloat(r1),
        rel2=ffloat(r2),
        joint_max=ffloat(max(r1, r2)),
        joint_rms=ffloat(mp.sqrt((r1 * r1 + r2 * r2) / 2)),
        joint_mean=ffloat((r1 + r2) / 2),
        agreement_rel=ffloat(abs(d1 - d2) / abs(delta_gap)),
    )


def frozen_Qs(masses: Dict[str, mp.mpf]) -> Tuple[mp.mpf, mp.mpf]:
    Q1 = (4 / mp.pi) * masses["pi"]
    Q2 = mp.sqrt(2 * mp.pi) * reduced(masses["mu"], masses["pi"])
    return Q1, Q2


# =============================================================================
# Controls
# =============================================================================


def random_delta_controls(target_gap: mp.mpf, d1: mp.mpf, d2: mp.mpf, n: int, band: Tuple[float, float], seed: int) -> Dict[str, object]:
    rng = random.Random(seed)
    r1 = float(abs(d1 - target_gap) / abs(target_gap))
    r2 = float(abs(d2 - target_gap) / abs(target_gap))
    jm = max(r1, r2)
    jr = math.sqrt((r1 * r1 + r2 * r2) / 2)
    both_count = 0
    jm_count = 0
    jr_count = 0
    rows: List[Dict[str, object]] = []
    lo, hi = band
    for i in range(n):
        fake = float(target_gap) * rng.uniform(lo, hi)
        fr1 = abs(float(d1) - fake) / abs(fake)
        fr2 = abs(float(d2) - fake) / abs(fake)
        fjm = max(fr1, fr2)
        fjr = math.sqrt((fr1 * fr1 + fr2 * fr2) / 2)
        if fr1 <= r1 and fr2 <= r2:
            both_count += 1
        if fjm <= jm:
            jm_count += 1
        if fjr <= jr:
            jr_count += 1
        rows.append({"i": i, "fake_delta": fake, "rel1": fr1, "rel2": fr2, "joint_max": fjm, "joint_rms": fjr})
    return {
        "rows": rows,
        "both_count": both_count,
        "jointmax_count": jm_count,
        "jointrms_count": jr_count,
        "both_p": empirical_p(both_count, n),
        "jointmax_p": empirical_p(jm_count, n),
        "jointrms_p": empirical_p(jr_count, n),
    }


def pseudo_masses(base: Dict[str, mp.mpf], family: str, jitter: float, rng: random.Random) -> Tuple[Dict[str, mp.mpf], bool]:
    """Generate pseudo masses. Returns (masses, diagnostic_flag).

    Main families break relevant ratios. Diagnostic families preserve too much structure.
    """
    me, mmu, mpi, mtau = [float(base[k]) for k in ("e", "mu", "pi", "tau")]
    diag = False

    def lnj():
        return math.exp(rng.uniform(-jitter, jitter))

    if family == "ratio_jitter_e_pi_only":
        me *= lnj(); mpi *= lnj()
    elif family == "ratio_jitter_pi_only":
        mpi *= lnj()
    elif family == "ratio_jitter_all":
        me *= lnj(); mmu *= lnj(); mpi *= lnj(); mtau *= lnj()
    elif family == "preserve_mu_pi_break_epi":
        # preserve mu/pi ratio but move electron/pion; still useful main control against e/pi accidental matching
        scale = lnj(); mmu *= scale; mpi *= scale; me *= lnj()
    elif family == "random_ratio_ordered_wide":
        # Draw ordered two-body ratios broadly around current scale.
        mmu = mmu
        rho = (mpi / mmu) * math.exp(rng.uniform(-2.0 * jitter, 2.0 * jitter))
        eps = (me / mmu) * math.exp(rng.uniform(-2.0 * jitter, 2.0 * jitter))
        mpi = mmu * rho
        me = mmu * eps
        mtau *= lnj()
    elif family == "affine_plus_small_ratio_noise":
        # Mostly common scaling plus small relative ratio noise. Main but weaker.
        scale = lnj(); me *= scale; mmu *= scale; mpi *= scale; mtau *= scale
        me *= math.exp(rng.uniform(-0.25*jitter, 0.25*jitter))
        mpi *= math.exp(rng.uniform(-0.25*jitter, 0.25*jitter))
    elif family == "preserve_e_pi_break_mupi":
        diag = True
        scale = lnj(); me *= scale; mpi *= scale; mmu *= lnj()
    elif family == "affine_scale_all":
        diag = True
        scale = lnj(); me *= scale; mmu *= scale; mpi *= scale; mtau *= scale
    else:
        raise ValueError(f"unknown pseudo-mass family: {family}")

    return {"e": mp.mpf(me), "mu": mp.mpf(mmu), "pi": mp.mpf(mpi), "tau": mp.mpf(mtau)}, diag


def pseudo_mass_controls(base_masses: Dict[str, mp.mpf], alpha_inv0: mp.mpf, delta_gap: mp.mpf, target_score: LawScore, families: Sequence[str], trials_per_family: int, jitter: float, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf], seed: int) -> Dict[str, object]:
    rng = random.Random(seed)
    rows: List[Dict[str, object]] = []
    fam_summary: List[Dict[str, object]] = []
    main_rows: List[Dict[str, object]] = []
    diag_rows: List[Dict[str, object]] = []

    for fam in families:
        fam_rows = []
        for i in range(trials_per_family):
            pm, diag = pseudo_masses(base_masses, fam, jitter, rng)
            Q1, Q2 = frozen_Qs(pm)
            sc = score_two_channel(delta_gap, Q1, Q2, pm, alpha_inv0, xs, ws)
            row = {
                "family": fam, "trial": i, "diagnostic": diag,
                "me": ffloat(pm["e"]), "mmu": ffloat(pm["mu"]), "mpi": ffloat(pm["pi"]),
                "rho_pi_mu": ffloat(pm["pi"] / pm["mu"]),
                **asdict(sc),
            }
            rows.append(row); fam_rows.append(row)
            (diag_rows if diag else main_rows).append(row)

        both = sum(1 for r in fam_rows if r["rel1"] <= target_score.rel1 and r["rel2"] <= target_score.rel2)
        jm = sum(1 for r in fam_rows if r["joint_max"] <= target_score.joint_max)
        fam_summary.append({
            "family": fam,
            "diagnostic": bool(fam_rows[0]["diagnostic"]) if fam_rows else False,
            "rows": len(fam_rows),
            "both_p": empirical_p(both, len(fam_rows)),
            "jointmax_p": empirical_p(jm, len(fam_rows)),
            "median_joint": statistics.median([r["joint_max"] for r in fam_rows]) if fam_rows else None,
        })

    def counts_for(rs):
        n = len(rs)
        both = sum(1 for r in rs if r["rel1"] <= target_score.rel1 and r["rel2"] <= target_score.rel2)
        jm = sum(1 for r in rs if r["joint_max"] <= target_score.joint_max)
        jr = sum(1 for r in rs if r["joint_rms"] <= target_score.joint_rms)
        return {"rows": n, "both_count": both, "jointmax_count": jm, "jointrms_count": jr,
                "both_p": empirical_p(both, n) if n else None,
                "jointmax_p": empirical_p(jm, n) if n else None,
                "jointrms_p": empirical_p(jr, n) if n else None}

    return {"rows": rows, "family_summary": fam_summary, "main": counts_for(main_rows), "diagnostic": counts_for(diag_rows)}


def local_mass_ratio_grid(base_masses: Dict[str, mp.mpf], alpha_inv0: mp.mpf, delta_gap: mp.mpf, target_score: LawScore, span: float, steps: int, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf]) -> Dict[str, object]:
    if steps < 3 or steps % 2 == 0:
        raise ValueError("local grid steps must be odd and >=3")
    half = steps // 2
    rows: List[Dict[str, object]] = []
    best = None
    for ia in range(steps):
        de = -span + 2 * span * ia / (steps - 1)
        for ib in range(steps):
            dpi = -span + 2 * span * ib / (steps - 1)
            pm = dict(base_masses)
            # vary e/mu and pi/mu while keeping mu fixed
            pm["e"] = base_masses["e"] * mp.e ** de
            pm["pi"] = base_masses["pi"] * mp.e ** dpi
            Q1, Q2 = frozen_Qs(pm)
            sc = score_two_channel(delta_gap, Q1, Q2, pm, alpha_inv0, xs, ws)
            row = {"dlog_e_mu": de, "dlog_pi_mu": dpi, "is_observed": ia == half and ib == half, **asdict(sc)}
            rows.append(row)
            if best is None or row["joint_max"] < best["joint_max"]:
                best = row
    n = len(rows)
    jm_count = sum(1 for r in rows if r["joint_max"] <= target_score.joint_max)
    both_count = sum(1 for r in rows if r["rel1"] <= target_score.rel1 and r["rel2"] <= target_score.rel2)
    offdiag = [r for r in rows if not r["is_observed"]]
    offdiag_count = sum(1 for r in offdiag if r["joint_max"] <= target_score.joint_max)
    return {"rows": rows, "best": best, "jointmax_p": empirical_p(jm_count, n), "both_p": empirical_p(both_count, n), "offdiag_joint_p": empirical_p(offdiag_count, len(offdiag))}


# =============================================================================
# Mechanism blocks
# =============================================================================


def scheme_robustness(base_masses: Dict[str, mp.mpf], alpha_inv0: mp.mpf, delta_gap: mp.mpf, Q1: mp.mpf, Q2: mp.mpf, target_score: LawScore, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf]) -> Dict[str, object]:
    schemes = ["exact_spacelike", "asymptotic_minus_5over3", "asymptotic_log_only", "small_q_series"]
    loop_variants = [
        (("e",), ("e", "mu"), "baseline"),
        (("e",), ("e", "mu", "tau"), "q2_tau_added"),
        (("e",), ("e",), "q2_e_only"),
        (("e", "mu"), ("e", "mu"), "q1_mu_added"),
    ]
    rows: List[Dict[str, object]] = []
    for scheme in schemes:
        for q1_loops, q2_loops, label in loop_variants:
            sc = score_two_channel(delta_gap, Q1, Q2, base_masses, alpha_inv0, xs, ws, scheme, q1_loops, q2_loops)
            rows.append({"scheme": scheme, "loop_label": label, "q1_loops": "+".join(q1_loops), "q2_loops": "+".join(q2_loops), **asdict(sc)})
    physical = [r for r in rows if r["scheme"] != "small_q_series"]
    exact_rows = [r for r in rows if r["scheme"] == "exact_spacelike"]
    exact_sorted = sorted(exact_rows, key=lambda r: r["joint_max"])
    baseline = next(r for r in rows if r["scheme"] == "exact_spacelike" and r["loop_label"] == "baseline")
    baseline_rank = 1 + [r["loop_label"] for r in exact_sorted].index("baseline")
    asymp = next(r for r in rows if r["scheme"] == "asymptotic_minus_5over3" and r["loop_label"] == "baseline")
    logonly = next(r for r in rows if r["scheme"] == "asymptotic_log_only" and r["loop_label"] == "baseline")
    return {
        "rows": sorted(rows, key=lambda r: (r["scheme"] == "small_q_series", r["joint_max"])),
        "baseline_exact_rank": baseline_rank,
        "asymp_ratio": asymp["joint_max"] / target_score.joint_max,
        "logonly_ratio": logonly["joint_max"] / target_score.joint_max,
        "baseline_exact_best": baseline_rank == 1,
    }


def threshold_structure(base_masses: Dict[str, mp.mpf], alpha_inv0: mp.mpf, delta_gap: mp.mpf, Q1: mp.mpf, Q2: mp.mpf, target_score: LawScore, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf]) -> Dict[str, object]:
    comps: List[Dict[str, object]] = []
    for label, Q, lepton in [("Q1", Q1, "e"), ("Q2", Q2, "e"), ("Q2", Q2, "mu")]:
        m = base_masses[lepton]
        exact = pi_kernel_component("exact_spacelike", Q, m, alpha_inv0, xs, ws)
        logv = pi_kernel_component("asymptotic_log_only", Q, m, alpha_inv0, xs, ws)
        asymp = pi_kernel_component("asymptotic_minus_5over3", Q, m, alpha_inv0, xs, ws)
        resid = exact - asymp
        comps.append({
            "channel": label, "lepton": lepton, "Q_over_m": ffloat(Q/m), "two_m_over_Q": ffloat(2*m/Q),
            "exact": ffloat(exact), "log_only": ffloat(logv), "asymp_minus_5over3": ffloat(asymp),
            "finite_mass_residual": ffloat(resid), "residual_over_exact": ffloat(resid/exact) if exact != 0 else None,
        })

    # hybrid ablations for muon component only; electron exact/asymp/log variants
    exact_e_Q1 = pi_kernel_component("exact_spacelike", Q1, base_masses["e"], alpha_inv0, xs, ws)
    exact_e_Q2 = pi_kernel_component("exact_spacelike", Q2, base_masses["e"], alpha_inv0, xs, ws)
    mu_exact = pi_kernel_component("exact_spacelike", Q2, base_masses["mu"], alpha_inv0, xs, ws)
    mu_asymp = pi_kernel_component("asymptotic_minus_5over3", Q2, base_masses["mu"], alpha_inv0, xs, ws)
    mu_log = pi_kernel_component("asymptotic_log_only", Q2, base_masses["mu"], alpha_inv0, xs, ws)
    mu_resid = mu_exact - mu_asymp

    def score_manual(d1: mp.mpf, d2: mp.mpf, label: str) -> Dict[str, object]:
        r1, r2 = relerr(d1, delta_gap), relerr(d2, delta_gap)
        return {"hybrid": label, "delta1": ffloat(d1), "delta2": ffloat(d2), "rel1": ffloat(r1), "rel2": ffloat(r2),
                "joint_max": ffloat(max(r1, r2)), "joint_rms": ffloat(mp.sqrt((r1*r1+r2*r2)/2)), "agreement_rel": ffloat(abs(d1-d2)/abs(delta_gap))}

    hybrids = [
        score_manual(exact_e_Q1, exact_e_Q2 + mu_exact, "baseline_all_exact"),
        score_manual(exact_e_Q1, exact_e_Q2 + mu_asymp, "electron_exact_mu_asymp"),
        score_manual(exact_e_Q1, exact_e_Q2 + mu_log, "electron_exact_mu_log"),
        score_manual(exact_e_Q1, exact_e_Q2 + mu_resid, "electron_exact_mu_residual_only"),
    ]
    best_hybrid = min(hybrids, key=lambda r: r["joint_max"])
    q2_rel_mu_asymp = relerr(exact_e_Q2 + mu_asymp, delta_gap)
    q2_rel_mu_log = relerr(exact_e_Q2 + mu_log, delta_gap)
    return {
        "components": comps,
        "hybrids": sorted(hybrids, key=lambda r: r["joint_max"]),
        "best_hybrid": best_hybrid["hybrid"],
        "q2_mu_exact": ffloat(mu_exact),
        "q2_mu_asymp": ffloat(mu_asymp),
        "q2_mu_log": ffloat(mu_log),
        "q2_mu_residual": ffloat(mu_resid),
        "q2_mu_residual_over_exact": ffloat(mu_resid / mu_exact),
        "q2_exact_to_mu_asymp_degradation": ffloat(q2_rel_mu_asymp / target_score.rel2),
        "q2_exact_to_mu_log_degradation": ffloat(q2_rel_mu_log / target_score.rel2),
    }


def derivative_log_r(func: Callable[[mp.mpf], mp.mpf], r: mp.mpf, order: int) -> mp.mpf:
    # derivative with respect to u=log(r)
    u0 = mp.log(r)
    return mp.diff(lambda u: func(mp.e ** u), u0, order)


def threshold_stationarity(base_masses: Dict[str, mp.mpf], alpha_inv0: mp.mpf, delta_gap: mp.mpf, Q2: mp.mpf, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf]) -> Dict[str, object]:
    rstar = Q2 / base_masses["mu"]
    me_over_mu = base_masses["e"] / base_masses["mu"]

    def mu_exact_of_r(r: mp.mpf) -> mp.mpf:
        return pi_kernel_component("exact_spacelike", r * base_masses["mu"], base_masses["mu"], alpha_inv0, xs, ws)

    def q2_channel_of_r(r: mp.mpf) -> mp.mpf:
        Q = r * base_masses["mu"]
        return pi_kernel_component("exact_spacelike", Q, base_masses["e"], alpha_inv0, xs, ws) + pi_kernel_component("exact_spacelike", Q, base_masses["mu"], alpha_inv0, xs, ws)

    value = mu_exact_of_r(rstar)
    d1 = derivative_log_r(mu_exact_of_r, rstar, 1)
    d2 = derivative_log_r(mu_exact_of_r, rstar, 2)
    ch = q2_channel_of_r(rstar)
    ch_d1 = derivative_log_r(q2_channel_of_r, rstar, 1)
    ch_d2 = derivative_log_r(q2_channel_of_r, rstar, 2)
    return {
        "rstar": ffloat(rstar),
        "two_over_rstar": ffloat(2/rstar),
        "mu_exact_value": ffloat(value),
        "mu_exact_Dlog1": ffloat(d1),
        "mu_exact_Dlog2": ffloat(d2),
        "abs_Dlog1_over_value": ffloat(abs(d1)/abs(value)),
        "abs_Dlog2_over_value": ffloat(abs(d2)/abs(value)),
        "full_q2_value": ffloat(ch),
        "full_q2_Dlog1": ffloat(ch_d1),
        "full_q2_Dlog2": ffloat(ch_d2),
        "full_q2_abs_Dlog1_over_target": ffloat(abs(ch_d1)/abs(delta_gap)),
        "one_percent_r_shift_rel_effect": ffloat(mp.mpf("0.01")*abs(ch_d1)/abs(delta_gap)),
        "near_stationary": bool(abs(d1) <= mp.mpf("0.05")*abs(value)),
        "near_inflection": bool(abs(d2) <= mp.mpf("0.05")*abs(value)),
    }


def mass_ratio_origin(base_masses: Dict[str, mp.mpf], alpha_inv0: mp.mpf, delta_gap: mp.mpf, target_score: LawScore, xs: Sequence[mp.mpf], ws: Sequence[mp.mpf], rho_span: float, rho_steps: int, ratio2d_span: float, ratio2d_steps: int, synthetic_trials: int, synthetic_band: float, seed: int) -> Dict[str, object]:
    rho_obs = base_masses["pi"] / base_masses["mu"]
    rows_rho: List[Dict[str, object]] = []
    half = rho_steps // 2
    for i in range(rho_steps):
        d = -rho_span + 2*rho_span*i/(rho_steps-1)
        pm = dict(base_masses)
        pm["pi"] = base_masses["mu"] * rho_obs * mp.e**d
        Q1, Q2 = frozen_Qs(pm)
        sc = score_two_channel(delta_gap, Q1, Q2, pm, alpha_inv0, xs, ws)
        rows_rho.append({"i": i, "dlog_rho": d, "rho": ffloat(pm["pi"]/pm["mu"]), "rho_over_obs": ffloat(mp.e**d), "is_observed": i == half, **asdict(sc)})
    rho_best = min(rows_rho, key=lambda r: r["joint_max"])
    rho_p = empirical_p(sum(1 for r in rows_rho if r["joint_max"] <= target_score.joint_max), len(rows_rho))

    grid = local_mass_ratio_grid(base_masses, alpha_inv0, delta_gap, target_score, ratio2d_span, ratio2d_steps, xs, ws)

    # reduced-mass alternatives
    mu, pi, e = base_masses["mu"], base_masses["pi"], base_masses["e"]
    alt_Q2s = [
        ("sqrt2pi_reduced_mu_pi_TARGET", mp.sqrt(2*mp.pi)*reduced(mu, pi), True),
        ("sqrt2pi_mu", mp.sqrt(2*mp.pi)*mu, False),
        ("sqrt2pi_pi", mp.sqrt(2*mp.pi)*pi, False),
        ("sqrt2pi_arithmetic_mean_mu_pi", mp.sqrt(2*mp.pi)*(mu+pi)/2, False),
        ("sqrt2pi_geometric_mean_mu_pi", mp.sqrt(2*mp.pi)*mp.sqrt(mu*pi), False),
        ("sqrt2pi_harmonic_mean_mu_pi", mp.sqrt(2*mp.pi)*(2*mu*pi/(mu+pi)), False),
        ("sqrt2pi_mass_difference", mp.sqrt(2*mp.pi)*abs(pi-mu), False),
        ("sqrt2pi_red_e_pi_DIAG", mp.sqrt(2*mp.pi)*reduced(e, pi), False),
        ("sqrt2pi_red_e_mu_DIAG", mp.sqrt(2*mp.pi)*reduced(e, mu), False),
    ]
    Q1, Q2 = frozen_Qs(base_masses)
    rm_rows = []
    for label, q2, target in alt_Q2s:
        sc = score_two_channel(delta_gap, Q1, q2, base_masses, alpha_inv0, xs, ws)
        rm_rows.append({"label": label, "target": target, **asdict(sc)})
    rm_rows = sorted(rm_rows, key=lambda r: r["joint_max"])
    rm_rank = 1 + [r["label"] for r in rm_rows].index("sqrt2pi_reduced_mu_pi_TARGET")

    # synthetic controls: draw rho near observed, keep mu fixed, update pi; also vary e/mu sometimes.
    rng = random.Random(seed)
    synth_rows = []
    for i in range(synthetic_trials):
        pm = dict(base_masses)
        pm["pi"] = base_masses["pi"] * mp.e ** rng.uniform(-synthetic_band, synthetic_band)
        pm["e"] = base_masses["e"] * mp.e ** rng.uniform(-synthetic_band, synthetic_band)
        Q1s, Q2s = frozen_Qs(pm)
        sc = score_two_channel(delta_gap, Q1s, Q2s, pm, alpha_inv0, xs, ws)
        synth_rows.append({"i": i, "rho_over_obs": ffloat((pm["pi"]/pm["mu"])/rho_obs), **asdict(sc)})
    synth_j = sum(1 for r in synth_rows if r["joint_max"] <= target_score.joint_max)
    synth_b = sum(1 for r in synth_rows if r["rel1"] <= target_score.rel1 and r["rel2"] <= target_score.rel2)

    # sensitivity around rho observed via central finite difference in log-rho for joint_rms.
    h = mp.mpf("1e-4")
    def joint_rms_at_d(d: mp.mpf) -> mp.mpf:
        pm = dict(base_masses); pm["pi"] = base_masses["pi"] * mp.e**d
        Q1t, Q2t = frozen_Qs(pm)
        sc = score_two_channel(delta_gap, Q1t, Q2t, pm, alpha_inv0, xs, ws)
        return mp.mpf(str(sc.joint_rms))
    deriv = (joint_rms_at_d(h) - joint_rms_at_d(-h))/(2*h)

    return {
        "rho_obs": ffloat(rho_obs),
        "rho_rows": rows_rho,
        "rho_best": rho_best,
        "rho_p_joint": rho_p,
        "ratio2d": grid,
        "reduced_mass_rows": rm_rows,
        "reduced_mass_rank": rm_rank,
        "synthetic_rows": synth_rows,
        "synthetic_joint_p": empirical_p(synth_j, synthetic_trials),
        "synthetic_both_p": empirical_p(synth_b, synthetic_trials),
        "dlog_joint_rms": ffloat(deriv),
        "abs_dlog_joint_rms_over_value": abs(ffloat(deriv)) / target_score.joint_rms,
    }


# =============================================================================
# CLI / main
# =============================================================================


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="ROT RH / Xi-QED shared-law master evidence audit (GitHub-ready).")
    p.add_argument("--dps", type=int, default=100)
    p.add_argument("--alpha0-inv", type=str, default="137.035999177")
    p.add_argument("--me-MeV", type=str, default="0.51099895000")
    p.add_argument("--mmu-MeV", type=str, default="105.6583755")
    p.add_argument("--mtau-MeV", type=str, default="1776.86")
    p.add_argument("--mpi-charged-MeV", type=str, default="139.57039")
    p.add_argument("--quad-n", type=int, default=256)
    p.add_argument("--quad-n-list", type=str, default="128,256,512")
    p.add_argument("--random-delta-controls", type=int, default=2000)
    p.add_argument("--random-delta-band", type=str, default="0.5,1.8")
    p.add_argument("--pseudo-mass-trials", type=int, default=512)
    p.add_argument("--pseudo-mass-families", type=str, default="ratio_jitter_e_pi_only,ratio_jitter_pi_only,ratio_jitter_all,preserve_mu_pi_break_epi,random_ratio_ordered_wide,affine_plus_small_ratio_noise,preserve_e_pi_break_mupi,affine_scale_all")
    p.add_argument("--pseudo-mass-jitter", type=float, default=0.25)
    p.add_argument("--local-mass-grid-span", type=float, default=0.30)
    p.add_argument("--local-mass-grid-steps", type=int, default=41)
    p.add_argument("--rho-scan-log-span", type=float, default=0.50)
    p.add_argument("--rho-scan-steps", type=int, default=401)
    p.add_argument("--ratio2d-log-span", type=float, default=0.35)
    p.add_argument("--ratio2d-steps", type=int, default=41)
    p.add_argument("--synthetic-ratio-trials", type=int, default=2000)
    p.add_argument("--synthetic-ratio-log-band", type=float, default=0.50)
    p.add_argument("--seed", type=int, default=314159)
    p.add_argument("--interesting-rel-err", type=float, default=1e-3)
    p.add_argument("--joint-rel-err", type=float, default=1e-3)
    p.add_argument("--max-joint-p", type=float, default=0.05)
    p.add_argument("--strong-joint-p", type=float, default=0.01)
    p.add_argument("--scheme-degradation-min", type=float, default=25.0)
    p.add_argument("--threshold-residual-ratio-min", type=float, default=0.05)
    p.add_argument("--sensitivity-min", type=float, default=0.02)
    p.add_argument("--protocol-only", action="store_true")
    p.add_argument("--out-prefix", type=str, default="xi_qed_shared_law_master_evidence")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    t0 = time.time()
    mp.mp.dps = args.dps

    protocol = {
        "audit": "rot_rh_xi_qed_shared_law_master_evidence_audit",
        "version": "1.0.0",
        "search_disabled": True,
        "frozen_law": {
            "A0": "2*pi/K2",
            "gap": "alpha0_inv - A0",
            "Q1": "(4/pi)*m_pi_charged",
            "Q1_loops": ["e"],
            "Q2": "sqrt(2*pi)*m_mu*m_pi/(m_mu+m_pi)",
            "Q2_loops": ["e", "mu"],
            "scheme": "exact_spacelike_one_loop_on_shell",
        },
        "args": vars(args),
    }
    protocol_hash = sha256_json(protocol)
    protocol["protocol_hash"] = protocol_hash

    write_json(f"{args.out_prefix}_master_protocol.json", protocol)

    banner("ROT RH / XI-QED SHARED LAW — MASTER EVIDENCE AUDIT")
    print(f"time                    : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"dps                     : {args.dps}")
    print(f"alpha(0)^-1             : {args.alpha0_inv}")
    print(f"masses MeV              : e={args.me_MeV}, mu={args.mmu_MeV}, tau={args.mtau_MeV}, pi_ch={args.mpi_charged_MeV}")
    print("frozen law              : Q1=(4/pi)m_pi, loop=e ; Q2=sqrt(2pi)*red(mu,pi), loop=e+mu")
    print("search                  : DISABLED — master evidence only; no coefficient/formula/additive-scale discovery")
    print(f"protocol hash           : {protocol_hash}")
    print(f"out_prefix              : {args.out_prefix}")
    print(SUBLINE)

    if args.protocol_only:
        print("Protocol-only mode: wrote protocol JSON and exiting.")
        return 0

    alpha_inv0 = mp.mpf(args.alpha0_inv)
    masses = {
        "e": mp.mpf(args.me_MeV),
        "mu": mp.mpf(args.mmu_MeV),
        "tau": mp.mpf(args.mtau_MeV),
        "pi": mp.mpf(args.mpi_charged_MeV),
    }

    xs, ws = gauss_legendre_nodes(args.quad_n, args.dps)

    section("COMPUTING XI CURVATURE COUNT")
    xi = compute_xi_invariants(args.dps)
    K2, A0, R4 = xi["K2"], xi["A0"], xi["R4"]
    gap = alpha_inv0 - A0
    print(f"K2 = -d2 log Xi(0)       : {mp.nstr(K2, 24)}")
    print(f"A0 = 2*pi/K2             : {mp.nstr(A0, 18)}")
    print(f"R4 = L4/K2^2             : {mp.nstr(R4, 18)}")
    print(f"delta alpha^-1 - A0      : {mp.nstr(gap, 18)}")
    print(f"relative gap             : {mp.nstr(gap/alpha_inv0, 15)}")

    section("FROZEN LAW BASELINE SCORE")
    Q1, Q2 = frozen_Qs(masses)
    score = score_two_channel(gap, Q1, Q2, masses, alpha_inv0, xs, ws)
    print(f"Q1=(4/pi)m_pi            : {score.Q1:.12f} MeV")
    print(f"Q2=sqrt2pi*red(mu,pi)    : {score.Q2:.12f} MeV")
    print(f"Q2/m_mu                  : {ffloat(Q2/masses['mu']):.15f}")
    print(f"Delta1 e-loop            : {score.delta1:.13g} rel={score.rel1:.6e}")
    print(f"Delta2 e+mu-loop         : {score.delta2:.13g} rel={score.rel2:.6e}")
    print(f"joint max/rms/mean rel   : {score.joint_max:.6e} / {score.joint_rms:.6e} / {score.joint_mean:.6e}")
    print(f"prediction agreement rel : {score.agreement_rel:.6e}")
    forbidden_Q = (masses["pi"] + masses["e"]) * (4/mp.pi)
    forbidden_delta = delta_for_loops(forbidden_Q, masses, ("e",), "exact_spacelike", alpha_inv0, xs, ws)
    forbidden_rel = relerr(forbidden_delta, gap)
    print(f"forbidden mixed diagnostic: Q=(m_pi+m_e)*4/pi = {ffloat(forbidden_Q):.12f} MeV rel={ffloat(forbidden_rel):.12e}")

    section("QUADRATURE STABILITY")
    stability_rows = []
    quad_list = [int(x.strip()) for x in args.quad_n_list.split(",") if x.strip()]
    for qn in quad_list:
        xq, wq = gauss_legendre_nodes(qn, args.dps)
        scq = score_two_channel(gap, Q1, Q2, masses, alpha_inv0, xq, wq)
        row = {"quad_n": qn, **asdict(scq)}
        stability_rows.append(row)
        print(f"quad={qn:4d} rel1={scq.rel1:.12e} rel2={scq.rel2:.12e} joint={scq.joint_max:.12e} agree={scq.agreement_rel:.12e}")
    stability_span = max(r["joint_max"] for r in stability_rows) - min(r["joint_max"] for r in stability_rows)

    section("RANDOM-DELTA CONTROLS")
    band = tuple(float(x.strip()) for x in args.random_delta_band.split(","))
    rd = random_delta_controls(gap, mp.mpf(str(score.delta1)), mp.mpf(str(score.delta2)), args.random_delta_controls, band, args.seed + 11)
    print(f"random both-individual p : {rd['both_p']:.8g} ({rd['both_count']}/{args.random_delta_controls})")
    print(f"random joint-max p       : {rd['jointmax_p']:.8g} ({rd['jointmax_count']}/{args.random_delta_controls})")
    print(f"random joint-rms p       : {rd['jointrms_p']:.8g} ({rd['jointrms_count']}/{args.random_delta_controls})")

    section("PSEUDO-MASS CONTROLS")
    families = [x.strip() for x in args.pseudo_mass_families.split(",") if x.strip()]
    pmc = pseudo_mass_controls(masses, alpha_inv0, gap, score, families, args.pseudo_mass_trials, args.pseudo_mass_jitter, xs, ws, args.seed + 22)
    for fs in pmc["family_summary"]:
        print(f"{fs['family']:<35s} diag={str(fs['diagnostic']):<5s} both_p={fs['both_p']:.6g} jointmax_p={fs['jointmax_p']:.6g} med_joint={fs['median_joint']:.3e}")
    print(SUBLINE)
    print(f"main pseudo rows          : {pmc['main']['rows']}")
    print(f"main pseudo both p        : {pmc['main']['both_p']:.8g} ({pmc['main']['both_count']}/{pmc['main']['rows']})")
    print(f"main pseudo joint-max p   : {pmc['main']['jointmax_p']:.8g} ({pmc['main']['jointmax_count']}/{pmc['main']['rows']})")

    section("LOCAL MASS-RATIO GRID CONTROLS")
    grid = local_mass_ratio_grid(masses, alpha_inv0, gap, score, args.local_mass_grid_span, args.local_mass_grid_steps, xs, ws)
    print(f"grid steps                : {args.local_mass_grid_steps} x {args.local_mass_grid_steps} = {args.local_mass_grid_steps**2}")
    print(f"grid joint-max p          : {grid['jointmax_p']:.8g}")
    print(f"grid both p               : {grid['both_p']:.8g}")
    print(f"grid offdiag joint p      : {grid['offdiag_joint_p']:.8g}")
    print(f"best grid dlog(e/mu),dlog(pi/mu): {grid['best']['dlog_e_mu']:.6g}, {grid['best']['dlog_pi_mu']:.6g}")

    section("SCHEME ROBUSTNESS")
    sr = scheme_robustness(masses, alpha_inv0, gap, Q1, Q2, score, xs, ws)
    print(f"exact-spacelike baseline loop rank  : {sr['baseline_exact_rank']} / 4")
    print(f"asymp -5/3 baseline joint ratio     : {sr['asymp_ratio']:.6g}")
    print(f"asymp log-only baseline joint ratio : {sr['logonly_ratio']:.6g}")
    print("Top scheme rows:")
    for i, r in enumerate(sr["rows"][:8], 1):
        print(f"{i:02d}. joint={r['joint_max']:.6e} rel1={r['rel1']:.3e} rel2={r['rel2']:.3e} scheme={r['scheme']} loops=({r['q1_loops']},{r['q2_loops']}) {r['loop_label']}")

    section("THRESHOLD STRUCTURE MECHANISM")
    th = threshold_structure(masses, alpha_inv0, gap, Q1, Q2, score, xs, ws)
    for c in th["components"]:
        print(f"{c['channel']} {c['lepton']:>2s}: Q/m={c['Q_over_m']:.6g} exact={c['exact']:.12g} asymp={c['asymp_minus_5over3']:.12g} resid={c['finite_mass_residual']:.12g} resid/exact={c['residual_over_exact']:.6g}")
    print(f"Q2 mu residual/exact                : {th['q2_mu_residual_over_exact']:.12g}")
    print(f"Q2 exact->mu-asymp degradation      : {th['q2_exact_to_mu_asymp_degradation']:.6g}x")
    print(f"Q2 exact->mu-log degradation        : {th['q2_exact_to_mu_log_degradation']:.6g}x")
    print(f"best hybrid                         : {th['best_hybrid']}")

    section("THRESHOLD STATIONARITY / SENSITIVITY")
    st = threshold_stationarity(masses, alpha_inv0, gap, Q2, xs, ws)
    print(f"r*=Q2/m_mu                          : {st['rstar']:.15f}")
    print(f"exact |Dlog1|/|value|               : {st['abs_Dlog1_over_value']:.6g}")
    print(f"exact |Dlog2|/|value|               : {st['abs_Dlog2_over_value']:.6g}")
    print(f"full Q2 |Dlog1|/target              : {st['full_q2_abs_Dlog1_over_target']:.12g}")
    print(f"1% r-shift approx rel effect        : {st['one_percent_r_shift_rel_effect']:.12g}")
    print(f"near stationary / inflection        : {st['near_stationary']} / {st['near_inflection']}")

    section("MASS-RATIO ORIGIN")
    mro = mass_ratio_origin(masses, alpha_inv0, gap, score, xs, ws, args.rho_scan_log_span, args.rho_scan_steps, args.ratio2d_log_span, args.ratio2d_steps, args.synthetic_ratio_trials, args.synthetic_ratio_log_band, args.seed + 33)
    print(f"rho*=m_pi/m_mu                      : {mro['rho_obs']:.15f}")
    print(f"rho scan best rho/obs               : {mro['rho_best']['rho_over_obs']:.12g}")
    print(f"rho scan p joint                    : {mro['rho_p_joint']:.8g}")
    print(f"ratio2d p joint                     : {mro['ratio2d']['jointmax_p']:.8g}")
    print(f"synthetic joint/both p              : {mro['synthetic_joint_p']:.8g} / {mro['synthetic_both_p']:.8g}")
    print(f"target reduced-mass rank            : {mro['reduced_mass_rank']} / {len(mro['reduced_mass_rows'])}")
    print(f"|Dlog joint_rms|/value              : {mro['abs_dlog_joint_rms_over_value']:.6g}")

    # Evidence flags
    passes = {
        "baseline_rel": score.rel1 <= args.interesting_rel_err and score.rel2 <= args.interesting_rel_err and score.joint_max <= args.joint_rel_err,
        "stability": stability_span < 1e-7,
        "random": rd["jointmax_p"] <= args.max_joint_p and rd["both_p"] <= args.max_joint_p,
        "pseudo": pmc["main"]["jointmax_p"] <= args.max_joint_p and pmc["main"]["both_p"] <= args.max_joint_p,
        "local_grid": grid["jointmax_p"] <= args.max_joint_p and grid["offdiag_joint_p"] <= args.max_joint_p,
        "scheme_exact_best": sr["baseline_exact_best"] and sr["asymp_ratio"] >= args.scheme_degradation_min,
        "threshold_exact_kernel": th["best_hybrid"] == "baseline_all_exact" and th["q2_mu_residual_over_exact"] >= args.threshold_residual_ratio_min and th["q2_exact_to_mu_asymp_degradation"] >= args.scheme_degradation_min,
        "threshold_sensitive": st["full_q2_abs_Dlog1_over_target"] >= args.sensitivity_min,
        "mass_ratio_origin": mro["rho_p_joint"] <= args.max_joint_p and mro["ratio2d"]["jointmax_p"] <= args.max_joint_p and mro["synthetic_joint_p"] <= args.max_joint_p and mro["reduced_mass_rank"] == 1,
    }
    all_core = all(passes.values())
    strong = all_core and rd["jointmax_p"] <= args.strong_joint_p and pmc["main"]["jointmax_p"] <= args.strong_joint_p and mro["synthetic_joint_p"] <= args.strong_joint_p
    if strong:
        global_flag = "XI_QED_SHARED_LAW_MASTER_EVIDENCE_HARDENED_STRONG"
    elif all_core:
        global_flag = "XI_QED_SHARED_LAW_MASTER_EVIDENCE_CONSOLIDATED_INTERESTING"
    elif passes["baseline_rel"]:
        global_flag = "XI_QED_SHARED_LAW_MASTER_EVIDENCE_BASELINE_ONLY"
    else:
        global_flag = "XI_QED_SHARED_LAW_MASTER_EVIDENCE_NOT_CONFIRMED"

    evidence_rows = [
        {"block": "baseline", "metric": "joint_max", "value": score.joint_max, "pass": passes["baseline_rel"]},
        {"block": "stability", "metric": "joint_span", "value": stability_span, "pass": passes["stability"]},
        {"block": "random_delta", "metric": "jointmax_p", "value": rd["jointmax_p"], "pass": passes["random"]},
        {"block": "random_delta", "metric": "both_p", "value": rd["both_p"], "pass": passes["random"]},
        {"block": "pseudo_mass", "metric": "main_jointmax_p", "value": pmc["main"]["jointmax_p"], "pass": passes["pseudo"]},
        {"block": "pseudo_mass", "metric": "main_both_p", "value": pmc["main"]["both_p"], "pass": passes["pseudo"]},
        {"block": "local_grid", "metric": "jointmax_p", "value": grid["jointmax_p"], "pass": passes["local_grid"]},
        {"block": "scheme", "metric": "asymp_degradation_ratio", "value": sr["asymp_ratio"], "pass": passes["scheme_exact_best"]},
        {"block": "threshold", "metric": "q2_mu_residual_over_exact", "value": th["q2_mu_residual_over_exact"], "pass": passes["threshold_exact_kernel"]},
        {"block": "stationarity", "metric": "full_q2_Dlog1_over_target", "value": st["full_q2_abs_Dlog1_over_target"], "pass": passes["threshold_sensitive"]},
        {"block": "mass_ratio", "metric": "rho_scan_p", "value": mro["rho_p_joint"], "pass": passes["mass_ratio_origin"]},
        {"block": "mass_ratio", "metric": "reduced_mass_rank", "value": mro["reduced_mass_rank"], "pass": passes["mass_ratio_origin"]},
    ]

    section("XI-QED SHARED LAW MASTER EVIDENCE SUMMARY")
    print(f"global flag                         : {global_flag}")
    print(f"protocol hash                       : {protocol_hash}")
    print(f"A0                                  : {mp.nstr(A0, 18)}")
    print(f"delta needed                        : {mp.nstr(gap, 18)}")
    print(f"Q1,Q2                               : {score.Q1:.12f} MeV, {score.Q2:.12f} MeV")
    print(f"baseline rel1, rel2                 : {score.rel1:.12e}, {score.rel2:.12e}")
    print(f"baseline joint max/rms              : {score.joint_max:.12e}, {score.joint_rms:.12e}")
    print(f"random both / jointmax p            : {rd['both_p']:.8g} / {rd['jointmax_p']:.8g}")
    print(f"main pseudo both / jointmax p       : {pmc['main']['both_p']:.8g} / {pmc['main']['jointmax_p']:.8g}")
    print(f"local grid jointmax / offdiag p     : {grid['jointmax_p']:.8g} / {grid['offdiag_joint_p']:.8g}")
    print(f"scheme asymp degradation            : {sr['asymp_ratio']:.6g}x")
    print(f"Q2 mu residual/exact                : {th['q2_mu_residual_over_exact']:.12g}")
    print(f"Q2 exact->mu-asymp degradation      : {th['q2_exact_to_mu_asymp_degradation']:.6g}x")
    print(f"threshold sensitivity Dlog/target   : {st['full_q2_abs_Dlog1_over_target']:.12g}")
    print(f"rho scan / ratio2d / synthetic p    : {mro['rho_p_joint']:.8g} / {mro['ratio2d']['jointmax_p']:.8g} / {mro['synthetic_joint_p']:.8g}")
    print(f"target reduced-mass rank            : {mro['reduced_mass_rank']} / {len(mro['reduced_mass_rows'])}")
    print(f"forbidden mixed diagnostic rel      : {ffloat(forbidden_rel):.12e}")
    print(f"stability span                      : {stability_span:.3e}")
    print("passes:")
    for k, v in passes.items():
        print(f"  {k:<25s}: {v}")
    print(SUBLINE)

    # Write outputs
    prefix = args.out_prefix
    write_csv(f"{prefix}_master_evidence_table.csv", evidence_rows)
    write_csv(f"{prefix}_master_stability.csv", stability_rows)
    write_csv(f"{prefix}_master_random_delta_controls.csv", rd["rows"])
    write_csv(f"{prefix}_master_pseudo_controls.csv", pmc["rows"])
    write_csv(f"{prefix}_master_pseudo_family_summary.csv", pmc["family_summary"])
    write_csv(f"{prefix}_master_local_mass_grid.csv", grid["rows"])
    write_csv(f"{prefix}_master_scheme_rows.csv", sr["rows"])
    write_csv(f"{prefix}_master_threshold_components.csv", th["components"])
    write_csv(f"{prefix}_master_threshold_hybrids.csv", th["hybrids"])
    write_csv(f"{prefix}_master_rho_scan.csv", mro["rho_rows"])
    write_csv(f"{prefix}_master_ratio2d_scan.csv", mro["ratio2d"]["rows"])
    write_csv(f"{prefix}_master_reduced_mass_diagnostics.csv", mro["reduced_mass_rows"])
    write_csv(f"{prefix}_master_synthetic_ratio_controls.csv", mro["synthetic_rows"])

    summary = {
        "global_flag": global_flag,
        "protocol_hash": protocol_hash,
        "K2": K2,
        "A0": A0,
        "R4": R4,
        "delta_needed": gap,
        "relative_gap": gap/alpha_inv0,
        "Q1": score.Q1,
        "Q2": score.Q2,
        "baseline_score": asdict(score),
        "forbidden_mixed_diagnostic": {"Q": ffloat(forbidden_Q), "delta": ffloat(forbidden_delta), "rel": ffloat(forbidden_rel)},
        "random_delta": {k: v for k, v in rd.items() if k != "rows"},
        "pseudo_mass_main": pmc["main"],
        "pseudo_mass_diagnostic": pmc["diagnostic"],
        "local_grid": {"jointmax_p": grid["jointmax_p"], "both_p": grid["both_p"], "offdiag_joint_p": grid["offdiag_joint_p"], "best": grid["best"]},
        "scheme": {k: v for k, v in sr.items() if k != "rows"},
        "threshold": {k: v for k, v in th.items() if k not in ("components", "hybrids")},
        "stationarity": st,
        "mass_ratio_origin": {k: v for k, v in mro.items() if k not in ("rho_rows", "ratio2d", "reduced_mass_rows", "synthetic_rows")},
        "passes": passes,
        "stability_span": stability_span,
    }
    write_json(f"{prefix}_master_summary.json", summary)
    meta = {"runtime_sec": time.time() - t0, "argv": sys.argv, "python": sys.version, "platform": sys.platform, "files_written_prefix": prefix}
    write_json(f"{prefix}_master_meta.json", meta)

    print("Files written")
    for suffix in [
        "master_protocol.json", "master_summary.json", "master_meta.json", "master_evidence_table.csv",
        "master_stability.csv", "master_random_delta_controls.csv", "master_pseudo_controls.csv",
        "master_pseudo_family_summary.csv", "master_local_mass_grid.csv", "master_scheme_rows.csv",
        "master_threshold_components.csv", "master_threshold_hybrids.csv", "master_rho_scan.csv",
        "master_ratio2d_scan.csv", "master_reduced_mass_diagnostics.csv", "master_synthetic_ratio_controls.csv",
    ]:
        print(f"  {prefix}_{suffix}")
    print(LINE)
    print(f"Total runtime                   : {time.time() - t0:.2f}s")
    print(LINE)
    print(f"AUDIT FLAG: {global_flag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
