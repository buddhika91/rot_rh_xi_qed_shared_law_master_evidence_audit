#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT RH / XI-QED SHARED LAW — PREREGISTERED LOCK AUDIT
======================================================

Purpose
-------
This is the no-search, preregistered lock test for the current best alpha-gap branch.
It freezes the law before scoring:

    A0     = 2*pi/K2, with K2 = - d^2/dgamma^2 log Xi(1/2+i gamma)|_{gamma=0}
    Delta  = alpha(0)^-1 - A0

    Q1     = (4/pi) * m_pi
    loop1  = electron only

    Q2     = sqrt(2*pi) * red(m_mu, m_pi)
           = sqrt(2*pi) * m_mu*m_pi/(m_mu+m_pi)
    loop2  = electron + muon

No coefficient menu.
No formula search.
No additive mixed-mass external Q scales.
No grammar-best candidate can determine pass/fail.

The audit tests:
  1. Frozen-law closeness to Delta.
  2. Quadrature stability.
  3. Random target-delta controls.
  4. Pseudo-mass controls, separated into main ratio-changing families and diagnostic degeneracies.
  5. Local mass-ratio grid around the real point.
  6. Optional protocol-only mode writing a JSON preregistration file and hash.

Interpretation
--------------
This script can support only a limited claim:

    The frozen two-scale law is an interesting numerical/physical signal if it
    beats random/pseudo controls.

It cannot prove a physical derivation of alpha, QED, or RH.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import sys
import time
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple, Any

import numpy as np
import mpmath as mp


# ----------------------------- formatting -----------------------------

def line(ch: str = "=") -> None:
    print(ch * 120)


def fmt(x: float, nd: int = 15) -> str:
    if not math.isfinite(float(x)):
        return str(x)
    return f"{float(x):.{nd}g}"


def safe_float(x: Any) -> float:
    """Convert mp/complex/string values that may be '(a + 0.0j)' to real float."""
    if isinstance(x, (float, int, np.floating, np.integer)):
        return float(x)
    if isinstance(x, complex):
        return float(x.real)
    if isinstance(x, mp.mpc):
        return float(mp.re(x))
    if isinstance(x, mp.mpf):
        return float(x)
    s = str(x).strip()
    try:
        return float(s)
    except Exception:
        pass
    try:
        z = complex(s.replace(" ", ""))
        return float(z.real)
    except Exception:
        pass
    # mpmath can parse a lot, but not Python complex text with j.
    if "j" in s:
        s2 = s.replace("j", "")
        # crude fallback for '(a+0.0)'
        s2 = s2.strip("()")
        if "+" in s2:
            return float(s2.split("+")[0])
        if "-" in s2[1:]:
            parts = s2.rsplit("-", 1)
            return float(parts[0])
    return float(mp.mpf(s))


# ----------------------------- Xi invariant -----------------------------

def xi_completed(s: mp.mpc) -> mp.mpc:
    return mp.mpf("0.5") * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def compute_xi_invariants(dps: int) -> Dict[str, str]:
    mp.mp.dps = dps

    def L(g):
        z = xi_completed(mp.mpf("0.5") + 1j * g)
        return mp.log(z)

    # L(gamma) = L0 + L2 gamma^2/2! + L4 gamma^4/4! + ...
    # K2 = -L''(0), because log Xi decays with gamma near zero.
    L2 = mp.diff(L, mp.mpf("0"), 2)
    L4 = mp.diff(L, mp.mpf("0"), 4)
    K2 = -mp.re(L2)
    R4 = mp.re(L4) / (K2 ** 2)
    A0 = 2 * mp.pi / K2
    return {
        "K2": mp.nstr(K2, 80),
        "A0": mp.nstr(A0, 80),
        "R4": mp.nstr(R4, 80),
    }


# ----------------------------- QED vacuum polarization -----------------------------

def gauss_legendre_nodes(n: int) -> Tuple[np.ndarray, np.ndarray]:
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (x + 1.0)
    wt = 0.5 * w
    return t.astype(np.float64), wt.astype(np.float64)


@dataclass
class Masses:
    me: float
    mmu: float
    mtau: float
    mpi: float


def delta_alpha_inv_spacelike_one_loop(Q: float, loop_masses: List[float], quad_n: int, nodes_cache: Dict[int, Tuple[np.ndarray, np.ndarray]]) -> float:
    """Exact one-loop spacelike vacuum-polarization shift in inverse alpha.

    Delta alpha^{-1}(Q) = sum_l (2/(pi)) int_0^1 dx x(1-x) log(1 + Q^2 x(1-x)/m_l^2)

    This convention matches the previous scripts in this sequence: alpha^{-1}(Q)=alpha^{-1}(0)-Delta.
    """
    if quad_n not in nodes_cache:
        nodes_cache[quad_n] = gauss_legendre_nodes(quad_n)
    x, w = nodes_cache[quad_n]
    xx = x * (1.0 - x)
    total = 0.0
    Q2 = float(Q) * float(Q)
    for m in loop_masses:
        z = Q2 * xx / (float(m) * float(m))
        total += float(np.sum(w * xx * np.log1p(z)))
    return (2.0 / math.pi) * total


def red_mass(a: float, b: float) -> float:
    return (a * b) / (a + b)


def score_frozen_law(delta_target: float, masses: Masses, quad_n: int, nodes_cache: Dict[int, Tuple[np.ndarray, np.ndarray]]) -> Dict[str, float]:
    Q1 = (4.0 / math.pi) * masses.mpi
    Q2 = math.sqrt(2.0 * math.pi) * red_mass(masses.mmu, masses.mpi)
    d1 = delta_alpha_inv_spacelike_one_loop(Q1, [masses.me], quad_n, nodes_cache)
    d2 = delta_alpha_inv_spacelike_one_loop(Q2, [masses.me, masses.mmu], quad_n, nodes_cache)
    rel1 = abs(d1 - delta_target) / abs(delta_target)
    rel2 = abs(d2 - delta_target) / abs(delta_target)
    joint_max = max(rel1, rel2)
    joint_rms = math.sqrt(0.5 * (rel1 * rel1 + rel2 * rel2))
    joint_mean = 0.5 * (rel1 + rel2)
    pred_agreement = abs(d2 - d1) / abs(delta_target)
    pred_agreement_internal = abs(d2 - d1) / max(abs(d1), abs(d2), 1e-300)
    both_individual = (rel1, rel2)
    return {
        "Q1_pion_4overpi": Q1,
        "Q2_reduced_sqrt2pi": Q2,
        "delta1_e_loop": d1,
        "delta2_emu_loop": d2,
        "rel1_pion": rel1,
        "rel2_reduced": rel2,
        "joint_max_rel": joint_max,
        "joint_rms_rel": joint_rms,
        "joint_mean_rel": joint_mean,
        "prediction_agreement_rel_to_delta": pred_agreement,
        "prediction_agreement_internal": pred_agreement_internal,
    }


# ----------------------------- pseudo mass generation -----------------------------

def jitter_log(x: float, sigma: float, rng: random.Random) -> float:
    return x * math.exp(rng.uniform(-sigma, sigma))


def generate_pseudo_masses(family: str, base: Masses, jitter: float, rng: random.Random) -> Tuple[Masses, bool]:
    """Return masses and whether family is diagnostic/degenerated rather than main."""
    me, mmu, mtau, mpi = base.me, base.mmu, base.mtau, base.mpi
    diagnostic = False

    if family == "ratio_jitter_e_pi_only":
        # Keep pion/muon/tau fixed; perturb e/pi ratio.
        me2 = jitter_log(me, jitter, rng)
        return Masses(me2, mmu, mtau, mpi), diagnostic

    if family == "ratio_jitter_pi_only":
        # Perturb pion absolute scale, keeping electron/muon/tau fixed.
        mpi2 = jitter_log(mpi, jitter, rng)
        return Masses(me, mmu, mtau, mpi2), diagnostic

    if family == "ratio_jitter_all":
        return Masses(jitter_log(me, jitter, rng), jitter_log(mmu, jitter, rng), jitter_log(mtau, jitter, rng), jitter_log(mpi, jitter, rng)), diagnostic

    if family == "preserve_mu_pi_break_epi":
        # Scale mu and pi together; break electron/pion ratio.
        s = math.exp(rng.uniform(-jitter, jitter))
        me2 = jitter_log(me, jitter, rng)
        return Masses(me2, mmu * s, mtau, mpi * s), diagnostic

    if family == "random_ratio_ordered_wide":
        # Draw ordered masses in broad windows preserving e << mu/pi ordering roughly.
        mpi2 = mpi * math.exp(rng.uniform(-2.0 * jitter, 2.0 * jitter))
        mmu2 = mpi2 * math.exp(rng.uniform(math.log(0.35), math.log(1.25)))
        me2 = mpi2 * math.exp(rng.uniform(math.log(0.0005), math.log(0.02)))
        mtau2 = max(mtau * math.exp(rng.uniform(-jitter, jitter)), 3.0 * max(mmu2, mpi2))
        return Masses(me2, mmu2, mtau2, mpi2), diagnostic

    if family == "affine_plus_small_ratio_noise":
        # Near scale degeneracy with small independent ratio perturbation. Main because ratios are not exactly fixed.
        s = math.exp(rng.uniform(-jitter, jitter))
        small = 0.15 * jitter
        return Masses(me * s * math.exp(rng.uniform(-small, small)),
                      mmu * s * math.exp(rng.uniform(-small, small)),
                      mtau * s * math.exp(rng.uniform(-small, small)),
                      mpi * s * math.exp(rng.uniform(-small, small))), diagnostic

    if family == "preserve_e_pi_break_mupi":
        # Diagnostic: pion/electron relation fixed so pion candidate is mechanically protected.
        diagnostic = True
        s = math.exp(rng.uniform(-jitter, jitter))
        mmu2 = jitter_log(mmu, jitter, rng)
        return Masses(me * s, mmu2, mtau, mpi * s), diagnostic

    if family == "affine_scale_all":
        # Diagnostic degeneracy: scales all masses together. Not a real ratio-changing control.
        diagnostic = True
        s = math.exp(rng.uniform(-jitter, jitter))
        return Masses(me * s, mmu * s, mtau * s, mpi * s), diagnostic

    raise ValueError(f"Unknown pseudo mass family: {family}")


def empirical_p(count_leq: int, n: int) -> float:
    return (count_leq + 1.0) / (n + 1.0)


# ----------------------------- CSV helpers -----------------------------

def write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        with open(path, "w", newline="", encoding="utf-8") as f:
            f.write("")
        return
    keys = []
    seen = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        for r in rows:
            wr.writerow(r)


def sha_json(obj: Any) -> str:
    data = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(data.encode("utf-8")).hexdigest()


# ----------------------------- main -----------------------------

def build_protocol(args: argparse.Namespace) -> Dict[str, Any]:
    return {
        "protocol_name": "ROT RH / Xi-QED shared law preregistered lock audit",
        "law_is_frozen": True,
        "no_coefficient_menu": True,
        "no_formula_search": True,
        "no_additive_mixed_mass_external_scales": True,
        "invariant": "A0 = 2*pi/K2, K2=-d2 log Xi(1/2+i gamma)|0",
        "gap": "Delta = alpha(0)^-1 - A0",
        "frozen_law": {
            "Q1": "(4/pi)*m_pi",
            "loop1": "electron only",
            "Q2": "sqrt(2*pi)*m_mu*m_pi/(m_mu+m_pi)",
            "loop2": "electron + muon",
        },
        "masses_MeV": {
            "me": args.me_MeV,
            "mmu": args.mmu_MeV,
            "mtau": args.mtau_MeV,
            "mpi_charged": args.mpi_charged_MeV,
        },
        "controls": {
            "random_delta_controls": args.random_delta_controls,
            "random_delta_band": args.random_delta_band,
            "pseudo_mass_trials": args.pseudo_mass_trials,
            "pseudo_mass_families": args.pseudo_mass_families,
            "pseudo_mass_jitter": args.pseudo_mass_jitter,
            "local_mass_grid_span": args.local_mass_grid_span,
            "local_mass_grid_steps": args.local_mass_grid_steps,
        },
        "thresholds": {
            "interesting_rel_err": args.interesting_rel_err,
            "strong_rel_err": args.strong_rel_err,
            "joint_rel_err": args.joint_rel_err,
            "max_joint_p": args.max_joint_p,
            "strong_joint_p": args.strong_joint_p,
        },
        "numeric": {
            "dps": args.dps,
            "quad_n": args.quad_n,
            "quad_n_list": args.quad_n_list,
        },
        "seed": args.seed,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="ROT RH / Xi-QED shared law preregistered lock audit")
    ap.add_argument("--dps", type=int, default=100)
    ap.add_argument("--alpha0-inv", type=float, default=137.035999177)
    ap.add_argument("--me-MeV", dest="me_MeV", type=float, default=0.51099895000)
    ap.add_argument("--mmu-MeV", dest="mmu_MeV", type=float, default=105.6583755)
    ap.add_argument("--mtau-MeV", dest="mtau_MeV", type=float, default=1776.86)
    ap.add_argument("--mpi-charged-MeV", dest="mpi_charged_MeV", type=float, default=139.57039)
    ap.add_argument("--quad-n", type=int, default=256)
    ap.add_argument("--quad-n-list", type=str, default="128,256,512")
    ap.add_argument("--random-delta-controls", type=int, default=2000)
    ap.add_argument("--random-delta-band", type=str, default="0.5,1.8")
    ap.add_argument("--pseudo-mass-trials", type=int, default=512)
    ap.add_argument("--pseudo-mass-families", type=str, default="ratio_jitter_e_pi_only,ratio_jitter_pi_only,ratio_jitter_all,preserve_mu_pi_break_epi,random_ratio_ordered_wide,affine_plus_small_ratio_noise,preserve_e_pi_break_mupi,affine_scale_all")
    ap.add_argument("--pseudo-mass-jitter", type=float, default=0.25)
    ap.add_argument("--local-mass-grid-span", type=float, default=0.30)
    ap.add_argument("--local-mass-grid-steps", type=int, default=41)
    ap.add_argument("--interesting-rel-err", type=float, default=1e-3)
    ap.add_argument("--strong-rel-err", type=float, default=1e-4)
    ap.add_argument("--joint-rel-err", type=float, default=1e-3)
    ap.add_argument("--max-joint-p", type=float, default=0.05)
    ap.add_argument("--strong-joint-p", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=1729)
    ap.add_argument("--protocol-only", action="store_true", help="Write preregistration protocol JSON/hash and exit before scoring")
    ap.add_argument("--protocol-file", type=str, default="", help="Optional protocol JSON to verify hash/use for run metadata")
    ap.add_argument("--out-prefix", type=str, default="xi_qed_shared_law_preregistered_lock")
    args = ap.parse_args()

    t0 = time.time()
    qlist = [int(x) for x in args.quad_n_list.split(",") if x.strip()]
    families = [x.strip() for x in args.pseudo_mass_families.split(",") if x.strip()]
    lo_s, hi_s = args.random_delta_band.split(",")
    random_delta_lo = float(lo_s)
    random_delta_hi = float(hi_s)

    protocol = build_protocol(args)
    protocol_hash = sha_json(protocol)
    protocol_path = f"{args.out_prefix}_lock_protocol.json"
    with open(protocol_path, "w", encoding="utf-8") as f:
        json.dump({"protocol_hash": protocol_hash, "protocol": protocol}, f, indent=2, sort_keys=True)

    line()
    print("ROT RH / XI-QED SHARED LAW — PREREGISTERED LOCK AUDIT")
    line()
    print(f"time                    : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"dps                     : {args.dps}")
    print(f"alpha(0)^-1             : {args.alpha0_inv:.15f}")
    print(f"masses MeV              : e={args.me_MeV}, mu={args.mmu_MeV}, tau={args.mtau_MeV}, pi_ch={args.mpi_charged_MeV}")
    print("frozen law              : Q1=(4/pi)m_pi, loop=e ; Q2=sqrt(2pi)*red(mu,pi), loop=e+mu")
    print("search                  : DISABLED — no coefficient menu, no formula search, no additive scales")
    print(f"protocol hash           : {protocol_hash}")
    print(f"out_prefix              : {args.out_prefix}")
    line("-")

    if args.protocol_file:
        with open(args.protocol_file, "r", encoding="utf-8") as f:
            loaded = json.load(f)
        loaded_hash = loaded.get("protocol_hash") or sha_json(loaded.get("protocol", loaded))
        print(f"loaded protocol file     : {args.protocol_file}")
        print(f"loaded protocol hash     : {loaded_hash}")
        if loaded_hash != protocol_hash:
            print("WARNING: loaded protocol hash differs from current command protocol hash")

    if args.protocol_only:
        print("Protocol-only mode: wrote protocol and exited before scoring.")
        print(f"protocol file            : {protocol_path}")
        line()
        return

    nodes_cache: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    base_masses = Masses(args.me_MeV, args.mmu_MeV, args.mtau_MeV, args.mpi_charged_MeV)

    line()
    print("COMPUTING XI CURVATURE COUNT")
    line()
    inv = compute_xi_invariants(args.dps)
    K2 = safe_float(inv["K2"])
    A0 = safe_float(inv["A0"])
    R4 = safe_float(inv["R4"])
    delta = args.alpha0_inv - A0
    rel_gap = delta / args.alpha0_inv
    print(f"K2 = -d2 log Xi(0)       : {fmt(K2, 21)}")
    print(f"A0 = 2*pi/K2             : {fmt(A0, 18)}")
    print(f"R4 = L4/K2^2             : {fmt(R4, 18)}")
    print(f"delta alpha^-1 - A0      : {fmt(delta, 18)}")
    print(f"relative gap             : {fmt(rel_gap, 13)}")

    line()
    print("FROZEN LAW SCORE")
    line()
    sc = score_frozen_law(delta, base_masses, args.quad_n, nodes_cache)
    print(f"Q1=(4/pi)m_pi            : {fmt(sc['Q1_pion_4overpi'], 15)} MeV")
    print(f"Q2=sqrt2pi*red(mu,pi)    : {fmt(sc['Q2_reduced_sqrt2pi'], 15)} MeV")
    print(f"Delta1 e-loop            : {fmt(sc['delta1_e_loop'], 15)} rel={sc['rel1_pion']:.6e}")
    print(f"Delta2 e+mu-loop         : {fmt(sc['delta2_emu_loop'], 15)} rel={sc['rel2_reduced']:.6e}")
    print(f"joint max/rms/mean rel   : {sc['joint_max_rel']:.6e} / {sc['joint_rms_rel']:.6e} / {sc['joint_mean_rel']:.6e}")
    print(f"prediction agreement rel : {sc['prediction_agreement_rel_to_delta']:.6e}")

    # forbidden diagnostic, printed but not used for pass/fail
    Q_forbid = (base_masses.mpi + base_masses.me) * 4.0 / math.pi
    d_forbid = delta_alpha_inv_spacelike_one_loop(Q_forbid, [base_masses.me], args.quad_n, nodes_cache)
    rel_forbid = abs(d_forbid - delta) / abs(delta)
    print(f"forbidden mixed diagnostic: Q=(m_pi+m_e)*4/pi = {fmt(Q_forbid, 15)} MeV rel={rel_forbid:.6e}")

    line()
    print("QUADRATURE STABILITY")
    line()
    stab_rows = []
    for qn in qlist:
        s2 = score_frozen_law(delta, base_masses, qn, nodes_cache)
        row = {"quad_n": qn, **s2}
        stab_rows.append(row)
        print(f"quad={qn:4d} rel1={s2['rel1_pion']:.12e} rel2={s2['rel2_reduced']:.12e} joint={s2['joint_max_rel']:.12e} agree={s2['prediction_agreement_rel_to_delta']:.12e}")
    stability_span = max(r["joint_max_rel"] for r in stab_rows) - min(r["joint_max_rel"] for r in stab_rows)

    # Random delta controls: frozen predicted deltas stay fixed; target delta varies.
    line()
    print("RANDOM-DELTA CONTROLS — FROZEN LAW ONLY")
    line()
    rng = random.Random(args.seed)
    rand_rows = []
    both_count = 0
    jointmax_count = 0
    jointrms_count = 0
    pred1 = sc["delta1_e_loop"]
    pred2 = sc["delta2_emu_loop"]
    real_joint = sc["joint_max_rel"]
    real_rms = sc["joint_rms_rel"]
    for i in range(args.random_delta_controls):
        dt = rng.uniform(random_delta_lo, random_delta_hi)
        r1 = abs(pred1 - dt) / abs(dt)
        r2 = abs(pred2 - dt) / abs(dt)
        jm = max(r1, r2)
        jr = math.sqrt(0.5 * (r1 * r1 + r2 * r2))
        both = (r1 <= sc["rel1_pion"] and r2 <= sc["rel2_reduced"])
        if both:
            both_count += 1
        if jm <= real_joint:
            jointmax_count += 1
        if jr <= real_rms:
            jointrms_count += 1
        rand_rows.append({"i": i, "delta_target": dt, "rel1": r1, "rel2": r2, "joint_max": jm, "joint_rms": jr, "both_individual_leq": int(both)})
        if (i + 1) % max(1, args.random_delta_controls // 8) == 0:
            print(f"  random-delta progress: {i+1}/{args.random_delta_controls}")
    random_both_p = empirical_p(both_count, args.random_delta_controls)
    random_jointmax_p = empirical_p(jointmax_count, args.random_delta_controls)
    random_jointrms_p = empirical_p(jointrms_count, args.random_delta_controls)
    print(f"random both-individual p : {random_both_p:.6g} ({both_count}/{args.random_delta_controls})")
    print(f"random joint-max p       : {random_jointmax_p:.6g} ({jointmax_count}/{args.random_delta_controls})")
    print(f"random joint-rms p       : {random_jointrms_p:.6g} ({jointrms_count}/{args.random_delta_controls})")

    # Pseudo controls
    line()
    print("PSEUDO-MASS CONTROLS — MAIN RATIO-CHANGING VS DIAGNOSTIC")
    line()
    pseudo_rows = []
    family_rows = []
    main_n = main_both = main_joint = main_rms = 0
    diag_n = diag_both = diag_joint = diag_rms = 0
    rng = random.Random(args.seed + 1)
    for fam in families:
        fam_rows = []
        fam_diag = None
        print(f"  pseudo family start: {fam}")
        for j in range(args.pseudo_mass_trials):
            pm, diag = generate_pseudo_masses(fam, base_masses, args.pseudo_mass_jitter, rng)
            fam_diag = diag
            ps = score_frozen_law(delta, pm, args.quad_n, nodes_cache)
            both = (ps["rel1_pion"] <= sc["rel1_pion"] and ps["rel2_reduced"] <= sc["rel2_reduced"])
            joint_leq = ps["joint_max_rel"] <= sc["joint_max_rel"]
            rms_leq = ps["joint_rms_rel"] <= sc["joint_rms_rel"]
            row = {
                "family": fam, "trial": j, "diagnostic": int(diag),
                "me": pm.me, "mmu": pm.mmu, "mtau": pm.mtau, "mpi": pm.mpi,
                **ps,
                "both_individual_leq": int(both),
                "joint_max_leq": int(joint_leq),
                "joint_rms_leq": int(rms_leq),
            }
            pseudo_rows.append(row)
            fam_rows.append(row)
            if diag:
                diag_n += 1; diag_both += int(both); diag_joint += int(joint_leq); diag_rms += int(rms_leq)
            else:
                main_n += 1; main_both += int(both); main_joint += int(joint_leq); main_rms += int(rms_leq)
            if (j + 1) % max(1, args.pseudo_mass_trials // 4) == 0:
                print(f"    {fam}: {j+1}/{args.pseudo_mass_trials}")
        n = len(fam_rows)
        fb = sum(r["both_individual_leq"] for r in fam_rows)
        fj = sum(r["joint_max_leq"] for r in fam_rows)
        fr = sum(r["joint_rms_leq"] for r in fam_rows)
        med_joint = float(np.median([r["joint_max_rel"] for r in fam_rows])) if fam_rows else float("nan")
        family_rows.append({
            "family": fam, "diagnostic": int(bool(fam_diag)), "n": n,
            "both_p": empirical_p(fb, n), "jointmax_p": empirical_p(fj, n), "jointrms_p": empirical_p(fr, n),
            "median_joint_max": med_joint,
        })
        print(f"{fam:35s} diag={bool(fam_diag)!s:5s} both_p={empirical_p(fb,n):.6g} jointmax_p={empirical_p(fj,n):.6g} med_joint={med_joint:.3e}")

    main_both_p = empirical_p(main_both, main_n) if main_n else float("nan")
    main_joint_p = empirical_p(main_joint, main_n) if main_n else float("nan")
    main_rms_p = empirical_p(main_rms, main_n) if main_n else float("nan")
    diag_both_p = empirical_p(diag_both, diag_n) if diag_n else float("nan")
    print("-" * 120)
    print(f"main pseudo rows          : {main_n}")
    print(f"main pseudo both p        : {main_both_p:.6g} ({main_both}/{main_n})")
    print(f"main pseudo joint-max p   : {main_joint_p:.6g} ({main_joint}/{main_n})")
    print(f"main pseudo joint-rms p   : {main_rms_p:.6g} ({main_rms}/{main_n})")
    print(f"diagnostic pseudo both p  : {diag_both_p:.6g}")

    # Local mass-ratio grid: perturb e/pi and mu/pi while anchored at pion scale.
    line()
    print("LOCAL MASS-RATIO GRID — JOINT LOCK")
    line()
    grid_rows = []
    steps = args.local_mass_grid_steps
    span = args.local_mass_grid_span
    shifts = np.linspace(-span, span, steps)
    grid_both = grid_joint = grid_rms = 0
    grid_off_n = grid_off_joint = 0
    for a in shifts:  # electron shift
        for b in shifts:  # muon shift
            pm = Masses(base_masses.me * math.exp(float(a)), base_masses.mmu * math.exp(float(b)), base_masses.mtau, base_masses.mpi)
            gs = score_frozen_law(delta, pm, args.quad_n, nodes_cache)
            both = (gs["rel1_pion"] <= sc["rel1_pion"] and gs["rel2_reduced"] <= sc["rel2_reduced"])
            joint_leq = gs["joint_max_rel"] <= sc["joint_max_rel"]
            rms_leq = gs["joint_rms_rel"] <= sc["joint_rms_rel"]
            grid_both += int(both)
            grid_joint += int(joint_leq)
            grid_rms += int(rms_leq)
            offdiag = not (abs(float(a)) < 1e-15 and abs(float(b)) < 1e-15)
            if offdiag:
                grid_off_n += 1
                grid_off_joint += int(joint_leq)
            grid_rows.append({"log_me_shift": float(a), "log_mmu_shift": float(b), **gs, "both_individual_leq": int(both), "joint_max_leq": int(joint_leq), "joint_rms_leq": int(rms_leq), "offdiag": int(offdiag)})
    grid_n = len(grid_rows)
    grid_both_p = empirical_p(grid_both, grid_n)
    grid_joint_p = empirical_p(grid_joint, grid_n)
    grid_rms_p = empirical_p(grid_rms, grid_n)
    grid_off_joint_p = empirical_p(grid_off_joint, grid_off_n) if grid_off_n else float("nan")
    grid_med = float(np.median([r["joint_max_rel"] for r in grid_rows]))
    grid_min = float(np.min([r["joint_max_rel"] for r in grid_rows]))
    print(f"grid steps               : {steps} x {steps} = {grid_n}")
    print(f"grid both p              : {grid_both_p:.6g} ({grid_both}/{grid_n})")
    print(f"grid joint-max p         : {grid_joint_p:.6g} ({grid_joint}/{grid_n})")
    print(f"grid joint-rms p         : {grid_rms_p:.6g} ({grid_rms}/{grid_n})")
    print(f"grid offdiag joint p     : {grid_off_joint_p:.6g} ({grid_off_joint}/{grid_off_n})")
    print(f"grid joint min/median    : {grid_min:.3e} / {grid_med:.3e}")

    # Pass/fail logic
    pass_rel = sc["joint_max_rel"] <= args.joint_rel_err
    pass_random = random_both_p <= args.max_joint_p and random_jointmax_p <= args.max_joint_p
    pass_pseudo = main_both_p <= args.max_joint_p and main_joint_p <= args.max_joint_p
    pass_grid = grid_joint_p <= args.max_joint_p and grid_off_joint_p <= args.max_joint_p
    pass_stability = stability_span <= 1e-7
    strong = (sc["joint_max_rel"] <= args.strong_rel_err and
              random_both_p <= args.strong_joint_p and
              main_both_p <= args.strong_joint_p and
              main_joint_p <= args.strong_joint_p and
              grid_joint_p <= args.strong_joint_p)

    if pass_rel and pass_random and pass_pseudo and pass_grid and pass_stability and strong:
        flag = "XI_QED_SHARED_LAW_PREREGISTERED_LOCK_HARDENED_STRONG"
    elif pass_rel and pass_random and pass_pseudo and pass_grid and pass_stability:
        flag = "XI_QED_SHARED_LAW_PREREGISTERED_LOCK_JOINT_INTERESTING"
    elif pass_rel and (pass_random or pass_pseudo):
        flag = "XI_QED_SHARED_LAW_PREREGISTERED_LOCK_WEAK_SIGNAL"
    else:
        flag = "XI_QED_SHARED_LAW_PREREGISTERED_LOCK_NOT_CONFIRMED"

    summary = {
        "global_flag": flag,
        "protocol_hash": protocol_hash,
        "K2": K2,
        "A0": A0,
        "R4": R4,
        "alpha0_inv": args.alpha0_inv,
        "delta_needed": delta,
        **sc,
        "forbidden_mixed_Q": Q_forbid,
        "forbidden_mixed_rel": rel_forbid,
        "stability_span": stability_span,
        "random_both_p": random_both_p,
        "random_jointmax_p": random_jointmax_p,
        "random_jointrms_p": random_jointrms_p,
        "main_pseudo_both_p": main_both_p,
        "main_pseudo_jointmax_p": main_joint_p,
        "main_pseudo_jointrms_p": main_rms_p,
        "diagnostic_pseudo_both_p": diag_both_p,
        "local_grid_both_p": grid_both_p,
        "local_grid_jointmax_p": grid_joint_p,
        "local_grid_jointrms_p": grid_rms_p,
        "local_grid_offdiag_jointmax_p": grid_off_joint_p,
        "passes": {
            "rel": pass_rel,
            "random": pass_random,
            "pseudo": pass_pseudo,
            "grid": pass_grid,
            "stability": pass_stability,
            "strong": strong,
        },
        "runtime_sec": time.time() - t0,
    }
    meta = {"protocol_hash": protocol_hash, "protocol": protocol, "summary": summary}

    prefix = args.out_prefix
    paths = {
        "protocol": protocol_path,
        "stability": f"{prefix}_pll_stability.csv",
        "random_delta": f"{prefix}_pll_random_delta_controls.csv",
        "pseudo_controls": f"{prefix}_pll_pseudo_controls.csv",
        "pseudo_family": f"{prefix}_pll_pseudo_family_summary.csv",
        "local_grid": f"{prefix}_pll_local_mass_grid.csv",
        "summary": f"{prefix}_pll_summary.json",
        "meta": f"{prefix}_pll_meta.json",
    }
    write_csv(paths["stability"], stab_rows)
    write_csv(paths["random_delta"], rand_rows)
    write_csv(paths["pseudo_controls"], pseudo_rows)
    write_csv(paths["pseudo_family"], family_rows)
    write_csv(paths["local_grid"], grid_rows)
    with open(paths["summary"], "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    with open(paths["meta"], "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)

    line()
    print("XI-QED SHARED LAW PREREGISTERED LOCK SUMMARY")
    line()
    print(f"global flag                    : {flag}")
    print(f"protocol hash                  : {protocol_hash}")
    print(f"A0                             : {fmt(A0, 18)}")
    print(f"delta needed                   : {fmt(delta, 18)}")
    print(f"Q1,Q2                          : {fmt(sc['Q1_pion_4overpi'], 12)} MeV, {fmt(sc['Q2_reduced_sqrt2pi'], 12)} MeV")
    print(f"rel1, rel2                     : {sc['rel1_pion']:.12e}, {sc['rel2_reduced']:.12e}")
    print(f"joint max/rms                  : {sc['joint_max_rel']:.12e}, {sc['joint_rms_rel']:.12e}")
    print(f"prediction agreement           : {sc['prediction_agreement_rel_to_delta']:.12e}")
    print(f"random both / jointmax p        : {random_both_p:.6g} / {random_jointmax_p:.6g}")
    print(f"main pseudo both / jointmax p   : {main_both_p:.6g} / {main_joint_p:.6g}")
    print(f"local grid jointmax / offdiag p : {grid_joint_p:.6g} / {grid_off_joint_p:.6g}")
    print(f"forbidden mixed diagnostic rel  : {rel_forbid:.12e}")
    print(f"stability span                  : {stability_span:.3e}")
    print(f"passes rel/random/pseudo/grid/stability/strong: {pass_rel} / {pass_random} / {pass_pseudo} / {pass_grid} / {pass_stability} / {strong}")
    print("-" * 120)
    print("Files written")
    for p in paths.values():
        print(f"  {p}")
    line()
    print(f"Total runtime: {time.time() - t0:.2f}s")
    print(f"AUDIT FLAG: {flag}")


if __name__ == "__main__":
    main()
