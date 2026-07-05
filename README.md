# Xi-QED Shared Law: Master Evidence Audit

This repository contains a frozen, search-disabled numerical audit of a proposed connection between the Riemann Xi curvature scale and finite-mass QED vacuum polarization.

The central claim tested here is not that the fine-structure constant has been derived, and not that the Riemann Hypothesis has been proved. The claim is narrower:

> A frozen two-scale law links the Xi curvature alpha-gap to exact spacelike one-loop QED vacuum-polarization shifts at two homogeneous physical mass scales.

The current master audit returns:

```text
AUDIT FLAG: XI_QED_SHARED_LAW_MASTER_EVIDENCE_HARDENED_STRONG
```

---

## 1. Frozen law

The audit begins with the completed Riemann Xi curvature invariant

```math
K_2 = -\left.\frac{d^2}{d\gamma^2}\log \Xi\left(\frac{1}{2}+i\gamma\right)\right|_{\gamma=0}.
```

From this, define the Xi curvature count

```math
A_0 = \frac{2\pi}{K_2}.
```

The numerical value obtained in the master audit is

```text
K2 = 0.0462099862308379415778676
A0 = 135.970291698259424
```

The observed input value is

```text
alpha(0)^-1 = 137.035999177
```

so the gap to be explained is

```math
\Delta = \alpha(0)^{-1} - A_0.
```

The master audit gives

```text
delta needed = 1.06570747874057561
```

The frozen two-scale law is

```math
Q_1 = \frac{4}{\pi}m_\pi,
```

with an electron-loop vacuum-polarization correction, and

```math
Q_2 = \sqrt{2\pi}\,\frac{m_\mu m_\pi}{m_\mu + m_\pi},
```

with an electron+muon-loop vacuum-polarization correction.

The corresponding scales are

```text
Q1 = 177.706539822109 MeV
Q2 = 150.735568582767 MeV
Q2/m_mu = 1.426631517562632
```

---

## 2. Baseline numerical result

The exact spacelike one-loop QED vacuum-polarization calculation gives

```text
Delta1 e-loop    = 1.064897819629
Delta2 e+mu-loop = 1.065861675182
```

against the target

```text
delta needed = 1.06570747874057561
```

The relative errors are

```text
rel1 = 7.597386037544e-04
rel2 = 1.446892739632e-04
```

with joint metrics

```text
joint max rel = 7.597386037544e-04
joint rms rel = 5.468718917784e-04
prediction agreement rel = 9.044278777176e-04
```

---

## 3. Master evidence table

The master audit consolidates the frozen-law result with random, pseudo-mass, local-grid, scheme, threshold, and mass-ratio-origin controls.

Current master evidence summary:

```text
random both / jointmax p            : 0.00049975012 / 0.00049975012
main pseudo both / jointmax p       : 0.00097624471 / 0.0048812236
local grid jointmax / offdiag p     : 0.0011890606 / 0.000594884
scheme asymp degradation            : 169.434x
Q2 mu residual/exact                : 3.82618278439
Q2 exact->mu-asymp degradation      : 889.622x
threshold sensitivity Dlog/target   : 0.255762884551
rho scan / ratio2d / synthetic p    : 0.0099502488 / 0.0011890606 / 0.00049975012
target reduced-mass rank            : 1 / 9
stability span                      : 1.038e-09
```

All major pass conditions returned true:

```text
passes:
  baseline_rel             : True
  stability                : True
  random                   : True
  pseudo                   : True
  local_grid               : True
  scheme_exact_best        : True
  threshold_exact_kernel   : True
  threshold_sensitive      : True
  mass_ratio_origin        : True
```

---

## 4. What is disabled

The master script is intentionally not a formula search.

Disabled by design:

```text
coefficient search
formula grammar search
additive mixed-mass external Q scales
post-hoc replacement of Q1 or Q2
new physical constants
new particle masses
```

The law tested is frozen before the controls are run.

---

## 5. Main script

The main script is

```text
rot_rh_xi_qed_shared_law_master_evidence_audit.py
```

Recommended repository layout:

```text
.
├── README.md
├── scripts/
│   └── rot_rh_xi_qed_shared_law_master_evidence_audit.py
└── results/
    ├── xi_qed_shared_law_master_evidence_master_summary.json
    ├── xi_qed_shared_law_master_evidence_master_evidence_table.csv
    ├── xi_qed_shared_law_master_evidence_master_stability.csv
    ├── xi_qed_shared_law_master_evidence_master_random_delta_controls.csv
    ├── xi_qed_shared_law_master_evidence_master_pseudo_controls.csv
    ├── xi_qed_shared_law_master_evidence_master_scheme_rows.csv
    ├── xi_qed_shared_law_master_evidence_master_threshold_components.csv
    ├── xi_qed_shared_law_master_evidence_master_threshold_hybrids.csv
    ├── xi_qed_shared_law_master_evidence_master_rho_scan.csv
    ├── xi_qed_shared_law_master_evidence_master_ratio2d_scan.csv
    ├── xi_qed_shared_law_master_evidence_master_reduced_mass_diagnostics.csv
    └── xi_qed_shared_law_master_evidence_master_synthetic_ratio_controls.csv
```

---

## 6. How to run

From the repository root:

```powershell
python scripts/rot_rh_xi_qed_shared_law_master_evidence_audit.py `
  --dps 100 `
  --alpha0-inv 137.035999177 `
  --me-MeV 0.51099895000 `
  --mmu-MeV 105.6583755 `
  --mtau-MeV 1776.86 `
  --mpi-charged-MeV 139.57039 `
  --quad-n 256 `
  --quad-n-list "128,256,512" `
  --random-delta-controls 2000 `
  --random-delta-band "0.5,1.8" `
  --pseudo-mass-trials 512 `
  --pseudo-mass-families "ratio_jitter_e_pi_only,ratio_jitter_pi_only,ratio_jitter_all,preserve_mu_pi_break_epi,random_ratio_ordered_wide,affine_plus_small_ratio_noise,preserve_e_pi_break_mupi,affine_scale_all" `
  --pseudo-mass-jitter 0.25 `
  --local-mass-grid-span 0.30 `
  --local-mass-grid-steps 41 `
  --rho-scan-log-span 0.50 `
  --rho-scan-steps 401 `
  --ratio2d-log-span 0.35 `
  --ratio2d-steps 41 `
  --synthetic-ratio-trials 2000 `
  --synthetic-ratio-log-band 0.50 `
  --seed 314159 `
  --interesting-rel-err 1e-3 `
  --joint-rel-err 1e-3 `
  --max-joint-p 0.05 `
  --strong-joint-p 0.01 `
  --scheme-degradation-min 25 `
  --threshold-residual-ratio-min 0.05 `
  --sensitivity-min 0.02 `
  --out-prefix results/xi_qed_shared_law_master_evidence
```

Protocol-only run:

```powershell
python scripts/rot_rh_xi_qed_shared_law_master_evidence_audit.py `
  --protocol-only `
  --out-prefix results/xi_qed_shared_law_master_evidence_protocol
```

The master run used the protocol hash

```text
445d92a9f5cb7e6bea491646b8479876cebb6b4a0883581088b1f12169ef68fb
```

---

## 7. Evidence modules

### 7.1 Quadrature stability

The baseline result is stable under quadrature sizes 128, 256, and 512:

```text
quad=128 joint=7.597395953648e-04
quad=256 joint=7.597386037544e-04
quad=512 joint=7.597385571152e-04
```

The reported stability span is

```text
1.038e-09
```

### 7.2 Random-delta controls

Random target gaps were tested against the frozen law.

```text
random both-individual p = 0.00049975012
random joint-max p       = 0.00049975012
```

### 7.3 Pseudo-mass controls

Pseudo-mass families perturb or replace mass ratios while preserving selected structural features.

Main pseudo-control result:

```text
main pseudo both p      = 0.00097624471
main pseudo joint-max p = 0.0048812236
```

### 7.4 Local mass-ratio grid

The local two-dimensional mass-ratio grid scans around the observed electron/muon and pion/muon ratios.

```text
grid joint-max p     = 0.0011890606
grid offdiag joint p = 0.000594884
best grid shift      = observed point
```

### 7.5 Scheme robustness

The exact spacelike one-loop scheme is essential. The baseline exact-spacelike loop assignment ranks first among exact-spacelike loop variants:

```text
exact-spacelike baseline loop rank = 1 / 4
```

Asymptotic approximations degrade the match strongly:

```text
asymp -5/3 baseline joint ratio     = 169.434x
asymp log-only baseline joint ratio = 267.389x
```

This indicates that the signal is not explained by generic logarithmic running of alpha.

### 7.6 Threshold structure

The Q2 channel depends on the finite-mass muon threshold structure.

At Q2:

```text
Q2/m_mu = 1.426631517562632
```

The muon contribution decomposes as

```text
Q2 mu exact contribution      = 0.0358923745156
Q2 mu asymptotic contribution = -0.101438410947
Q2 mu finite-mass residual    = 0.137330785462
Q2 mu residual/exact          = 3.82618278439
```

Replacing the exact muon kernel by asymptotic approximations damages the result:

```text
Q2 exact->mu-asymp degradation = 889.622x
Q2 exact->mu-log degradation   = 257.219x
```

### 7.7 Threshold sensitivity

The threshold stationarity test found that Q2/m_mu is not a simple stationary point or inflection point of the muon kernel.

Instead, it is a sensitive finite-mass threshold point:

```text
near stationary / inflection      : False / False
threshold sensitivity Dlog/target : 0.255762884551
```

### 7.8 Mass-ratio origin

The physical pion/muon mass ratio is locally special under the tested controls.

```text
rho = m_pi/m_mu = 1.320959075317224
```

The scan results are

```text
rho scan p joint      = 0.0099502488
ratio2d p joint       = 0.0011890606
synthetic ratio p     = 0.00049975012
target reduced-mass rank = 1 / 9
```

The reduced-mass structure tested is

```math
\operatorname{red}(m_\mu,m_\pi) = \frac{m_\mu m_\pi}{m_\mu + m_\pi}.
```

---

## 8. Forbidden diagnostic

The audit also reports a forbidden diagnostic scale:

```math
Q_{\mathrm{forbidden}} = \frac{4}{\pi}(m_\pi + m_e).
```

It gives a closer numerical match:

```text
forbidden mixed diagnostic rel = 3.207225521622e-05
```

However, this scale is excluded from the frozen law because it uses additive mixed-mass external scale construction. It is kept only as a diagnostic and is not counted as part of the claim.

---

## 9. Interpretation

The strongest current interpretation is:

> The Xi curvature alpha-gap is matched by a frozen, homogeneous, two-scale law only when QED is evaluated with the exact spacelike one-loop finite-mass vacuum-polarization kernel. The Q2 channel specifically depends on the finite-mass muon threshold structure and the physical pion/muon reduced-mass ratio.

This is stronger and more precise than saying that generic alpha running explains the gap.

The result should be described as a hardened numerical signal candidate requiring analytic explanation.

---

## 10. Limitations

This repository does not claim:

```text
a proof of the Riemann Hypothesis
a derivation of alpha from first principles
a derivation of the pion or muon mass
a complete physical theory of the fine-structure constant
a replacement for Standard Model renormalization
```

The open analytic questions are:

```text
Why should the Xi curvature count A0 be corrected by QED vacuum polarization?
Why should the correction use Q1=(4/pi)m_pi?
Why should the second scale use sqrt(2pi)*red(m_mu,m_pi)?
Why does the exact spacelike finite-mass kernel appear rather than an asymptotic scheme?
Can the two-scale law be derived from a common geometric or spectral principle?
```

---

## 11. Reproducibility notes

The run used:

```text
alpha(0)^-1       = 137.035999177
m_e               = 0.51099895000 MeV
m_mu              = 105.6583755 MeV
m_tau             = 1776.86 MeV
m_pi_charged      = 139.57039 MeV
dps               = 100
quad-n            = 256
quad-n-list       = 128,256,512
random controls   = 2000
pseudo trials     = 512 per family
synthetic trials  = 2000
seed              = 314159
```

The master evidence run returned:

```text
AUDIT FLAG: XI_QED_SHARED_LAW_MASTER_EVIDENCE_HARDENED_STRONG
```

---

## 12. Recommended citation wording

If referencing this repository, use careful language:

> This repository reports a frozen, search-disabled numerical audit of a proposed Xi-QED shared law. The current master run finds a hardened numerical signal connecting the Xi curvature alpha-gap to exact spacelike one-loop QED vacuum polarization at two homogeneous mass scales. The result is not a proof or derivation, but it survives multiple numerical controls and motivates analytic investigation.
