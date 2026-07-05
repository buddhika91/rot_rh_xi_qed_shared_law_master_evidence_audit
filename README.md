# ROT RH / Xi-QED Shared-Law Preregistered Lock Audit

This repository documents a numerical audit of a **preregistered, search-free two-scale law** connecting a canonical curvature invariant of the completed Riemann Xi function to a QED vacuum-polarization correction near the observed fine-structure constant.

The central object is not claimed to be a proof of the Riemann Hypothesis, nor a derivation of the fine-structure constant. The result should be read as a hardened numerical signal candidate: a frozen two-scale law survives random-delta, pseudo-mass, local-grid, and quadrature-stability controls.

---

## 1. Core idea

Let \(\Xi(s)\) denote the completed Riemann Xi function. On the critical line, define the curvature invariant

```math
K_2 = -\left.\frac{d^2}{d\gamma^2}\log\Xi\left(\frac12+i\gamma\right)\right|_{\gamma=0}.
```

From this curvature, define the Xi curvature count

```math
A_0 = \frac{2\pi}{K_2}.
```

Numerically, the audit obtains

```math
K_2 = 0.0462099862308379405462,
```

and therefore

```math
A_0 = 135.970291698259416.
```

The observed low-energy inverse fine-structure constant is taken as

```math
\alpha(0)^{-1}=137.035999177.
```

So the gap being studied is

```math
\Delta = \alpha(0)^{-1} - A_0.
```

Numerically,

```math
\Delta = 1.0657074787405918.
```

The question tested in this repository is:

> Can this gap be reproduced by a frozen, physically constrained QED vacuum-polarization law, without coefficient search, formula search, or additive mixed-mass tuning?

---

## 2. The preregistered two-scale law

The audit freezes two homogeneous external momentum scales.

The first scale is a pion projection scale:

```math
Q_1 = \frac{4}{\pi}m_\pi.
```

The second scale is a reduced muon-pion scale:

```math
Q_2 = \sqrt{2\pi}\,\frac{m_\mu m_\pi}{m_\mu+m_\pi}.
```

Using the charged pion and muon masses supplied in the run,

```math
m_\pi = 139.57039\ \mathrm{MeV},
```

```math
m_\mu = 105.6583755\ \mathrm{MeV},
```

these become

```math
Q_1 = 177.706539822\ \mathrm{MeV},
```

```math
Q_2 = 150.735568583\ \mathrm{MeV}.
```

The locked law then evaluates exact one-loop QED vacuum-polarization shifts:

```math
\Delta_1 = \Delta\alpha_e^{-1}(Q_1),
```

using the electron loop only, and

```math
\Delta_2 = \Delta\alpha_{e+\mu}^{-1}(Q_2),
```

using electron plus muon loops.

The audit compares both predictions to the Xi-alpha gap

```math
\Delta = \alpha(0)^{-1}-\frac{2\pi}{K_2}.
```

---

## 3. What is deliberately disabled

The main script is intentionally restrictive.

It disables:

- coefficient search;
- formula grammar search;
- additive mixed-mass external scales;
- arbitrary zero-shell corrections;
- direct fitting of the fine-structure constant;
- use of Riemann zeros in the QED matching law.

The external scales are frozen before scoring:

```math
Q_1 = \frac{4}{\pi}m_\pi,
```

```math
Q_2 = \sqrt{2\pi}\,\frac{m_\mu m_\pi}{m_\mu+m_\pi}.
```

The script still prints a forbidden diagnostic scale,

```math
Q_{\mathrm{forbidden}} = \frac{4}{\pi}(m_\pi+m_e),
```

because it is numerically closer, but it is not allowed to determine the pass/fail result. It is excluded because it uses additive mixed-mass tuning in the external scale.

---

## 4. Main audit result

The current run returns:

```text
AUDIT FLAG: XI_QED_SHARED_LAW_PREREGISTERED_LOCK_JOINT_INTERESTING
```

The protocol hash is:

```text
5b0b6d40a43635e7fc9c8a9d8d9c81b890f0d58e6de502c094f6519badf3c9ee
```

The frozen-law score is:

```math
Q_1 = 177.706539822\ \mathrm{MeV},
```

```math
Q_2 = 150.735568583\ \mathrm{MeV}.
```

The two predicted vacuum-polarization shifts are:

```math
\Delta_1 = 1.06489781962867,
```

```math
\Delta_2 = 1.06586167518193.
```

The target Xi-alpha gap is:

```math
\Delta = 1.0657074787405918.
```

Relative errors:

```math
\frac{|\Delta_1-\Delta|}{\Delta}=7.597386037697\times 10^{-4},
```

```math
\frac{|\Delta_2-\Delta|}{\Delta}=1.446892739479\times 10^{-4}.
```

Joint metrics:

```math
\max(\mathrm{rel}_1,\mathrm{rel}_2)=7.597386037697\times10^{-4},
```

```math
\mathrm{RMS}(\mathrm{rel}_1,\mathrm{rel}_2)=5.468718917871\times10^{-4}.
```

Prediction agreement:

```math
\frac{|\Delta_2-\Delta_1|}{\Delta}=9.044278777176\times10^{-4}.
```

---

## 5. Control tests

The script performs several controls to attack the result.

### 5.1 Random-delta controls

Random target gaps are sampled from the configured band. The frozen two-scale law is then tested against each random gap.

Reported result:

```text
random both / jointmax p : 0.00049975 / 0.00149925
```

### 5.2 Pseudo-mass controls

The script perturbs or replaces the physical mass hierarchy using several pseudo-mass families. Main ratio-changing pseudo-mass families are separated from diagnostic degenerate families.

Reported result:

```text
main pseudo both / jointmax p : 0.00065083 / 0.00715913
```

### 5.3 Local mass-grid controls

The script scans a local grid around the physical mass ratios and tests whether nearby perturbations beat the real point.

Reported result:

```text
local grid jointmax / offdiag p : 0.00118906 / 0.000594884
```

### 5.4 Quadrature stability

The audit checks numerical stability across Gauss-Legendre quadrature sizes:

```text
quad-n-list = 128,256,512
```

The reported stability span is:

```text
stability span : 1.038e-09
```

---

## 6. How to run

From PowerShell, run:

```powershell
python rot_rh_xi_qed_shared_law_preregistered_lock_audit.py `
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
  --interesting-rel-err 1e-3 `
  --strong-rel-err 1e-4 `
  --joint-rel-err 1e-3 `
  --max-joint-p 0.05 `
  --strong-joint-p 0.01 `
  --out-prefix xi_qed_shared_law_preregistered_lock
```

For a protocol-only run:

```powershell
python rot_rh_xi_qed_shared_law_preregistered_lock_audit.py `
  --protocol-only `
  --out-prefix xi_qed_shared_law_preregistered_lock_protocol
```

---

## 7. Output files

A typical run writes:

```text
xi_qed_shared_law_preregistered_lock_lock_protocol.json
xi_qed_shared_law_preregistered_lock_pll_stability.csv
xi_qed_shared_law_preregistered_lock_pll_random_delta_controls.csv
xi_qed_shared_law_preregistered_lock_pll_pseudo_controls.csv
xi_qed_shared_law_preregistered_lock_pll_pseudo_family_summary.csv
xi_qed_shared_law_preregistered_lock_pll_local_mass_grid.csv
xi_qed_shared_law_preregistered_lock_pll_summary.json
xi_qed_shared_law_preregistered_lock_pll_meta.json
```

---

## 8. Interpretation

The strongest careful statement supported by the current audit is:

> A preregistered, search-free two-scale law links the Xi curvature gap \(\alpha(0)^{-1}-2\pi/K_2\) to exact one-loop QED vacuum-polarization shifts at two homogeneous scales, \(Q_1=(4/\pi)m_\pi\) and \(Q_2=\sqrt{2\pi}\,m_\mu m_\pi/(m_\mu+m_\pi)\). The locked law survives random-delta, pseudo-mass, local-grid, and quadrature-stability controls as a joint-interesting signal candidate.

This should not be overstated. The audit does not prove that the Riemann Xi function derives the fine-structure constant. It shows that a specific, preregistered, search-free two-scale law produces a statistically interesting match to the Xi-alpha gap under the tested controls.

---

## 9. Negative results and excluded routes

This repository should be read together with the negative results that motivated the locked law.

The following routes are not claimed as successful:

- direct equality \(A_0=\alpha(0)^{-1}\);
- flexible zero-shell correction formulas;
- broad formula grammar discovery;
- additive mixed-mass external scales;
- unconstrained threshold naturalness searches;
- GUE/RH proof claims.

The surviving object is narrower:

```math
A_0 = \frac{2\pi}{K_2}
```

and the surviving signal candidate is the locked QED correction law:

```math
\alpha(0)^{-1}-A_0 \approx \Delta\alpha_e^{-1}\left(\frac{4}{\pi}m_\pi\right)
```

jointly with

```math
\alpha(0)^{-1}-A_0 \approx \Delta\alpha_{e+\mu}^{-1}\left(\sqrt{2\pi}\,\frac{m_\mu m_\pi}{m_\mu+m_\pi}\right).
```

---

## 10. Repository status

Current status:

```text
XI_QED_SHARED_LAW_PREREGISTERED_LOCK_JOINT_INTERESTING
```

Scientific status:

```text
Interesting numerical signal candidate; not a proof or derivation.
```

