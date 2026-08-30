# -*- coding: utf-8 -*-
"""Loader / sanity check for the Quantum Arithmetic Benchmark Dataset v1.
Run:  python3 loader_demo.py
Loads each Hamiltonian, checks tridiagonal symmetry, computes eigenvalues and
basic spectral diagnostics (gap ratio, participation proxy).
"""
import json, os, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SLUGS = ["zeta_J50", "zeta_J100", "beta_J50", "delta_J22", "product_J100"]

def load(slug):
    d = os.path.join(HERE, slug)
    H = np.loadtxt(os.path.join(d, "hamiltonian.csv"), delimiter=",")
    meta = json.load(open(os.path.join(d, "metadata.json")))
    return H, meta

def gap_ratio(eigvals_gamma):
    g = np.sort(eigvals_gamma)
    d = np.diff(g)
    r = np.minimum(d[:-1], d[1:]) / np.maximum(d[:-1], d[1:])
    return float(np.mean(r))

def main():
    print(f"{'slug':14s} {'N':>4s} {'sym':>4s} {'trid':>5s} {'posdef':>7s} "
          f"{'<r> locked-like':>15s}")
    for slug in SLUGS:
        H, meta = load(slug)
        N = H.shape[0]
        sym = np.allclose(H, H.T, atol=1e-12)
        off = H - np.diag(np.diag(H)) - np.diag(np.diag(H, 1), 1) - np.diag(np.diag(H, -1), -1)
        trid = np.allclose(off, 0, atol=1e-12)
        ev = np.linalg.eigvalsh(H)
        posdef = ev.min() > 0
        # gamma = 1/sqrt(lambda); locked-like = the low-lying (largest lambda) band
        ev_desc = np.sort(ev)[::-1]
        n_lock = meta.get("locked_levels", max(1, N // 2))
        gammas = 1.0 / np.sqrt(np.maximum(ev_desc[:n_lock], 1e-300))
        r = gap_ratio(gammas) if n_lock >= 6 else float("nan")
        print(f"{slug:14s} {N:4d} {str(sym):>4s} {str(trid):>5s} {str(posdef):>7s} "
              f"{r:15.4f}   (Poisson .386 / GOE .536 / GUE .600)")
    print("\nMetadata example (zeta_J50):")
    _, meta = load("zeta_J50")
    print(json.dumps({k: v for k, v in meta.items()
                      if k not in ("alpha_diag", "b_offdiag")}, indent=1)[:900])

if __name__ == "__main__":
    main()
