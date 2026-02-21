#!/usr/bin/env python3
"""
简要分析并绘制 residual_history.dat 的收敛曲线。
用法: python plot_residual.py [data/residual_history.dat]
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    fpath = sys.argv[1] if len(sys.argv) > 1 else "data/residual_history.dat"
    # Robust load: handle inconsistent columns (e.g. corrupted/incomplete lines)
    rows = []
    with open(fpath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 5:
                try:
                    rows.append([float(p) for p in parts])
                except ValueError:
                    pass
    if not rows:
        print("Empty file or no data.")
        return
    data = np.array(rows)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    step, t, L1, L2, ratio = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]

    # 简要分析
    print("=== Residual History Analysis ===")
    print(f"Steps recorded: {len(step)}")
    print(f"Initial L1: {L1[0]:.6e}  L2: {L2[0]:.6e}  ratio: {ratio[0]:.6e}")
    print(f"Final   L1: {L1[-1]:.6e}  L2: {L2[-1]:.6e}  ratio: {ratio[-1]:.6e}")
    if len(step) > 1:
        decay = np.exp(np.log(L1[-1]/L1[0]) / (step[-1] - step[0])) if L1[0] > 1e-30 else 0
        print(f"L1 decay factor per step (approx): {decay:.4f}")

    # 绘图
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1, ax2 = axes

    ax1.semilogy(step, L1, "b-o", markersize=4, label="L1")
    ax1.semilogy(step, L2, "r-s", markersize=4, label="L2")
    ax1.set_ylabel("Residual norm")
    ax1.legend()
    ax1.grid(True, which="both", alpha=0.3)

    ax2.plot(step, ratio, "g-o", markersize=4)
    ax2.set_xlabel("Step")
    ax2.set_ylabel("L1 ratio (R/R0)")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out = fpath.replace(".dat", ".png")
    plt.savefig(out, dpi=150)
    print(f"Saved: {out}")
    plt.show()

if __name__ == "__main__":
    main()
