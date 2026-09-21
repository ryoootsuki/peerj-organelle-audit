#!/usr/bin/env python3
"""Render a lightweight workflow overview without Graphviz."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


def main() -> None:
    out = Path("results/combined/figures")
    out.mkdir(parents=True, exist_ok=True)
    labels = [
        "Pinned manifests\n+ local FASTQ index",
        "Local references\n+ SeqKit cleanup",
        "Contig removal\n+ interval masking",
        "fastp + pair repair\n(run-level parallel)",
        "Independent mapping\nBAM-less MAPQ audit",
        "Subset selection\n(predefined rules)",
        "BAM/depth/variants\n(selected runs only)",
        "Tables, figures,\nrevision materials",
    ]
    fig, ax = plt.subplots(figsize=(15, 3.3))
    ax.set_xlim(0, len(labels) * 2)
    ax.set_ylim(0, 2)
    ax.axis("off")
    for i, label in enumerate(labels):
        x = i * 2 + 0.15
        box = FancyBboxPatch((x, 0.65), 1.55, 0.7, boxstyle="round,pad=0.04", linewidth=1.2, facecolor="white")
        ax.add_patch(box)
        ax.text(x + 0.775, 1.0, label, ha="center", va="center", fontsize=9)
        if i < len(labels) - 1:
            ax.add_patch(FancyArrowPatch((x + 1.56, 1.0), (x + 1.96, 1.0), arrowstyle="->", mutation_scale=12))
    ax.text(0.15, 1.65, "Same code path: Silene (nuclear + chloroplast + mitochondrion) and chicken (nuclear + mitochondrion)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "Figure_1_workflow_overview.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / "Figure_1_workflow_overview.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] workflow overview written: {out}")

if __name__ == "__main__":
    main()
