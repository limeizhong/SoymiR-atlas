#!/usr/bin/env python3
"""Compose final Figure 2 from annotation-category and confidence panels."""

from pathlib import Path
import os
import subprocess

from PIL import Image, ImageChops
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_MODULE = SCRIPT_DIR.parents[1]

PANEL_A_SCRIPT = FIGURE_MODULE / "scripts/Figure2/plot_figure2_from_supplementary_tables.py"
PANEL_A = FIGURE_MODULE / "results/intermediate_figures/misc/miRNA_annotation_source_annotation_category_100pct_bar.png"
PANEL_B_SCRIPT = FIGURE_MODULE / "scripts/Figure2/plot_status_star_distribution.py"
PANEL_B = FIGURE_MODULE / "results/intermediate_figures/misc/status_star_distribution.png"
OUTPUT = FIGURE_MODULE / "results/final_figures/Figure_2.png"

def white_background(image: Image.Image) -> Image.Image:
    image = image.convert("RGBA")
    background = Image.new("RGBA", image.size, "white")
    background.alpha_composite(image)
    return background.convert("RGB")


def resize_to_width(image: Image.Image, width: int) -> Image.Image:
    if image.width == width:
        return image
    height = round(image.height * width / image.width)
    return image.resize((width, height), Image.Resampling.LANCZOS)


def trim_white_border(image: Image.Image, padding: int = 55) -> Image.Image:
    """Trim outer white margins while preserving a small padding."""
    image = image.convert("RGB")
    diff = ImageChops.difference(image, Image.new("RGB", image.size, "white"))
    bbox = diff.getbbox()
    if bbox is None:
        return image
    left = max(bbox[0] - padding, 0)
    top = max(bbox[1] - padding, 0)
    right = min(bbox[2] + padding, image.width)
    bottom = min(bbox[3] + padding, image.height)
    return image.crop((left, top, right, bottom))


def detect_x_axis_range(image: Image.Image) -> tuple[int, int]:
    """Return x-start and x-end for the longest dark horizontal axis line."""
    arr = np.asarray(image.convert("RGB"))
    dark = (arr[:, :, 0] < 45) & (arr[:, :, 1] < 45) & (arr[:, :, 2] < 45)
    best = (0, image.width - 1, 0)
    for y in range(dark.shape[0]):
        xs = np.where(dark[y])[0]
        if len(xs) < image.width * 0.15:
            continue
        span = int(xs.max() - xs.min())
        if span > best[2]:
            best = (int(xs.min()), int(xs.max()), span)
    return best[0], best[1]


def align_panel_x_axis(panel: Image.Image, target_width: int, target_axis: tuple[int, int]) -> Image.Image:
    """Scale and place a panel so its x-axis endpoints match target_axis."""
    source_axis = detect_x_axis_range(panel)
    source_span = max(source_axis[1] - source_axis[0], 1)
    target_span = max(target_axis[1] - target_axis[0], 1)
    scale = target_span / source_span
    scaled_w = round(panel.width * scale)
    scaled_h = round(panel.height * scale)
    scaled = panel.resize((scaled_w, scaled_h), Image.Resampling.LANCZOS)
    offset_x = round(target_axis[0] - source_axis[0] * scale)
    out = Image.new("RGB", (target_width, scaled_h), "white")
    out.paste(scaled, (offset_x, 0))
    return out


def main() -> None:
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(FIGURE_MODULE / ".matplotlib_cache"))
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    # Sync Panel B data with latest annotation workflow before plotting
    UPDATE_SCRIPT = FIGURE_MODULE / "scripts/Figure2/update_confidence_star_data.py"
    subprocess.run(["python3", str(UPDATE_SCRIPT)], check=True, env=env)
    subprocess.run(["python3", str(PANEL_A_SCRIPT)], check=True, env=env)
    subprocess.run(["python3", str(PANEL_B_SCRIPT)], check=True, env=env)

    panel_a = trim_white_border(white_background(Image.open(PANEL_A)))
    panel_b = trim_white_border(white_background(Image.open(PANEL_B)))
    panel_b = align_panel_x_axis(panel_b, panel_a.width, detect_x_axis_range(panel_a))

    margin_x = 180
    margin_top = 80
    panel_gap = 80
    margin_bottom = 100
    canvas_w = panel_a.width + margin_x * 2
    canvas_h = margin_top + panel_a.height + panel_gap + panel_b.height + margin_bottom

    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")

    x0 = margin_x
    y_a = margin_top
    y_b = y_a + panel_a.height + panel_gap

    canvas.paste(panel_a, (x0, y_a))
    canvas.paste(panel_b, (x0, y_b))

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(OUTPUT, dpi=(600, 600))
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
