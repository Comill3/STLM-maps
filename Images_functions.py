"""Functions to create images of a grid with one line colored each time and to concatenate images horizontally or vertically."""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

mpl.rcParams["font.family"] = "Arial"


def create_grid(
    grid_size, image_size, step_linescan, missing_spectra, saving_path, SVG
):
    """Create 32 images of a grid with one line colored each time

    Args:
        grid_size (int): Number of points per line.
        image_size (int): Size of one line in nm.
        step_linescan (float): step between 2 spectra in nm.
        missing_spectra (list): List of missing spectra.
        saving_path (str): Saving path.
        SVG (bool): If True, save the images as SVG files, otherwise as PNG files.
    """

    for i in range(grid_size[0]):
        fig = plt.figure(1, figsize=[5, 5])
        ax = fig.add_subplot(111)

        # Return evenly spaced numbers over a specified interval
        X = np.linspace(
            0.5 * step_linescan, (grid_size[1] - 0.5) * step_linescan, grid_size[1]
        )

        for g in range(grid_size[1]):
            Y = np.empty(grid_size[1])
            Y.fill((g + 0.5) * step_linescan)

            if g == i:
                for ii in range(grid_size[1]):
                    couleur = ii / grid_size[1]
                    cmap = plt.cm.rainbow

                    X1 = (ii + 0.5) * step_linescan
                    Y1 = Y[0]
                    plt.plot(X1, Y1, ".r", markersize=6, color=cmap(couleur))
            else:
                plt.plot(X, Y, ".k", markersize=4)

        for k in missing_spectra:
            l = k // grid_size[1]
            c = k % grid_size[1]
            plt.plot(
                (c + 0.5) * step_linescan, (l + 0.5) * step_linescan, ".w", markersize=4
            )
        if SVG:
            savename = os.path.join(saving_path, f"line{i}.svg")
        else:
            savename = os.path.join(saving_path, f"line{i}.png")

        plt.rcParams.update({"font.size": 8})
        ax.xaxis.tick_top()
        ax.yaxis.tick_right()
        ax.xaxis.set_label_coords(0.5, 1.14)
        ax.yaxis.set_label_coords(1.12, 0.5)
        ax.set_aspect("equal", adjustable="box")
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=12)
        plt.xlabel("Position (nm)", fontsize=14)
        plt.ylabel("Position (nm)", fontsize=14, rotation=-90)
        plt.xlim(image_size, 0)
        plt.ylim(step_linescan * grid_size[0], 0)

        plt.savefig(savename, dpi=300, bbox_inches="tight")
        plt.show()


def get_concat_h(im1, im2, wm):
    """Concatenates two images horizontally with an optional width margin (wm) and returns the resulting image.

    Args:
        im1 (Image): 1/2 image to concatenate.
        im2 (Image): 2/2 image to concatenate.
        wm (float): Width margin, which specifies the gap or margin between the two concatenated images.

    Returns:
        Image: Resulting image.
    """
    dst = Image.new(
        "RGB", (im1.width + wm, max(im1.height, im2.height)), color=(255, 255, 255, 0)
    )
    dst.paste(im1, (0, 0))
    print(im1.width)
    dst.paste(im2, (im1.width + wm - im2.width, 0))
    return dst


def get_concat_v(im1, im2, hm):
    """Concatenates two images vertically with an optional height margin (hm) and returns the resulting image.

    Args:
        im1 (Image): 1/2 image to concatenate.
        im2 (Image): 2/2 image to concatenate.
        hm (float): Height margin, which specifies the gap or margin between the two concatenated images.

    Returns:
        Image: Resulting image.
    """
    dst = Image.new(
        "RGB", (max(im1.width, im2.width), im1.height + hm), color=(255, 255, 255, 0)
    )
    dst.paste(im1, (350, 0))
    dst.paste(im2, (0, im1.height + hm - im2.height))
    return dst


def join_images(data_path, grid_path, mode, number_line, name, grid_size):
    """Concatenate images horizontally or vertically.

    Args:
        data_path (str): Spectra images path.
        grid_path (str): Grid images path.
        mode (str): Concatenation mode ('horizontal' or 'vertical').
        number_line (_type_):
        name (_type_): _description_
        grid_size (_type_): _description_
    """

    foldcontents = os.listdir(data_path)
    foldername = "Gif_images_" + mode
    if foldername not in foldcontents:
        os.makedirs(os.path.join(data_path, foldername))
    gif_image_path = os.path.join(data_path, foldername)

    ensemble_images = []

    data_all = []
    grid_all = []

    wm = 0
    hm = 0

    for g in range(number_line):
        gridname = f"line{g}.png"
        dataname = name + f"_from{g * grid_size[1]}to{(g + 1) * grid_size[1]}.png"

        data = Image.open(os.path.join(data_path, dataname))
        grid = Image.open(os.path.join(grid_path, gridname))
        data_all.append(data)
        grid_all.append(grid)

        w = data.width
        h = data.height

        if w > wm:
            wm = w
        if h > hm:
            hm = h

    for i, _ in enumerate(data_all):
        # for i in range(len(data_all)):

        data = data_all[i]
        grid = grid_all[i]

        if mode == "horizontal":
            image = get_concat_h(grid, data, wm)
            print(image.size)
            (w, h) = image.size
            image = image.resize((w // 2, h // 2), Image.ANTIALIAS)
            print(image.size)
        elif mode == "vertical":
            image = get_concat_v(grid, data, hm)
            print(image.size)
            image = image.resize((1076, 1398), Image.ANTIALIAS)
            print(image.size)

        image.save(os.path.join(gif_image_path, name + f"_{i}.png"))

        ensemble_images.append(image)

    foldcontents = os.listdir(data_path)
    gif_foldername = "Gif_" + mode
    if gif_foldername not in foldcontents:
        os.makedirs(os.path.join(data_path, gif_foldername))
    gif_path = os.path.join(data_path, gif_foldername)

    gifname = name + "_gif.gif"

    ensemble_images[0].save(
        os.path.join(gif_path, gifname),
        save_all=True,
        append_images=ensemble_images[1:],
        duration=1000,
        loop=0,
    )
