import os
import numpy as np
import matplotlib.pyplot as plt
import odl
from odl.applications import tomo

# absolute path, so that it does not change when the working directory changes
fig_dir = os.path.abspath('figures')
os.makedirs(fig_dir, exist_ok=True)


def save_figure(name):
    plt.savefig(os.path.join(fig_dir, f'{name}.svg'), bbox_inches='tight')


def show_images(images, titles, clim=(0.98, 1.06), filename=None):
    # The default window shows the inner ellipses, the skull (value 2) is saturated
    vmin, vmax = clim if clim is not None else (None, None)
    fig, axes = plt.subplots(1, len(images), figsize=(4 * len(images), 4))
    for ax, image, title in zip(np.atleast_1d(axes), images, titles):
        ax.imshow(image.asarray().T, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.axis('off')
    plt.tight_layout()
    if filename is not None:
        save_figure(filename)
    plt.show()


def show_slices(x, title, clim=None, filename=None):
    # axial, coronal and sagittal slices through the middle of the volume
    volume = x.asarray()
    nx, ny, nz = volume.shape
    slices = [volume[:, :, nz // 2], volume[:, ny // 2, :], volume[nx // 2, :, :]]
    vmin, vmax = clim if clim is not None else (None, None)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, image, name in zip(axes, slices, ['Axial', 'Coronal', 'Sagittal']):
        ax.imshow(image.T, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(name)
        ax.axis('off')
    fig.suptitle(title)
    plt.tight_layout()
    if filename is not None:
        save_figure(filename)
    plt.show()


def axial_slice(x):
    # middle axial slice of a volume, as an element of the corresponding 2D space
    space = x.space
    slice_space = odl.uniform_discr(space.min_pt[:2], space.max_pt[:2], space.shape[:2], dtype='float32')
    volume = x.asarray()
    return slice_space.element(volume[:, :, volume.shape[2] // 2])


def add_noise(data, relative_noise, seed=0):
    noise = odl.phantom.white_noise(data.space, seed=seed)
    noise = float(relative_noise * data.norm() / noise.norm()) * noise
    return data + noise, noise.norm()


def relative_error(x, x_true):
    return (x - x_true).norm() / x_true.norm()


def huber_reconstruction(ray_trafo, data, lam, eps, niter, op_norms, x_true=None, x0=None):
    lam = float(lam)
    grad = odl.Gradient(ray_trafo.domain)
    data_fit = 0.5 * odl.functionals.L2NormSquared(ray_trafo.range).translated(data) * ray_trafo
    regularizer = lam * odl.functionals.Huber(grad.range, eps) * grad
    objective = data_fit + regularizer

    ray_trafo_norm, grad_norm = op_norms
    L = ray_trafo_norm**2 + lam * grad_norm**2 / eps
    step = 1.9 / L

    errors = []
    def callback(x):
        if x_true is not None:
            errors.append(relative_error(x, x_true))

    x = ray_trafo.domain.zero() if x0 is None else x0.copy()
    odl.solvers.steepest_descent(objective, x, line_search=step, maxiter=niter, callback=callback)
    return x, errors


def parker_weights(ray_trafo):
    weights = tomo.parker_weighting(ray_trafo)
    if ray_trafo.domain.ndim == 2:
        # In 2D, ODL's Parker weights are mirrored along the detector axis, so we flip them back
        weights = ray_trafo.range.element(np.flip(weights.asarray(), axis=1).copy())
    return weights
