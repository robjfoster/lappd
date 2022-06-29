import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import (RBFInterpolator, RectBivariateSpline, bisplev,
                               bisplrep, interp2d)
from scipy.ndimage import gaussian_filter
from skimage.feature import (blob_dog, blob_doh, blob_log, match_template,
                             peak_local_max)
from sklearn.cluster import MeanShift


def detect_peaks(matrix, threshold=5, min_distance=10, exclude_border=False):
    # neighborhood = generate_binary_structure(2, 2)
    # image_max = maximum_filter(matrix, size=2, mode="constant")
    coordinates = peak_local_max(
        matrix, threshold_abs=threshold, min_distance=min_distance, exclude_border=exclude_border)
    return coordinates


def match(matrix, template):
    result = match_template(matrix, template)
    ij = np.unravel_index(np.argmax(result), result.shape)
    x, y = ij[::-1]
    plt.imshow(matrix, aspect="auto")
    plt.scatter(x, y, c="red")
    plt.show()


def cluster(matrix):
    # bandwidth = estimate_bandwidth(leftmatrix, quantile=0.2, n_samples=500)
    ms = MeanShift()
    ms.fit(leftmatrix)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    print("number of estimated clusters : %d" % n_clusters_)


def find_blobs(matrix, threshold=5):
    print("Finding peaks")
    blobs = blob_log(matrix, threshold=threshold)
    if len(blobs) > 0:
        print(blobs)
        plt.imshow(matrix, interpolation='none',
                   aspect='auto', cmap=cm.coolwarm)
        for blob in blobs:
            plt.scatter(blob[1], blob[2], c="red")
        plt.show()
    return blobs


def interp_matrix(matrix, startx=0, stopx=1014, starty=0, stopy=28, interpfactor=10):
    step = 1
    x = np.arange(startx, stopx, step)
    y = np.arange(starty, stopy, step)
    newx = np.arange(startx, stopx-1, step/interpfactor)
    newy = np.arange(starty, stopy-1, step/interpfactor)
    fn2 = interp2d(x, y, matrix, kind="cubic")
    interped = fn2(newx, newy)
    return interped


def gauss_interp(matrix, startx=0, stopx=1014, starty=0, stopy=28, interpfactor=10):
    step = 1
    x = np.arange(startx, stopx, step)
    y = np.arange(starty, stopy, step)
    newx = np.arange(startx, stopx-1, step/interpfactor)
    newy = np.arange(starty, stopy-1, step/interpfactor)
    xx, yy = np.meshgrid(x, y)
    newxx, newyy = np.meshgrid(newx, newy)
    sparse_grid = np.stack([xx.ravel(), yy.ravel()], -1)
    new_grid = np.stack([newxx.ravel(), newyy.ravel()], -1)
    gausfn = RBFInterpolator(sparse_grid, matrix.ravel(),
                             smoothing=0.5, kernel="gaussian", neighbors=10, epsilon=2.0)
    interped = gausfn(new_grid).reshape(newxx.shape)
    return interped


def rbs_interp(matrix, startx=0, stopx=1014, starty=0, stopy=28, interpfactor=10):
    step = 1
    x = np.arange(startx, stopx, step)
    y = np.arange(starty, stopy, step)
    newx = np.arange(startx, stopx-1, step/interpfactor)
    newy = np.arange(starty, stopy-1, step/interpfactor)
    xx, yy = np.meshgrid(x, y)
    newxx, newyy = np.meshgrid(newx, newy)
    rbsfn = RectBivariateSpline(y, x, matrix, s=100)
    interped = rbsfn(newy, newx)
    return interped


def bispl_interp(matrix, startx=0, stopx=1014, starty=0, stopy=28, interpfactor=10):
    step = 1
    x = np.arange(startx, stopx, step)
    y = np.arange(starty, stopy, step)
    z = matrix.ravel()
    xx, yy = np.meshgrid(x, y)
    bisplfn = bisplrep(yy.ravel(), xx.ravel(), z, kx=3, ky=3, s=2700.0)
    newx = np.arange(startx, stopx-1, step/interpfactor)
    newy = np.arange(starty, stopy-1, step/interpfactor)
    return bisplev(newy, newx, bisplfn).T


def qspline_interp(matrix, startx=0, stopx=1014, starty=0, stopy=28, interpfactor=10):
    from scipy.signal import cspline2d
    newsignal = cspline2d(matrix, 3)
    return newsignal


def image_gaussian(matrix):
    image = gaussian_filter(matrix, 1)
    import matplotlib.pyplot as plt
    plt.imshow(image, aspect="auto")
    plt.show()
    breakpoint()
    return image


def gaussMM(matrix, hits):
    from sklearn.mixture import GaussianMixture as GM
    gmm = GM(n_components=len(hits))
    newdata = matrix.reshape(matrix.shape[0] * matrix.shape[1], 1)
    gmm.fit(newdata)
    cluster = gmm.predict(newdata)
    cluster = cluster.reshape(matrix.shape[0], matrix.shape[1])
    plt.imshow(gmm.means_[cluster], aspect="auto")
    plt.show()
    breakpoint()


def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x, y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp(- (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                                     + c*((y-yo)**2)))
    return g.ravel()


if __name__ == "__main__":
    import sys

    import matplotlib.pyplot as plt
    from matplotlib import cm

    from ..lappdevent import LAPPDEvent

    base_dir = sys.argv[1]
    stripnumber = sys.argv[2]
    for levent, lpulses in LAPPDEvent.search_strip(stripnumber, base_dir):
        leftmatrix = levent.leftmatrix * -1
        rightmatrix = levent.rightmatrix * -1
        leftinterp = interp_matrix(leftmatrix)
        rightinterp = interp_matrix(rightmatrix)
        leftpeaks = detect_peaks2(leftinterp)
        rightpeaks = detect_peaks2(rightinterp)
        leftpeaks2 = detect_peaks2(leftinterp, threshold=10)
        rightpeaks2 = detect_peaks2(rightinterp, threshold=10)
        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(leftinterp, aspect="auto")
        ax[1].imshow(rightinterp, aspect="auto")
        ax[0].set_title("Left")
        ax[1].set_title("Right")
        # plt.imshow(leftinterp, aspect="auto")
        for peak in leftpeaks:
            ax[0].scatter(peak[1], peak[0], c="red", s=5)
        for peak in leftpeaks2:
            ax[0].scatter(peak[1], peak[0], c="green", marker="x", s=5)
        for peak in rightpeaks:
            ax[1].scatter(peak[1], peak[0], c="red", s=5)
        for peak in rightpeaks2:
            ax[1].scatter(peak[1], peak[0], c="green", marker="x", s=5)
        plt.show()
        levent.plot_both()
        breakpoint()
        # blobs = find_blobs(newvals2)
        # newvals = fn(newXnewY).reshape(10130, -1).T
        # breakpoint()
        # plt.imshow(newvals, interpolation='nearest',
        #            aspect='auto', cmap=cm.coolwarm)
        # plt.show()
        # newx = np.arange(0, 1013, 0.1)
        # newy = np.arange(0, 27, 0.1)
        # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        # newnewx, newnewy = np.meshgrid(newx, newy)
        # surf = ax.plot_surface(newnewy, newnewx, leftinterp, cmap=cm.coolwarm,
        #                        rstride=10, cstride=10, edgecolor="none")
        # plt.show()
        # ax.scatter(newx[int(blobs[0][1])], newy[int(blobs[0][2])], c="red")
        # fig.colorbar(surf, shrink=0.5, aspect=5)
        # plt.show()
        # plt.imshow(leftmatrix, interpolation='linear',
        #            aspect='auto', cmap=cm.coolwarm)
        # plt.colorbar()
        # plt.show()
        # breakpoint()
