import numpy as np
from skimage.restoration import wiener, richardson_lucy
from scipy.signal import fftconvolve


def wiener_deconv(matrix, psf, balance=0.1):
    deconved = wiener(matrix, psf, balance)
    return deconved


def do_wiener(matrix: np.ndarray, psf: np.ndarray, balance: float = 0.1, plot=False):
    deconved = wiener(matrix, psf, 1.0/(0.63**2 / np.std(matrix)**2))
    #deconved = richardson_lucy(matrix, psf, num_iter=25)
    if plot:
        import matplotlib.pyplot as plt
        plt.imshow(deconved, aspect="auto")
        plt.show()
    return deconved * matrix


if __name__ == "__main__":
    import sys

    import matplotlib.pyplot as plt
    from matplotlib import cm

    from ..lappdevent import LAPPDEvent
    from .interpolation import interp_matrix

    fancy_plot = False
    template = np.load("template.npy")
    template2d = np.zeros((28, 1014))
    template2d[13, 492:522] = template * 0.15 * 0.15
    template2d[14, 492:522] = template * 0.15
    template2d[15, 492:522] = template * 0.15 * 0.15
    template2d *= -1
    print(template)
    breakpoint()
    newx = np.arange(0, 1013, 0.1)
    newy = np.arange(0, 27, 0.1)
    newnewx, newnewy = np.meshgrid(newx, newy)
    base_dir = sys.argv[1]
    stripnumber = int(sys.argv[2])

    for levent in LAPPDEvent.itr_all(base_dir):
        if levent.event_no < 101:
            continue
        left = levent.leftmatrix * -1
        deconved = wiener_deconv(left, template2d, balance=0.1)
        interpleft = interp_matrix(left)
        interpdeconv = interp_matrix(deconved)
        if fancy_plot:
            fig = plt.figure()
            ax0 = fig.add_subplot(331)
            ax1 = fig.add_subplot(332)
            ax2 = fig.add_subplot(333, projection="3d")
            ax3 = fig.add_subplot(334)
            ax4 = fig.add_subplot(335)
            ax5 = fig.add_subplot(336, projection="3d")
            ax6 = fig.add_subplot(337)
            ax7 = fig.add_subplot(338)
            ax8 = fig.add_subplot(339, projection="3d")
            ax0.imshow(left, aspect="auto", interpolation="none")
            ax1.imshow(interpleft, aspect="auto", interpolation="none")
            ax2.plot_surface(newnewy, newnewx, interpleft, cmap=cm.coolwarm,
                             rstride=10, cstride=10, edgecolor="none")
            ax3.imshow(deconved, aspect="auto", interpolation="none")
            ax4.imshow(interpdeconv, aspect="auto", interpolation="none")
            ax5.plot_surface(newnewy, newnewx, interpdeconv, cmap=cm.coolwarm,
                             rstride=10, cstride=10, edgecolor="none")
            ax6.imshow(np.where(deconved < -0.0, -deconved *
                                left, 0), aspect="auto", interpolation="none")
            ax7.imshow(np.where(interpdeconv < -0.0, -interpdeconv *
                                interpleft, 0), aspect="auto", interpolation="none")
            ax8.plot_surface(newnewy, newnewx, np.where(interpdeconv < -0.0, -interpdeconv *
                                                        interpleft, 0), cmap=cm.coolwarm,
                             rstride=10, cstride=10, edgecolor="none")
            ax0.set_title("Original")
            ax1.set_title("Original, interpolated")
            ax2.set_title("Original, interpolated")
            ax3.set_title("Wiener deconvolution")
            ax4.set_title("Wiener deconvolution, interpolated")
            ax5.set_title("Wiener deconvolution, interpolated")
            ax6.set_title("Inversed")
            ax7.set_title("Inversed, interpolated")
            ax8.set_title("Inversed, interpolated")
        else:
            fig, ax = plt.subplots(2, 3)
            base_img = ax[0, 0].imshow(
                left, aspect="auto", interpolation="none")
            ax[0, 1].imshow(deconved, aspect="auto", interpolation="none")
            # ax[0, 2].imshow(np.where(deconved < -0.0, -deconved * left, 0),
            #                 aspect="auto", interpolation="none")
            ax[0, 2].imshow(-deconved * left,
                            aspect="auto", interpolation="none")
            base_deconv = ax[1, 0].imshow(
                interpleft, aspect="auto", interpolation="none")
            ax[1, 1].imshow(interpdeconv, aspect="auto", interpolation="none")
            # ax[1, 2].imshow(np.where(interpdeconv < -0.0, -interpdeconv *
            #                          interpleft, 0), aspect="auto", interpolation="none")
            ax[1, 2].imshow(-interpdeconv * interpleft,
                            aspect="auto", interpolation="none")
            ax[0, 0].set_title("Signal")
            ax[0, 1].set_title("Deconvolved")
            ax[0, 2].set_title("Inversed")
            ax[1, 0].set_title("Signal, interpolated")
            ax[1, 1].set_title("Deconvoled, interpolated")
            ax[1, 2].set_title("Inversed, interpolated")
            fig.colorbar(base_deconv, ax=ax.ravel().tolist(),
                         orientation="horizontal", fraction=0.046)
        plt.show()
        breakpoint()
