from collections import namedtuple
from math import sqrt
from typing import List, Tuple

import numpy as np
from skimage.feature import peak_local_max

from .utils import lappdcfg as cfg
from .utils.pdf import AmplPDF, LocPDF, TimePDF


Hit = namedtuple("Hit", "x y")
Hit3D = namedtuple("Hit3D", "x y z")
Pair = namedtuple("Pair", "left right")

x = np.arange(0, 1014, 1)
y = np.arange(0, 8, 1)
xx, yy = np.meshgrid(x, y)
# newx = np.arange(0, 1013, 0.1)
# newy = np.arange(0, 27, 0.1)
# newnewx, newnewy = np.meshgrid(newx, newy)


timepdf = TimePDF(cfg.TIMEPDFSIGMA, cfg.TIMEPDFSTART, cfg.TIMEPDFEND, 1000)
locpdf = LocPDF(cfg.LOCPDFSIGMA, cfg.LOCPDFSTART, cfg.LOCPDFEND, 1000)
amplpdf = AmplPDF(cfg.AMPLPDFSIGMA, cfg.AMPLPFDFSTART, cfg.AMPLPDFEND, 1000)


def x_to_t(x, interpfactor=10.0):
    return (x / interpfactor) * cfg.NS_PER_SAMPLE


def y_to_loc(y, interpfactor=10.0):
    strip = y / interpfactor
    return 200 - ((cfg.STRIPWIDTH * (strip+1)) + cfg.STRIPSPACING)


def z_to_ampl(z):
    return z


def transverse_position(lefttime, righttime):
    timedelta = lefttime - righttime
    pos = timedelta / cfg.MAXDELTA * (cfg.STRIPLENGTH / 2.0)
    timedeltaerror = sqrt((cfg.TIMEJITTER)**2 + (cfg.TIMEJITTER)**2)
    poserror = pos * sqrt((timedeltaerror / timedelta) **
                          2 + (cfg.MAXDELTAERROR / cfg.MAXDELTA)**2)
    return pos, abs(poserror)


def detect_peaks(matrix, threshold=5, min_distance=10, exclude_border=False):
    # neighborhood = generate_binary_structure(2, 2)
    # image_max = maximum_filter(matrix, size=2, mode="constant")
    coordinates = peak_local_max(
        matrix, threshold_abs=threshold, min_distance=min_distance, exclude_border=exclude_border)
    peaks = []
    for coord in coordinates:
        peaks.append(Hit(coord[1], coord[0]))
    return peaks


def match_peaks(leftpeaks: np.ndarray,
                rightpeaks: np.ndarray,
                min_like: float = 1e-15
                ) -> List[Tuple[int]]:
    if len(leftpeaks) == 0 or len(rightpeaks) == 0:
        return [], [[], []]
    likelihood = []
    for lpeak in leftpeaks:
        # Convert x, y, z to t, loc, ampl
        lpeakt = x_to_t(lpeak.x)
        lpeakloc = y_to_loc(lpeak.y)
        lpeakampl = z_to_ampl(leftcleaninterp[lpeak.y, lpeak.x])
        thislike = []
        for rpeak in rightpeaks:
            rpeakt = x_to_t(rpeak.x)
            rpeakloc = y_to_loc(rpeak.y)
            rpeakampl = z_to_ampl(rightcleaninterp[rpeak.y, rpeak.x])
            # Calculate the left right deltas
            deltax = np.abs(lpeakt - rpeakt)
            deltay = np.abs(lpeakloc - rpeakloc)
            deltaz = np.abs(lpeakampl - rpeakampl)
            # Draw likelihood values from PDFs
            likex = timepdf.sample(deltax)
            likey = locpdf.sample(deltay)
            likez = amplpdf.sample(deltaz)
            # Total likelihood is product of each individual likelihood
            thislike.append(likex * likey * likez)
        likelihood.append(thislike)
    likelihood = np.asarray(likelihood)
    min_like = 1e-15
    # Copy of likelihood matrix that will be *modified* in the loop.
    # likelihood is left unmodified
    mod_likelihood = np.copy(likelihood)
    pairs = []
    while True:
        # Get index of max value in matrix
        pair = np.unravel_index(
            np.argmax(mod_likelihood), mod_likelihood.shape)
        # Find value of likelihood at this index
        maxlike = mod_likelihood[pair[0], pair[1]]
        # Compare to our likelihood limit
        if maxlike < min_like:
            print(f"Last likelihood value: {maxlike}")
            break
        # If it passes, keep it
        thispair = Pair(leftpeaks[pair[0]], rightpeaks[pair[1]])
        # pairs.append((leftpeaks[pair[1]], rightpeaks[pair[0]]))
        pairs.append(thispair)
        # Remove the row and column from contention because these have been paired off.
        mod_likelihood[pair[0], :] = -9999
        mod_likelihood[:, pair[1]] = -9999
    unmatched = np.nonzero(mod_likelihood > -1)
    left_unmatched = []
    for i in unmatched[0]:
        left_unmatched.append(leftpeaks[i])
    right_unmatched = []
    for i in unmatched[1]:
        right_unmatched.append(rightpeaks[i])
    return pairs, (left_unmatched, right_unmatched)


strip_positions = y_to_loc(np.arange(0, 28, 1), interpfactor=1)

if __name__ == "__main__":
    import sys

    import matplotlib.pyplot as plt
    from matplotlib import cm

    from .lappdevent import LAPPDEvent
    from .utils.interpolation import interp_matrix
    from .utils.wiener import do_wiener

    fancy_plot = True
    base_dir = sys.argv[1]
    stripnumber = int(sys.argv[2])
    for levent in LAPPDEvent.itr_all(base_dir):
        if levent.event_no < 301:
            continue
        if np.max(levent.leftmatrix) > 70 or np.max(levent.rightmatrix) > 70:
            continue
        template = np.load("template2d.npy")
        print(f"Event number: {levent.event_no}")
        left = levent.leftmatrix
        right = levent.rightmatrix
        thismin = min((np.min(left), np.min(right)))
        thismax = max((np.max(left), np.max(right)))
        leftclean = do_wiener(left, template)
        rightclean = do_wiener(right, template)
        # leftdeconinterp = interp_matrix(leftdecon)
        # rightdeconinterp = interp_matrix(rightdecon)
        # leftinterp = interp_matrix(left)
        # rightinterp = interp_matrix(right)
        # leftclean = leftdeconinterp * leftinterp
        # rightclean = rightdeconinterp * rightinterp
        # leftpeaks = detect_peaks2(leftclean)
        # rightpeaks = detect_peaks2(rightclean)
        leftcleaninterp = interp_matrix(leftclean)
        rightcleaninterp = interp_matrix(rightclean)
        leftpeaks = detect_peaks(
            leftcleaninterp, threshold=6, exclude_border=(3, 3))
        rightpeaks = detect_peaks(
            rightcleaninterp, threshold=6, exclude_border=(3, 3))
        # likelihood = []
        # if len(leftpeaks) == 0 or len(rightpeaks) == 0:
        #     continue
        # for lpeak in leftpeaks:
        #     lpeakt = x_to_t(lpeak[1])
        #     lpeakloc = y_to_loc(lpeak[0])
        #     lpeakampl = z_to_ampl(leftcleaninterp[lpeak[0], lpeak[1]])
        #     thislike = []
        #     for rpeak in rightpeaks:
        #         rpeakt = x_to_t(rpeak[1])
        #         rpeakloc = y_to_loc(rpeak[0])
        #         rpeakampl = z_to_ampl(rightcleaninterp[rpeak[0], rpeak[1]])
        #         deltax = np.abs(lpeakt - rpeakt)
        #         deltay = np.abs(lpeakloc - rpeakloc)
        #         deltaz = np.abs(lpeakampl - rpeakampl)
        #         likex = time_like(deltax, pdfsum=time_sum)
        #         likey = loc_like(deltay, pdfsum=loc_sum)
        #         likez = ampl_like(deltaz, pdfsum=ampl_sum)
        #         thislike.append(likex * likey * likez)
        #     likelihood.append(thislike)
        # likelihood = np.asarray(likelihood)
        # min_like = 1e-15
        # # Copy of likelihood matrix that will be *modified* in the loop.
        # # likelihood is left unmodified
        # mod_likelihood = np.copy(likelihood)
        # pairs = []
        # while True:
        #     # Get index of max value in matrix
        #     pair = np.unravel_index(
        #         np.argmax(mod_likelihood), mod_likelihood.shape)
        #     # Find value of likelihood at this index
        #     maxlike = mod_likelihood[pair[0], pair[1]]
        #     # Compare to our likelihood limit
        #     if maxlike < min_like:
        #         print(f"Last likelihood value: {maxlike}")
        #         break
        #     # If it passes, keep it
        #     pairs.append(pair)
        #     # Remove the row and column from contention because these have been paired off.
        #     mod_likelihood[pair[0], :] = -9999
        #     mod_likelihood[:, pair[1]] = -9999
        pairs, unmatched = match_peaks(leftpeaks, rightpeaks)
        # Plotting
        if fancy_plot:
            fig = plt.figure()
            fig.suptitle(f"Event {levent.event_no}")
            ax0 = fig.add_subplot(231, projection="3d")
            ax1 = fig.add_subplot(232)
            ax2 = fig.add_subplot(233)
            ax3 = fig.add_subplot(234, projection="3d")
            ax4 = fig.add_subplot(235)
            ax5 = fig.add_subplot(236)
            base_img = ax1.imshow(
                left[0:8], aspect="auto", interpolation="none", vmin=thismin, vmax=thismax)
            ax4.imshow(right[0:8], aspect="auto",
                       interpolation="none", vmin=thismin, vmax=thismax)
            ax2.imshow(leftcleaninterp[0:80],
                       aspect="auto", interpolation="none", vmin=thismin, vmax=thismax)
            ax5.imshow(rightcleaninterp[0:80],
                       aspect="auto", interpolation="none", vmin=thismin, vmax=thismax)
            surf1 = ax0.plot_surface(yy, xx, left[0:8],
                                     rstride=1, cstride=1, vmin=thismin, vmax=thismax, cmap=cm.coolwarm)
            surf2 = ax3.plot_surface(yy, xx, right[0:8],
                                     rstride=1, cstride=1, vmin=thismin, vmax=thismax, cmap=cm.coolwarm)
            ax0.set_title("Left, observed signal")
            ax3.set_title("Right, observed signal")
            ax1.set_title("Left, observed signal")
            ax4.set_title("Right, observed signal")
            ax2.set_title("Left, deconvolved, interpolated")
            ax5.set_title("Right, deconvolved, interpolated")
            # To set the axis ticks:
            # ax1.set_xticks(np.arange(0, 1014, 100))
            # ax1.set_xticklabels(np.arange(0, 210, 20))
            fig.colorbar(base_img, ax=[ax0, ax1, ax2, ax3, ax4, ax5],
                         orientation="horizontal", fraction=0.016)
        else:
            fig, ax = plt.subplots(2, 2)
            base_img = ax[0, 0].imshow(
                left[0:8], aspect="auto", interpolation="none")
            ax[0, 1].imshow(right[0:8], aspect="auto", interpolation="none")
            ax[1, 0].imshow(leftcleaninterp[0:80],
                            aspect="auto", interpolation="none")
            ax[1, 1].imshow(rightcleaninterp[0:80],
                            aspect="auto", interpolation="none")
            ax[0, 0].set_title("Left, observed signal")
            ax[0, 1].set_title("Right, observed signal")
            ax[1, 0].set_title("Left, deconvolved, interpolated")
            ax[1, 1].set_title("Right, deconvolved, interpolated")
            fig.colorbar(base_img, ax=ax.ravel().tolist(),
                         orientation="horizontal", fraction=0.046)
            fig.suptitle(f"Event {levent.event_no}")
            # plt.imshow(leftinterp, aspect="auto")
            for peak in leftpeaks:
                if peak[0] < 3:
                    continue
                ax[1, 0].scatter(peak[1], peak[0], c="red", s=5, marker="x")
            for peak in rightpeaks:
                if peak[0] < 3:
                    continue
                ax[1, 1].scatter(peak[1], peak[0], c="red", s=5, marker="x")
        hits = []
        hiterrs = []
        # Plot the matched peaks as pink circles and number them
        for i, pair in enumerate(pairs):
            # lpairx = leftpeaks[pair[0]][1]
            # lpairy = leftpeaks[pair[0]][0]
            # rpairx = rightpeaks[pair[1]][1]
            # rpairy = rightpeaks[pair[1]][0]
            xpos, xposerr = transverse_position(
                x_to_t(pair.left.x), x_to_t(pair.right.x))
            ypos = (y_to_loc(pair.left.y) + y_to_loc(pair.right.y)) / 2.0
            hits.append(Hit(xpos, ypos))
            hiterrs.append(Hit(xposerr, 5))
            ax2.scatter(pair.left.x, pair.left.y, c="pink", s=25, marker="o")
            ax5.scatter(pair.right.x, pair.right.y, c="pink", s=25, marker="o")
            # , c="green", s=4, marker="x")
            ax2.annotate(str(i), xy=(pair.left.x, pair.left.y), c="white")
            # , c="green", s=4, marker="x")
            ax5.annotate(str(i), xy=(pair.right.x, pair.right.y), c="white")
        lpair = [i[0] for i in pairs]
        rpair = [i[1] for i in pairs]
        # Plot the unmatched peaks as red crosses
        for lhit in unmatched[0]:
            ax2.scatter(lhit.x, lhit.y, c="red", s=25, marker="x")
        for rhit in unmatched[1]:
            ax5.scatter(rhit.x, rhit.y, c="red", s=25, marker="x")
        # for i, peak in enumerate(leftpeaks):
        #     if i in lpair:
        #         continue
        #     ax2.scatter(peak[1], peak[0], c="red", s=25, marker="x")
        # for i, peak in enumerate(rightpeaks):
        #     if i in rpair:
        #         continue
        #     ax5.scatter(peak[1], peak[0], c="red", s=25, marker="x")
        plt.show()
        for i, (xpos, ypos) in enumerate(hits):
            plt.annotate(str(i), xy=(xpos+2.5, ypos+2.5), c="black")
        plt.errorbar([hit.x for hit in hits], [hit.y for hit in hits], xerr=[hiterr.x for hiterr in hiterrs],
                     yerr=[hiterr.y for hiterr in hiterrs], marker="o", markersize=2.5, linestyle="None", capsize=2.5)
        for i in strip_positions:
            plt.axhline(i, c="purple", alpha=0.2)
        plt.ylim(0, 200)
        plt.xlim(-150, 150)
        plt.xlabel("Horizontal position (mm)")
        plt.ylabel("Vertical position (mm)")
        plt.show()
        # breakpoint()
