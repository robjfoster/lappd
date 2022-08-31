from collections import namedtuple
from math import sqrt
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from skimage.feature import peak_local_max

import lappd.utils.sigutils as su
from lappd.utils import lappdcfg as cfg
from lappd.utils.interpolation import gaussMM, image_gaussian
from lappd.utils.lognormal import lngaus2D
from lappd.utils.pdf import AmplPDF, LocPDF, TimePDF

Hit = namedtuple("Hit", "x y")
Hit3D = namedtuple("Hit3D", "x y z")
RecoHit = namedtuple("RecoHit", "recox recoy recot x y z")
Pair = namedtuple("Pair", "left right")

x = np.arange(0, cfg.NSAMPLES - cfg.NREMOVEDSAMPLES, 1)
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
    return 100 - ((cfg.STRIPWIDTH * (strip+1)) + (cfg.STRIPSPACING/2.0))


def z_to_ampl(z):
    return z


def transverse_position(lefttime, righttime):
    timedelta = lefttime - righttime
    pos = timedelta / cfg.MAXDELTA * (cfg.STRIPLENGTH / 2.0)
    #pos = timedelta * (cfg.STRIPVELOCITY*cfg.SOL) / 2.0
    timedeltaerror = sqrt((cfg.TIMEJITTER)**2 + (cfg.TIMEJITTER)**2)
    poserror = pos * sqrt((timedeltaerror / timedelta) **
                          2 + (cfg.MAXDELTAERROR / cfg.MAXDELTA)**2)
    return pos, abs(poserror)


def cfd_timestamp(matrix: np.ndarray, hit: Hit3D, times=None):
    strippulse = matrix[hit.y]
    timestamp, status = su.cfd(
        strippulse*-1, 0.2, userpeak=hit.x, times=times, plot=True)
    return timestamp, status


def get_centroid(matrix: np.ndarray, hit: Hit3D, width: int = 10) -> float:
    if hit.y+width >= matrix.shape[0]:
        upper = matrix.shape[0]-1 - hit.y
    else:
        upper = width
    if hit.y-width < 0:
        lower = -hit.y
    else:
        lower = width
    slicedmatrix = matrix[hit.y-lower:hit.y+upper, hit.x]
    # if we assume the hit profile across strips is gaussian...
    centroid = np.mean(slicedmatrix)
    return hit.y - upper + centroid


def get_centroid_combined(
    leftmatrix: np.ndarray,
    rightmatrix: np.ndarray,
    lefthit: Hit3D,
    righthit: Hit3D,
    width: int = 10
) -> float:
    if lefthit.y+10 >= leftmatrix.shape[0]:
        leftupper = leftmatrix.shape[0]-1
    else:
        leftupper = width
    if lefthit.y-10 < 0:
        leftlower = 0
    else:
        leftlower = width
    if righthit.y+10 >= rightmatrix.shape[0]:
        rightupper = rightmatrix.shape[0]-1
    else:
        rightupper = width
    if righthit.y-10 < 0:
        rightlower = 0
    else:
        rightlower = width
    lower = max(leftlower, rightlower)
    upper = min(leftupper, rightupper)
    leftslicedmatrix = leftmatrix[lefthit.y-lower:lefthit.y+upper, lefthit.x]
    rightslicedmatrix = rightmatrix[righthit.y -
                                    lower:righthit.y+upper, righthit.x]
    slicedmatrix = (leftslicedmatrix + rightslicedmatrix) / 2.0
    centroid = np.mean(slicedmatrix)
    return ((lefthit.y + righthit.y) / 2.0) - upper + centroid


def detect_peaks(matrix, threshold=3, min_distance=5, exclude_border=False):
    # neighborhood = generate_binary_structure(2, 2)
    # image_max = maximum_filter(matrix, size=2, mode="constant")
    coordinates = peak_local_max(
        matrix, threshold_abs=threshold, min_distance=min_distance, exclude_border=exclude_border)
    peaks = []
    for coord in coordinates:
        peaks.append(Hit3D(coord[1], coord[0], matrix[coord[0], coord[1]]))
    return peaks


def plot_hitmap(hits, bins=150):
    heatmap, xedges, yedges = np.histogram2d([hit.recox for hit in hits], [
                                             hit.recoy for hit in hits], bins=bins, range=[[-150, 150], [-100, 100]])
    extent = [-150, 150, 0, 200]
    plt.clf()
    plt.imshow(heatmap.T, extent=extent, origin='lower')
    plt.colorbar()
    plt.title("BGO + Sr90 masked 0-160 mm in x")
    plt.xlabel("Reconstructed x position (mm)")
    plt.ylabel("Reconstructed y position (mm)")
    plt.show()


def match_peaks(leftpeaks: np.ndarray,
                rightpeaks: np.ndarray,
                leftcleaninterp,
                rightcleaninterp,
                min_like: float = 1e-15,
                interpfactor: float = 10
                ) -> List[Tuple[int]]:
    if len(leftpeaks) == 0 or len(rightpeaks) == 0:
        return [], [[], []]
    likelihood = []
    for lpeak in leftpeaks:
        # Convert x, y, z to t, loc, ampl
        lpeakt = x_to_t(lpeak.x, interpfactor=interpfactor)
        lpeakloc = y_to_loc(lpeak.y, interpfactor=interpfactor)
        lpeakampl = z_to_ampl(
            leftcleaninterp[lpeak.y, lpeak.x])
        thislike = []
        for rpeak in rightpeaks:
            rpeakt = x_to_t(rpeak.x, interpfactor=interpfactor)
            rpeakloc = y_to_loc(rpeak.y, interpfactor=interpfactor)
            rpeakampl = z_to_ampl(
                rightcleaninterp[rpeak.y, rpeak.x])
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
    min_like = cfg.MINLIKE
    # Copy of likelihood matrix that will be *modified* in the loop.
    # likelihood is left unmodified
    mod_likelihood = np.copy(likelihood)
    pairs = []
    while True:
        # Get index of max value in matrix
        pair = np.unravel_index(
            np.argmax(mod_likelihood), mod_likelihood.shape)
        # Find value of likelihood at this index
        max_like = mod_likelihood[pair[0], pair[1]]
        # print(f'Max likelihood: {max_like}')
        # Compare to our likelihood limit
        if max_like < min_like:
            # print(f"Last likelihood value: {max_like}")
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
    #print(f"Matched {len(pairs)} pairs")
    left_unmatched = list(set(left_unmatched))
    right_unmatched = list(set(right_unmatched))
    # if len(left_unmatched) == 0 and len(right_unmatched) == 0:
    #     print("All pairs matched")
    # else:
    #     print(f"{len(left_unmatched)} left hit unmatched")
    #     print(f"{len(right_unmatched)} right hit unmatched")
    return pairs, (left_unmatched, right_unmatched)


strip_positions = y_to_loc(np.arange(0, 28, 1), interpfactor=1)

if __name__ == "__main__":
    import sys

    from matplotlib import cm

    from lappd.lappdevent import LAPPDEvent
    from lappd.utils import interpolation
    from lappd.utils.wiener import do_wiener
    all_hits = []
    all_times = []
    fancy_plot = False
    unfancy_plot = False
    base_dir = sys.argv[1]
    stripnumber = int(sys.argv[2])
    eventcount = 0
    for levent in LAPPDEvent.itr_all_raw(base_dir, trigger=None):
        if su.threshold(levent.leftmatrix, 70, polarity="positive") \
                or su.threshold(levent.rightmatrix, 70, polarity="positive") \
                or not su.threshold(levent.leftmatrix, 2, polarity="positive") \
                or not su.threshold(levent.rightmatrix, 2, polarity="positive"):
            continue
        # if levent.event_no > 5000:
        #    break
        breakpoint()
        template = np.load("template2d.npy")
        temp1d = np.load("template.npy")
        #print(f"Event number: {levent.event_no}")
        eventcount += 1
        left = levent.leftmatrix
        right = levent.rightmatrix
        thismin = min((np.min(left), np.min(right)))
        thismax = max((np.max(left), np.max(right)))
        leftclean = do_wiener(left, template)
        rightclean = do_wiener(right, template)
        leftcleaninterp = interpolation.interp_matrix(leftclean)
        rightcleaninterp = interpolation.interp_matrix(rightclean)
        leftpeaks = detect_peaks(
            leftcleaninterp, threshold=cfg.MINHEIGHT, exclude_border=(3, 3))
        rightpeaks = detect_peaks(
            rightcleaninterp, threshold=cfg.MINHEIGHT, exclude_border=(3, 3))
        pairs, unmatched = match_peaks(
            leftpeaks, rightpeaks, leftcleaninterp, rightcleaninterp)
        # leftpeakscoarse = detect_peaks(
        #     leftclean, threshold=cfg.MINHEIGHT, exclude_border=(1, 1))
        # rightpeakscoarse = detect_peaks(
        #     rightclean, threshold=cfg.MINHEIGHT, exclude_border=(1, 1))
        # pairscoarse, unmatchedcoarse = match_peaks(
        #     leftpeakscoarse, rightpeakscoarse, leftclean, rightclean, interpfactor=1)
        if len(pairs) == 0 and len(unmatched[0]) == 0 and len(unmatched[1]) == 0:
            continue
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
            ax1.set_xticks(np.arange(0, 1014, 100))
            ax1.set_xticklabels(np.arange(0, 210, 20))
            ax1.set_yticks(np.arange(0, 8, 1))
            ax1.set_yticklabels(np.arange(14, 6, -1))
            ax4.set_xticks(np.arange(0, 1014, 100))
            ax4.set_xticklabels(np.arange(0, 210, 20))
            ax4.set_yticks(np.arange(0, 8, 1))
            ax4.set_yticklabels(np.arange(14, 6, -1))
            ax2.set_xticks(np.arange(0, 10140, 1000))
            ax2.set_xticklabels(np.arange(0, 210, 20))
            ax2.set_yticks(np.arange(0, 80, 10))
            ax2.set_yticklabels(np.arange(14, 6, -1))
            ax5.set_xticks(np.arange(0, 10140, 1000))
            ax5.set_xticklabels(np.arange(0, 210, 20))
            ax5.set_yticks(np.arange(0, 80, 10))
            ax5.set_yticklabels(np.arange(14, 6, -1))
            ax1.set_xlabel("Time (ns)")
            ax1.set_ylabel("Stripnumber")
            ax4.set_xlabel("Time (ns)")
            ax4.set_ylabel("Stripnumber")
            ax2.set_xlabel("Time (ns)")
            ax2.set_ylabel("Stripnumber")
            ax5.set_xlabel("Time (ns)")
            ax5.set_ylabel("Stripnumber")
            fig.colorbar(base_img, ax=[ax0, ax1, ax2, ax3, ax4, ax5],
                         orientation="horizontal", fraction=0.016)
        elif unfancy_plot:
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
                         orientation="horizontal", fraction=0.046, anchor=(1.0, 0.0))
            fig.suptitle(f"Event {levent.event_no}")
            for peak in leftpeaks:
                if peak[0] < 3:
                    continue
                ax[1, 0].scatter(peak[1], peak[0], c="red", s=5, marker="x")
            for peak in rightpeaks:
                if peak[0] < 3:
                    continue
                ax[1, 1].scatter(peak[1], peak[0], c="red", s=5, marker="x")
        hits = []
        times = []
        hiterrs = []
        for i, pair in enumerate(pairs):
            # lpairx = leftpeaks[pair[0]][1]
            # lpairy = leftpeaks[pair[0]][0]
            # rpairx = rightpeaks[pair[1]][1]
            # rpairy = rightpeaks[pair[1]][0]
            leftcfd, leftstatus = cfd_timestamp(leftcleaninterp, pair.left)
            rightcfd, rightstatus = cfd_timestamp(rightcleaninterp, pair.right)
            try:
                trigger_time = levent.trigevent.pulses[0].left.cfpeak
            except:
                trigger_time = 0
            event_time = x_to_t(((leftcfd + rightcfd) / 2.0)) - trigger_time
            if leftstatus is False or rightstatus is False:
                print("Skipped due to CFD failure")
                continue
            xpos, xposerr = transverse_position(
                x_to_t(leftcfd), x_to_t(rightcfd))
            ypos = (y_to_loc(pair.left.y) + y_to_loc(pair.right.y)) / 2.0
            ypos = y_to_loc(get_centroid(leftcleaninterp, pair.left) +
                            get_centroid(rightcleaninterp, pair.right) / 2.0)
            recohit = RecoHit(xpos, ypos, 0, (pair.left.x + pair.right.x) /
                              2.0, (pair.left.y + pair.right.y) / 2.0, (pair.left.z + pair.right.z) / 2.0)
            hits.append(recohit)
            hiterrs.append(Hit(xposerr, 5))
            times.append(event_time)
        if fancy_plot:
            # Plot the matched peaks as pink circles and number them
            for i, pair in enumerate(pairs):
                # xpos, xposerr = transverse_position(
                #     x_to_t(pair.left.x), x_to_t(pair.right.x))
                # ypos = (y_to_loc(pair.left.y) + y_to_loc(pair.right.y)) / 2.0
                # hits.append(Hit(xpos, ypos))
                # hiterrs.append(Hit(xposerr, 5))
                ax2.scatter(pair.left.x, pair.left.y,
                            c="pink", s=25, marker="o")
                ax5.scatter(pair.right.x, pair.right.y,
                            c="pink", s=25, marker="o")
                # , c="green", s=4, marker="x")
                ax2.annotate(str(i), xy=(pair.left.x, pair.left.y), c="white")
                # , c="green", s=4, marker="x")
                ax5.annotate(str(i), xy=(
                    pair.right.x, pair.right.y), c="white")
            lpair = [i[0] for i in pairs]
            rpair = [i[1] for i in pairs]
            # Plot the unmatched peaks as red crosses
            for lhit in unmatched[0]:
                ax2.scatter(lhit.x, lhit.y, c="red", s=25, marker="x")
            for rhit in unmatched[1]:
                ax5.scatter(rhit.x, rhit.y, c="red", s=25, marker="x")
            plt.show()
        for i, (xpos, ypos, _, _, _) in enumerate(hits):
            plt.annotate(str(i), xy=(xpos+2.5, ypos+2.5), c="black")
        if fancy_plot:
            plt.errorbar([hit.recox for hit in hits], [hit.recoy for hit in hits], xerr=[hiterr.x for hiterr in hiterrs],
                         yerr=[hiterr.y for hiterr in hiterrs], marker="o", markersize=2.5, linestyle="None", capsize=2.5)
            for i in strip_positions:
                plt.axhline(i, c="purple", alpha=0.2)
            plt.ylim(0, 200)
            plt.xlim(-150, 150)
            plt.xlabel("Horizontal position (mm)")
            plt.ylabel("Vertical position (mm)")
            plt.show()
        if hits:
            all_hits += hits
        if times:
            all_times += times
    print("Total hits: ", len(all_hits))
    print("Total events analysed: ", eventcount)
    np.save("allhits.npy", all_hits)
    plot_hitmap(all_hits, bins=100)
    breakpoint()
