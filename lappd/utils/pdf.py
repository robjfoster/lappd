import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from lappd.utils import lappdcfg as cfg


class PDF:

    def __init__(self, var, start, end, npoints) -> None:
        self.var = var
        self.start = start
        self.end = end
        self.npoints = npoints
        self.xs = np.linspace(self.start, self.end, self.npoints)
        self.pdf = self.func(self.xs)
        self.pdfsum = np.sum(self.pdf)
        self.pdf = self.pdf / self.pdfsum  # Sum of PDF should be 1

    def func(self, value):
        raise NotImplementedError(
            "Method must be implemented from a derived class")

    def sample(self, value):
        return self.func(value)

    def plot(self):
        plt.plot(self.xs, self.pdf)
        plt.show()


class TimePDF(PDF):

    def __init__(self, var, start, end, step) -> None:
        super().__init__(var, start, end, step)

    def func(self, value):
        return (1 - (0.5 * (1 + erf((value - (cfg.STRIPLENGTH / (cfg.STRIPVELOCITY * cfg.SOL))) / (self.var * np.sqrt(2))))))


class LocPDF(PDF):

    def __init__(self, var, start, end, npoints) -> None:
        super().__init__(var, start, end, npoints)

    def func(self, value):
        return ((1 / (self.var * np.sqrt(2 * np.pi))) * np.exp(-(value)**2 / (2 * self.var**2)))


class AmplPDF(PDF):

    def __init__(self, var, start, end, npoints) -> None:
        super().__init__(var, start, end, npoints)

    def func(self, value):
        return ((1 / (self.var * np.sqrt(2 * np.pi))) * np.exp(-(value)**2 / (2 * self.var**2)))


class ChargeSharePDF(PDF):

    def __init__(self, var, start, end, npoints) -> None:
        super().__init__(var, start, end, npoints)

    def func(self, value):
        return ((1 / (self.var * np.sqrt(2 * np.pi))) * np.exp(-(value)**2 / (2 * self.var**2)))


if __name__ == "__main__":
    timepdf = TimePDF(0.1, 0, 2, 1000)
    locpdf = LocPDF(0.3, 0, 20, 1000)
    amplpdf = AmplPDF(5, 0, 100, 1000)
    timepdf.plot()
    locpdf.plot()
    amplpdf.plot()
    breakpoint()
