import ROOT as root


def roothist(data, bin_width, bin_min=None, bin_max=None, name="test1234"):
    if not bin_min:
        min_val = min(data)
    else:
        min_val = bin_min
    if not bin_max:
        max_val = max(data)
    else:
        max_val = bin_max
    nbin = int((max_val-min_val) / bin_width)
    th1d = root.TH1D(name, name, nbin, min_val, max_val)
    for value in data:
        th1d.Fill(value)
    return th1d
