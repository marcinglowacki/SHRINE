#!/usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
#from tqdm import tqdm

"""
Originally written by Timothy Perrett
Modified by Danica Scott to get profiles via dynamic spectra and allow
variable crop duration
"""

def _main():

    args = get_args()

    X = np.load(f"{args.label}_X_t_{args.DM}.npy", mmap_mode='r')
    Y = np.load(f"{args.label}_Y_t_{args.DM}.npy", mmap_mode='r')
    I = np.load(args.I, mmap_mode='r')

    # crop the data
    X, Y = do_crop(X, Y, I, args.crop_dur, args.force_peak)

    f0 = args.f0                      # centre frequency (MHz)
    bw = args.bw                      # bandwidth (MHz)
    dDM = args.dDM                    # DM step
    cDM = args.cDM                    # DM count
    dt = args.dt                      # time resolution to return in us
    if cDM == 0:
        DMs = np.arange(args.DM_low, args.DM_high, dDM)  # DM range
    else:
        DMs = np.linspace(args.DM_low, args.DM_high, cDM)
    
    I_profs = []

    nchan = X.shape[0]
    freqs = get_freqs(f0, bw, nchan)

    for DM in DMs:
        I_profs.append(do_DM(fft(X), fft(Y), DM, freqs, dt, bw))

    I_profs = np.array(I_profs)
    #I_profs = np.transpose(np.transpose(I_profs)[3000:8000])

    np.save(f"{args.label}_DMs.npy", DMs)
    np.save(f"{args.label}_I_{dt}us.npy", I_profs)

    np.savetxt(f"{args.label}_DMs.dat", DMs)
    np.savetxt(f"{args.label}_I_{dt}us.dat", I_profs)

    summary_file = open(f"{args.label}_profile_summaryfile.txt", "w")
    summary_file.write(f"//begin generate_profiles summary//\n/*\n")
    summary_file.write(f"X and Y files were of length {len(X)}\n")
    summary_file.write(f"DM range was from {DMs[0]} to {DMs[-1]}\n")
    summary_file.write(f"DM range was made up of {len(DMs)} steps of size {DMs[1]-DMs[0]}\n")
    summary_file.write(f"I profile was of shape {I_profs.shape}\n")
    summary_file.write(f"*/\n//end generate_profiles summary//\n\n\n")
    summary_file.close()

def get_args() -> ArgumentParser:
    """
    Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description = "Generates an I(DM,t) profile for a given FRB.", 
        formatter_class = ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-l", "--label",
        type = str,
        help = "FRB label"
    )
    parser.add_argument(
        "-d", "--DM",
        type = str,
        help = "Dispersion measure to which the input has been dedispersed in pc/cm3"
    )
    parser.add_argument(
        "-r", "--DM_range",
        type = float,
        help = "Range of delta DMs either side of 0 to profile in pc/cm3",
        default = 1
    )
    parser.add_argument(
        "-L", "--DM_low",
        type = float,
        help = "Lower bound of delta DM range",
        default = -5
    )
    parser.add_argument(
        "-H", "--DM_high",
        type = float,
        help = "Upper bound of delta DM range",
        default = 5
    )
    parser.add_argument(
        "--dDM",
        type = float,
        help = "Delta DM step size in pc/cm3",
        default = 0.1
    )
    parser.add_argument(
        "--cDM",
        type = int,
        help = "Number of DM steps to divide the range into",
        default = 0
    )
    parser.add_argument(
        "-t", "--dt",
        type = int,
        help = "Time resolution to return in us",
        default = 1
    )
    parser.add_argument(
        "-f", "--f0", 
        type = float,
        help = "Central frequency of observation in MHz"
    )
    parser.add_argument(
        "-b", "--bw",
        type = int,
        help = "Bandwidth of observation in MHz",
        default = 336
    )
    parser.add_argument(
        "--crop_dur",
        type = int,
        help = "Duration of the crop in ms",
        default=10
    )
    parser.add_argument(
        "-I",
        type = str,
        help = "I profile filename (for doing the crop)"
    )
    parser.add_argument(
        "--force_peak",
        type = int,
        help = "Force the peak to be at this time step (in ms)",
        default=None
    )
    return parser.parse_args()


def do_crop(X, Y, I, crop_dur_ms, force_peak=None):
    if force_peak is None:
        # find rough peak
        I_red = tscrunch(I[:,:-1].sum(axis=0), 1000)

        peak_1ms = np.argmax(I_red) # to nearest ms
    else:
        peak_1ms = force_peak

    peak_3ns = peak_1ms * 1000*336

    crop_dur_3ns = crop_dur_ms * 1000*336
    crop = slice(peak_3ns - crop_dur_3ns//2, peak_3ns + crop_dur_3ns//2)

    return X[crop].copy(), Y[crop].copy()


def do_DM(X_f, Y_f, DM, freqs, dt, bw):
    X_dd = ifft(dedisperse(X_f, DM, freqs))
    Y_dd = ifft(dedisperse(Y_f, DM, freqs))

    # modification: generate profiles via dynamic spectra so we can
    # normalise the channels and get rid of the sometimes-bad top channel
    X_ds = generate_dynspec(X_dd)
    Y_ds = generate_dynspec(Y_dd)

    I_ds = np.abs(X_ds)**2 + np.abs(Y_ds)**2

    # normalise channels
    for i in range(I_ds.shape[1]):
        I_ds[:,i] -= np.mean(I_ds[:,i])
        I_ds[:,i] /= np.std(I_ds[:,i])

    I_dd = I_ds[:,:-1].sum(axis=1)

    I_red = tscrunch(I_dd, dt)

    return I_red


def dedisperse(
    spec: np.ndarray, DM: float, freqs: np.ndarray
) -> np.ndarray:
    """
    Coherently dedisperse the given complex spectrum.

    Coherent dedispersion is performed by applying the inverse of the
    transfer function that acts on radiation as it travels through a
    charged medium. This is detailed in Lorimer & Kramer's Handbook of
    Pulsar Astronomy (2005, Cambridge University Press).

    In practice, this is a frequency-dependent rotation of the complex
    spectrum. None of the amplitudes are altered.

    :param spec: complex 1d-spectrum in a single polarisation
    :type spec: :class:`np.ndarray`
    :param DM: Dispersion measure to dedisperse to (pc/cm3)
    :type DM: float
    :param freqs: array of frequencies
    :type freqs: :class:`np.ndarray`
    :return: Coherently dedispersed complex spectrum
    :rtype: :class:`np.ndarray`
    """
    k_DM = 2.41e-4

    #f_ref = np.min(freqs)
    f_ref = np.median(freqs)

    dedisp_phases = np.exp(
        2j * np.pi * DM / k_DM * ((freqs - f_ref) ** 2 / f_ref ** 2 / freqs * 1e6)
    )

    spec *= dedisp_phases

    return spec


def get_freqs(f0: float, bw: float, nchan: int) -> np.ndarray:
    """
    Create array of frequencies.

    The returned array is the central frequency of `nchan` channels
    centred on `f0` with a bandwidth of `bw`.

    :param f0: Central frequency (arb. units, must be same as `bw`)
    :type f0: float
    :param bw: Bandwidth (arb. units, must be same as `f0`)
    :type bw: float
    :param nchan: Number of channels
    :type nchan: int
    :return: Central frequencies of `nchan` channels centred on `f0`
        over a bandwidth `bw`
    :rtype: :class:`np.ndarray`
    """
    fmin = f0 - bw / 2
    fmax = f0 + bw / 2

    chan_width = bw / nchan

    freqs = np.linspace(fmax, fmin, nchan, endpoint=False) + chan_width / 2

    return freqs


def tscrunch(a: np.ndarray, n: int) -> np.ndarray:
    """
    Reduces the time resolution of array `a` by a factor `n`.

    :param a: The array to scrunch
    :type a: :class:`np.ndarray`
    :param n: The factor by which to reduce `a`
    :type n: int
    :return: Array with reduced time resolution
    :rtype: :class:`np.ndarray`
    """

    a_red = []

    for i in range(int(a.shape[0] / n)):
        a_red.append(np.sum(a[i * n:(i + 1) * n], axis=0))

    a_red = np.array(a_red) / np.sqrt(n)

    return a_red


def generate_dynspec(
        t_ser: np.ndarray, nchan: int = 336,
    ) -> np.ndarray:
    """
    Creates a dynamic spectrum at the highest time resolution from the 
    given time series.

    :param t_ser: input time series of voltages
    :param nchan: number of frequency channels [Default = 336]
    :return: dynamic spectrum of voltages
    :rtype: :class:`np.ndarray`
    """
    dynspec = np.zeros(
        (int(t_ser.shape[0] / nchan), nchan), dtype=np.complex64
    )

    for i in range(int(t_ser.shape[0] / nchan)):
        dynspec[i, :] = np.fft.fft(t_ser[i * nchan : (i + 1) * nchan])

    return dynspec


if __name__ == "__main__":
    _main()
