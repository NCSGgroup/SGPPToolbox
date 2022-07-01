"""
    It provides several classes to decrease the correlation of the certain SHCs for the sake of removing the 'stripes'
of the time-variable gravity(TVG) filed.
    Currently, PnMm methods, sliding stable window de-correlation filter and sliding variable window de-correlation filters are available.
"""
import copy

import numpy as np

from pysrc.SHC import SHC


def dec1d(x: np.ndarray, y: np.ndarray, n):
    z = np.polyfit(x, y, n)
    p = np.poly1d(z)
    return y - p(x)


class PnMm:
    def __init__(self, n, m):
        self.poly_n = n
        self.start_m = m

    def ApplyTo(self, *shcs: SHC):
        nmax = shcs[0].nmax

        shcs_filtered = []
        for i in range(len(shcs)):
            this_shc = copy.deepcopy(shcs[i])
            for m in range(self.start_m, nmax + 1):
                length1 = len(shcs[i].Cnm2d[m::2, m])
                length2 = len(shcs[i].Cnm2d[m + 1::2, m])
                if length1 <= self.poly_n or length2 <= self.poly_n:
                    continue

                this_shc.Cnm2d[m::2, m] = dec1d(np.arange(length1), shcs[i].Cnm2d[m::2, m], self.poly_n)
                this_shc.Cnm2d[m + 1::2, m] = dec1d(np.arange(length2), shcs[i].Cnm2d[m + 1::2, m], self.poly_n)

                this_shc.Snm2d[m::2, m] = dec1d(np.arange(length1), shcs[i].Snm2d[m::2, m], self.poly_n)
                this_shc.Snm2d[m + 1::2, m] = dec1d(np.arange(length2), shcs[i].Snm2d[m + 1::2, m], self.poly_n)

            shcs_filtered.append(this_shc)

        if len(shcs) == 1:
            shcs_filtered = shcs_filtered[0]
        return shcs_filtered


class StableWindow:
    def __init__(self, n, m, window_len):
        self.poly_n = n
        self.start_m = m
        self.window_len = window_len

    def ApplyTo(self, *shcs: SHC):
        nmax = shcs[0].nmax

        shcs_filtered = []
        for i in range(len(shcs)):
            this_shc = copy.deepcopy(shcs[i])
            for m in range(self.start_m, nmax - 2 * self.poly_n):
                a = int((self.window_len - 1) / 2)

                for l0 in range(m, nmax + 1):
                    ll = np.array(range(l0 - 2 * a, l0 + 2 * a + 1, 2))
                    if np.min(ll) < m:
                        ll += (m - np.min(ll))
                    if np.max(ll) > nmax:
                        ll -= (np.max(ll) - nmax)

                    yC = np.array([shcs[i].Cnm2d[ll[n]][m] for n in range(len(ll))])
                    zC = np.polyfit(ll, yC, self.poly_n)
                    pC = np.poly1d(zC)
                    this_shc.Cnm2d[l0][m] = shcs[i].Cnm2d[l0][m] - pC(l0)

                    yS = np.array([shcs[i].Snm2d[ll[n]][m] for n in range(len(ll))])
                    zS = np.polyfit(ll, yS, self.poly_n)
                    pS = np.poly1d(zS)
                    this_shc.Snm2d[l0][m] = shcs[i].Snm2d[l0][m] - pS(l0)

            shcs_filtered.append(this_shc)

        if len(shcs) == 1:
            shcs_filtered = shcs_filtered[0]
        return shcs_filtered


class VariableWindow:
    def __init__(self, n, m, window_min, a, k):
        self.poly_n = n
        self.start_m = m
        self.window_a = a
        self.window_k = k
        self.window_min = window_min

    def ApplyTo(self, *shcs: SHC):
        nmax = shcs[0].nmax

        shcs_filtered = []
        for i in range(len(shcs)):
            this_shc = copy.deepcopy(shcs[i])

            for m in range(self.start_m, nmax + 1):
                window_len = max(self.window_a * np.exp(-m / self.window_k) + 1, self.window_min)
                if nmax - m + 1 < 2 * window_len:
                    continue

                a = int((window_len - 1) / 2)
                for l0 in range(m, nmax + 1):
                    ll = np.array(range(l0 - 2 * a, l0 + 2 * a + 1, 2))

                    degree_beyond_flag = np.min(ll) < m or np.max(ll) > nmax
                    while degree_beyond_flag:
                        if np.min(ll) < m:
                            ll += 2
                        if np.max(ll) > nmax:
                            ll -= 2
                        degree_beyond_flag = np.min(ll) < m or np.max(ll) > nmax

                    yC = np.array([shcs[i].Cnm2d[l][m] for l in ll])
                    zC = np.polyfit(ll, yC, self.poly_n)
                    pC = np.poly1d(zC)
                    this_shc.Cnm2d[l0][m] = shcs[i].Cnm2d[l0][m] - pC(l0)

                    yS = np.array([shcs[i].Snm2d[l][m] for l in ll])
                    zS = np.polyfit(ll, yS, self.poly_n)
                    pS = np.poly1d(zS)
                    this_shc.Snm2d[l0][m] = shcs[i].Snm2d[l0][m] - pS(l0)

            shcs_filtered.append(this_shc)

        if len(shcs) == 1:
            shcs_filtered = shcs_filtered[0]
        return shcs_filtered
