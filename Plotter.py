import astropy.time
import astropy.units as u
import hapsira.frames
import numpy as np
from astropy.time import Time, TimeDelta
from hapsira.bodies import Earth
from hapsira.core.maneuver import (
    bielliptic,
)
from hapsira.maneuver import Maneuver
from hapsira.plotting import OrbitPlotter
from hapsira.twobody import Orbit
from hapsira.util import time_range, norm
from matplotlib import pyplot as plt

from JPLQuery import JPLQuery
from Mechanics import Mechanics


class Plotter:
    rand_colors = ['red', 'green', 'blue', 'yellow', 'orange', 'violet']

    def __init__(self):
        pass

    @staticmethod
    def generate_orbit(ephemerides: dict, EPOCH: astropy.time.Time, r_min: float, r_max: float,
                       plane: hapsira.frames.Planes, true_anomaly: float):
        a = (r_min + r_max) / 2.0
        ecc = (r_max - r_min) / (r_max + r_min)

        return Orbit.from_classical(
            attractor=Earth,
            plane=plane,
            epoch=EPOCH,
            a=a * u.km,
            ecc=ecc * u.one,
            inc=ephemerides['inc'] * (np.pi / 180.0) * u.rad,
            raan=ephemerides['raan'] * (np.pi / 180.0) * u.rad,
            argp=ephemerides['argp'] * (np.pi / 180.0) * u.rad,
            nu=true_anomaly * u.rad
        )

    @staticmethod
    def plot_timed(body, start_time, end_time, delta, r1=100 * 1000.0 + Earth.R.value / 1000.0, show=False):
        # Get the JPLQuery object located at the earth, starting from the initial query times given.
        s = JPLQuery(body, "500@399", {'start': start_time, 'stop': end_time, 'step': '1d'})
        # Zoom in to the closest approach and let the epochs
        vals = s.zoom_to_closest_approach()

        hor = s.get_ephem(attractor=Earth,
                          epoch=time_range(start=vals['date'] - delta,
                                           end=vals['date'] + delta))

        ephemerides = s.get_orbital_elements(
            epochs_time_range=time_range(start=vals['date'] - delta,
                                         end=vals['date'] + delta),
            EPOCH=vals['date'])

        # Precalculate the initial orbit
        r_max = (vals['dist'] * u.km).to(u.m).value
        r_min = (r1 * u.km).to(u.m).value
        a = (r_min + r_max) / 2.0
        e = (r_max - r_min) / (r_min + r_max)

        dv = abs(Mechanics.vel(r_min, r_min) - Mechanics.vel(a, r_min))
        dt = TimeDelta(Mechanics.elliptical_period(a, e, r_min, r_max) * u.s)

        # first transfer from MEO to asteroid
        initial_orbit = Plotter.generate_orbit(ephemerides, vals['date'] - dt, r1, r1,
                                               hor.plane, 0)
        r, v = initial_orbit.rv()
        v = v.to(u.m / u.s)
        r = r.to(u.m)

        rendezvous_orbit = initial_orbit.apply_maneuver(
            Maneuver((0 * u.s, ((dv * v / norm(v)).value * (u.m / u.s))))
        ).propagate(dt)

        # Calculate dynamics of the return orbit
        dv = abs(Mechanics.vel(r_min, r_min) - Mechanics.vel(a, r_max))

        return_orbit = rendezvous_orbit.apply_maneuver(
            Maneuver((0 * u.s, ((dv * r / norm(r)).value * (u.m / u.s))))
        )

        return_time = return_orbit.time_to_anomaly(return_orbit.nu + np.pi * u.rad)
        return_orbit = return_orbit.propagate(return_time)

        # calculate adjustment orbit by computing a bi-elliptic transfer to a lower orbit.
        adjustment_orbit = []

        rv = return_orbit.rv()
        rv = (rv[0].to_value(u.m), rv[-1].to_value(u.m / u.s))

        # first orbit arbitrarily can be 20 percent bigger than initial
        dv_a, dv_b, dv_c, t_trans1, t_trans2 = bielliptic(Earth.k.value, r_min * 1.2, r_min, rv)

        man1 = Maneuver(
            (0 * u.s, (dv_a * u.m / u.s).decompose())
        )

        adjustment_orbit.append(return_orbit.apply_maneuver(man1).propagate(t_trans1 * u.s + t_trans2 * u.s))

        # second little transfer
        adjustment_orbit.append(adjustment_orbit[-1].apply_maneuver(
            Maneuver.hohmann(adjustment_orbit[-1], r_min * u.m)
        ))

        if show:
            frame = OrbitPlotter()
            frame.set_body_frame(Earth)

            frame.plot(initial_orbit, label="Initial Orbit", color="deepskyblue")
            frame.plot(rendezvous_orbit, label="Rendezvous Orbit", color="darkorange")
            frame.plot(return_orbit, label="Return Orbit", color="mediumseagreen")
            frame.plot(adjustment_orbit[0], label="Adjustment Orbit", color="gold")
            frame.plot(adjustment_orbit[1], label="Adjustment Orbit", color="gold")

            frame.plot_ephem(hor,
                             epoch=vals['date'], label=body,
                             color='red')

            plt.tight_layout()
            plt.show()


if __name__ == "__main__":
    Plotter.plot_timed("2006 WB",
                       Time.now().strftime("%Y-%m-%d 00:00:00"),
                       (Time.now() + TimeDelta(1 * u.year)).strftime("%Y-%m-%d 00:00:00"),
                       0.5 * u.day,
                       show=True)
