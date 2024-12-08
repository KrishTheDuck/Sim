from copy import copy

import astropy.time
import astropy.units as u
import hapsira.frames
import numpy as np
from astropy.time import TimeDelta, Time
from hapsira.bodies import Earth
from hapsira.maneuver import Maneuver
from hapsira.plotting import OrbitPlotter
from hapsira.twobody import Orbit
from hapsira.util import time_range, norm
from matplotlib import pyplot as plt

from JPLQuery import JPLQuery
from Mechanics import Mechanics


class Plotter:
    __colors = ["deepskyblue", "darkorange", "mediumseagreen", "gold", "blue", "green", "darkred", "cyan"]
    __ephemerides_search_range = 1 * u.day

    def __init__(self, body, location, start_epoch, stop_epoch):
        self.body = body
        self.location = location

        self.s = JPLQuery(self.body, self.location, {
            'start': start_epoch.strftime("%Y-%m-%d 00:00:00"),
            'stop': stop_epoch.strftime("%Y-%m-%d 00:00:00"),
            'step': '1d'
        })

        self.vals = self.s.zoom_to_closest_approach()

        self.hor = self.s.get_ephem(attractor=Earth,
                                    epoch=time_range(start=self.vals['date'] - Plotter.__ephemerides_search_range,
                                                     end=self.vals['date'] + Plotter.__ephemerides_search_range))

        self.ephemerides = self.s.get_orbital_elements(
            epochs_time_range=time_range(start=self.vals['date'] - Plotter.__ephemerides_search_range,
                                         end=self.vals['date'] + Plotter.__ephemerides_search_range),
            EPOCH=self.vals['date'])

    @staticmethod
    def generate_orbit(ephemerides: dict, EPOCH: astropy.time.Time, r_min: float, r_max: float,
                       plane: hapsira.frames.Planes, true_anomaly: float):
        a = (r_min + r_max) / 2.0
        ecc = (r_max - r_min) / (r_max + r_min)

        return Orbit.from_classical(
            attractor=Earth,
            plane=plane,
            epoch=EPOCH,
            a=a * u.m,
            ecc=ecc * u.one,
            inc=ephemerides['inc'] * (np.pi / 180.0) * u.rad,
            raan=ephemerides['raan'] * (np.pi / 180.0) * u.rad,
            argp=ephemerides['argp'] * (np.pi / 180.0) * u.rad,
            nu=true_anomaly * u.rad
        )

    def animate_timed(self, r_min=100 * 1000.0 * 1000.0 + Earth.R.to_value(u.m)):
        maneuvers = []
        r_max = (self.vals['dist'] * u.km).to_value(u.m)

        # Region Initial Orbit Preconditions
        a = (r_min + r_max) / 2.0
        e = (r_max - r_min) / (r_min + r_max)

        dv = abs(Mechanics.vel(r_min, r_min) - Mechanics.vel(a, r_min))
        dt = TimeDelta(Mechanics.elliptical_period(a, e, r_min, r_max) * u.s)
        # End
        # Region: first transfer from MEO to asteroid
        orbit = (Plotter.generate_orbit(self.ephemerides, self.vals['date'] - dt,
                                        r_min, r_min, self.hor.plane, 0))
        orb0 = copy(orbit)

        v = orbit.v.to(u.m / u.s)
        r = orbit.r.to(u.m)

        man = Maneuver((0 * u.s, Mechanics.vec_in_dir(dv, v) * (u.m / u.s)))
        maneuvers.append(man)
        orbit = orbit.apply_maneuver(man)
        dt = abs(orbit.time_to_anomaly(np.pi * u.rad))

        orbit = orbit.propagate_to_anomaly(np.pi * u.rad)
        # End
        # Region: Calculate dynamics of the return orbit
        dv = Mechanics.vel(r_min, r_min) - Mechanics.vel(a, r_max)
        dv_vector = Mechanics.vec_in_dir(dv, r) * (u.m / u.s)

        man = Maneuver((0 * u.s, dv_vector))
        maneuvers.append(Maneuver((dt.to(u.s), dv_vector)))

        orbit = orbit.apply_maneuver(man)
        return_time = orbit.time_to_anomaly(np.pi * u.rad + orbit.nu)
        orbit = orbit.propagate(return_time)

        # End
        # Region: adjustment
        # Velocity difference vector
        dv_vector = orbit.v.to(u.m / u.s) - v
        dv = norm(dv_vector).value

        # Scale thrust and construct maneuver
        maneuvers.append(Maneuver((return_time,
                                   (u.m / u.s) * Mechanics.vec_in_dir(-dv, dv_vector))))

        print(maneuvers)
        # End
        # Coalesce the maneuvers
        total_maneuvers = Maneuver(*[(m[0], m[1]) for man in maneuvers for m in man])

        complete = [orb0]+ orb0.apply_maneuver(total_maneuvers, intermediate=True)

        return complete, self.s.get_ephem(attractor=Earth,
                                          epoch=time_range(
                                              start=self.vals['date'] - dt,
                                              end=self.vals['date'] + dt)), total_maneuvers

    def orbits(self, r_min=200 * 1000 * 1000):
        r_max = (self.vals['dist'] * u.km).to_value(u.m)
        a = (r_min + r_max) / 2.0
        ecc = (r_max - r_min) / (r_min + r_max)
        orbits = []

        # Region: initial orbit
        dt = TimeDelta(Mechanics.elliptical_period(a, ecc, r_min, r_max) * u.s)

        orbits.append(Plotter.generate_orbit(self.ephemerides, self.vals['date'] - dt,
                                             r_min, r_min, self.hor.plane, 0))
        r = orbits[0].r.to(u.m)
        v = orbits[0].v.to(u.m / u.s)
        # End
        # Region: Second Orbit
        orbits.append(Plotter.generate_orbit(self.ephemerides, self.vals['date'] - dt,
                                             r_min, r_max, self.hor.plane, 0))
        # End
        # Region: Return orbit
        dv = abs(Mechanics.vel(r_min, r_min) - Mechanics.vel(a, r_max))

        v_ellipse = Mechanics.vec_in_dir(Mechanics.vel(a, r_max), -v) * (u.m / u.s)
        v_radial = Mechanics.vec_in_dir(dv, r) * (u.m / u.s)
        position = Mechanics.vec_in_dir(r_max, -r) * u.m

        orbits.append((orb := Orbit.from_vectors(r=position, v=(v_ellipse + v_radial).to(u.km / u.s),
                                                 epoch=self.vals['date'], plane=self.hor.plane,
                                                 attractor=Earth))
                      .propagate_to_anomaly(orb.nu + np.pi * u.rad))

        # End
        # Region: adj
        orbits.append(Plotter.generate_orbit(self.ephemerides, orbits[-1].epoch,
                                             r_min, r_min, self.hor.plane, 0)
                      .propagate(orbits[0].period / 2))
        # End
        # ! plot
        return orbits

    def plot(self, orbits, title, titles):
        frame = OrbitPlotter()
        frame.set_body_frame(Earth)

        for index, o in enumerate(orbits):
            frame.plot(o, label=f"{titles[index]}",
                       color=Plotter.__colors[index % len(Plotter.__colors)])

        frame.plot_ephem(self.hor, epoch=self.vals['date'],
                         label=f'{self.body}', color='red')
        plt.title(title)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    bodies = ['2021 WA5', '2022 YO1', '2007 XB23', '2024 UU3', '2024 PT5', '2020 XR']

    for i, body in enumerate(bodies):
        p = Plotter(body=body, location="500@399", start_epoch=Time.now(),
                    stop_epoch=Time.now() + 2 * u.year)

        ol, _, mans = p.animate_timed(r_min=200 * 1000 * 1000)
        print(len(ol))
        time = mans.get_total_time().to_value(u.day)
        cost = mans.get_total_cost().to_value(u.km / u.s)

        p.plot(ol, f"{body}, cost: {cost:.3} km/s, time: {time:.3} days",
               ['Initial Orbit', 'Rendezvous Orbit', 'Return Orbit After Shooting Ion Plume',
                'Corrective Orbit Back to Initial'])
