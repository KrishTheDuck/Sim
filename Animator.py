import astropy.units as u
import numpy as np
from astropy.time import Time
from hapsira.ephem import Ephem
from hapsira.util import time_range
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from Plotter import Plotter

p = Plotter(body="2024 PT5", location="500@399", start_epoch=Time.now(),
            stop_epoch=Time.now() + 2 * u.year, delta=2 * u.day)

total, asteroid, _ = p.animate_timed(r1=200 * 1000, plot=True)

positions = []
for i, orbit in enumerate(total[:-1]):
    # get time range
    tot_times = (total[i + 1].epoch - total[i].epoch).to_value(u.s)
    positions.append(Ephem
                     .from_orbit(orbit, epochs=time_range(start=total[i].epoch, end=total[i + 1].epoch))
                     .sample())

positions.append(Ephem.from_orbit(total[-1],
                                  epochs=time_range(
                                      start=total[-1].epoch,
                                      end=total[-1].epoch + total[-1].period)).sample())

# Extract coordinates from all trajectories
trajectories = positions
moon_traj = asteroid.sample()

print(trajectories)
print(moon_traj)

# Extract coordinates from all trajectories for satellite
coords = []
for traj in trajectories:
    x = traj.x.to_value(u.km)
    y = traj.y.to_value(u.km)
    z = traj.z.to_value(u.km)
    coords.append(np.column_stack((x, y, z)))

# Moon coordinates from your data
moon_coords = np.column_stack((moon_traj.x.to_value(u.km),
                               moon_traj.y.to_value(u.km),
                               moon_traj.z.to_value(u.km)))

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_aspect('equal')
ax.grid(True)
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')

# Plot Earth
earth = plt.Circle((0, 0), 6371, color='blue', alpha=0.3, label='Earth')
ax.add_artist(earth)

# Plot all trajectories at the start
colors = ['cyan', 'orange', 'green']
for i, coord in enumerate(coords):
    ax.plot(coord[:, 0], coord[:, 1], '--', color=colors[i], alpha=0.5, label=f'Trajectory {i + 1}')

# Plot Moon's trajectory
ax.plot(moon_coords[:, 0], moon_coords[:, 1], '--', color='red', alpha=0.7, label='Moon trajectory')

# Current position of satellite and Moon
satellite, = ax.plot([], [], 'ro', markersize=8, label='Satellite')
moon, = ax.plot([], [], 'ko', markersize=12, label='Moon')


def init():
    satellite.set_data([], [])
    moon.set_data([], [])
    return [satellite, moon]


def animate(frame):
    # Satellite position
    total_points = sum(len(c) for c in coords)
    points_per_traj = [len(c) for c in coords]
    current_frame = frame % total_points

    current_traj = 0
    local_frame = current_frame
    while local_frame >= points_per_traj[current_traj]:
        local_frame -= points_per_traj[current_traj]
        current_traj += 1

    satellite.set_data(coords[current_traj][local_frame, 0],
                       coords[current_traj][local_frame, 1])

    # Moon position from your data
    moon_frame = frame % len(moon_coords)
    moon.set_data(moon_coords[moon_frame, 0], moon_coords[moon_frame, 1])

    return [satellite, moon]


# Create animation
total_frames = sum(len(c) for c in coords)
anim = FuncAnimation(fig, animate, init_func=init,
                     frames=total_frames,
                     interval=50,
                     blit=True)

plt.tight_layout()
plt.legend()
plt.title('Orbital Transfer Animation')

# To save the animation
anim.save('orbital_transfer.gif', writer='pillow', fps=30)

plt.show()
