import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Plot the curve and its projection
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
plt.ion()


# Define the curve y = x^2
def curve_x_y(x):
    return 0.5 * x ** 2


def generate_points_in_triangle(n_points, x1, y1, x2, y2, x3, y3):
    points = []
    for _ in range(n_points):
        # Random barycentric coordinates
        r1 = np.sqrt(np.random.rand())
        r2 = np.random.rand()
        x = (1 - r1) * x1 + r1 * (1 - r2) * x2 + r1 * r2 * x3
        y = (1 - r1) * y1 + r1 * (1 - r2) * y2 + r1 * r2 * y3
        points.append((x, y))
    return np.array(points)


def transform(x, y, th):
    t = np.deg2rad(th)
    return x, y * np.sin(theta) ** 2, y * np.cos(theta) * np.sin(theta)


# Define the inclined plane
theta = -20
deflection_angle = -10

t1 = np.deg2rad(theta)

# Generate points for the curve
x = np.linspace(-1, 1, 100)
y = curve_x_y(x)
z = np.zeros_like(x)

normal_vector = np.array([0, np.cos(t1), np.sin(t1)])
xx, yy, zz = transform(x, y, theta)

# Define extrusion parameters
extrusion_length = 1.0  # Length of extrusion
delta = np.deg2rad(deflection_angle)

# Extrusion vector
v_x = 0
v_y = np.cos(delta)
v_z = np.sin(delta)

# Generate extruded surface points
extruded_x, extruded_y, extruded_z = [], [], []
n_extrusion_steps = 50
for step in np.linspace(0, extrusion_length, n_extrusion_steps):
    extruded_x.append(xx + step * v_x)
    extruded_y.append(yy + step * v_y)
    extruded_z.append(zz + step * v_z)

# Convert to numpy arrays for plotting
extruded_x = np.array(extruded_x)
extruded_y = np.array(extruded_y)
extruded_z = np.array(extruded_z)

# Plot the extruded surface
for i in range(len(extruded_x)):
    ax.plot(extruded_x[i], extruded_y[i], extruded_z[i], color='green', alpha=0.7)

# Show the surface boundary
ax.plot_surface(extruded_x, extruded_y, extruded_z, color='cyan', alpha=0.3)

plt.show()

ax.plot(x, y, z, color='red')
ax.plot(xx, yy, zz, color='blue')

ax.set_xlim([-1, 1])
ax.set_ylim([0, 1])
ax.set_zlim([-1, 1])

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()
input()
