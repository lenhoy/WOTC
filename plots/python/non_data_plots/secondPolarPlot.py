import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Set up the plot for a technical diagram
fig, ax = plt.subplots(figsize=(8, 7))

# Configure font to match standard academic/LaTeX style (Serif)
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "mathtext.fontset": "cm",  # Computer Modern (TeX-like)
})

# Define key coordinates
origin = np.array([0, 0])
# Length rho and angle psi (bearing from North)
rho = 6
psi_deg = 30  # Angle from North
psi_rad = np.radians(psi_deg)

# Calculate point p0 based on bearing psi
# In standard math, 0 is East. Here 0 is North.
# x (East) = rho * sin(psi)
# y (North) = rho * cos(psi)
p0_x = rho * np.sin(psi_rad)
p0_y = rho * np.cos(psi_rad)
p0 = np.array([p0_x, p0_y])

# Axis limits
ax_lim = 8
ax.set_xlim(-2, ax_lim)
ax.set_ylim(-1, ax_lim)
ax.set_aspect('equal')
ax.axis('off')  # Turn off the box frame

# --- 1. Draw Axes (North/East) ---
# North Axis
ax.arrow(0, 0, 0, ax_lim - 0.5, head_width=0.3, head_length=0.4, fc='black', ec='black', lw=1.5)
ax.text(-0.8, ax_lim - 0.5, "North", fontsize=18, ha='center')

# East Axis
ax.arrow(0, 0, ax_lim - 0.5, 0, head_width=0.3, head_length=0.4, fc='black', ec='black', lw=1.5)
ax.text(ax_lim - 0.5, -0.8, "East", fontsize=18, va='center')

# Origin Label (o_n) - lowercase o, subscript n
ax.text(-0.4, -0.4, r"$o_n$", fontsize=16)

# --- 2. Draw Object CO (Triangle) ---
# Triangle vertices near origin
tri_scale = 0.6
triangle = patches.Polygon([
    (0, 0), 
    (-0.5*tri_scale, -0.8*tri_scale), 
    (0.5*tri_scale, -0.8*tri_scale)
], closed=True, fill=False, edgecolor='black', lw=1.5)
ax.add_patch(triangle)
# Label CO (upright, no sub/super)
ax.text(0, -1.2*tri_scale - 0.2, "CO", fontsize=14, ha='center')

# --- 3. Draw Vector rho to p0 ---
ax.plot([0, p0_x], [0, p0_y], 'k-', lw=1.5)
# Label rho (midpoint)
ax.text(p0_x/2 + 0.3, p0_y/2, r"$\rho$", fontsize=16)

# Point p0
ax.plot(p0_x, p0_y, 'ko', markersize=5)
# Label p0 (lowercase p, subscript 0)
ax.text(p0_x - 0.5, p0_y + 0.2, r"$p_0$", fontsize=16)


# --- 4. Angles ---
# Psi (Heading from North)
# Arc from North (90 deg in std polar) to vector
arc_rad = 1.5
theta1 = 90 - psi_deg
theta2 = 90
arc_psi = patches.Arc((0, 0), arc_rad*2, arc_rad*2, theta1=theta1, theta2=theta2, color='black')
ax.add_patch(arc_psi)
# Psi Label
ax.text(0.5, 2.0, r"$\psi$", fontsize=16)

# Gamma (Angle at p0)
# This represents a local vertical reference usually
# Let's draw a vertical line at p0 for reference
ax.plot([p0_x, p0_x], [p0_y, p0_y + 2], 'k-', lw=0.5, alpha=0.6)
# Arc for gamma
arc_gamma = patches.Arc((p0_x, p0_y), 1.5, 1.5, theta1=0, theta2=150, color='black') # Arbitrary visual angle
ax.add_patch(arc_gamma)
ax.text(p0_x + 0.5, p0_y + 0.5, r"$\gamma$", fontsize=16)


# --- 5. Projections and Dimensions ---

# -- North Axis Projections --
# Tick marks
tick_len = 0.2
# x^n (at origin/start level) - technically 0 but diagram shows offset
ax.plot([-tick_len, tick_len], [0, 0], 'k-', lw=1.5) # x^n
ax.text(-0.6, 0, r"$x^n$", fontsize=16, va='center')

# x_0^n (at p0 level)
ax.plot([-tick_len, tick_len], [p0_y, p0_y], 'k-', lw=1.5) # x_0^n
ax.text(-0.6, p0_y, r"$x_0^n$", fontsize=16, va='center')

# Dimension line (North)
dim_x = -1.5
ax.annotate('', xy=(dim_x, 0), xytext=(dim_x, p0_y), arrowprops=dict(arrowstyle='<->', lw=1.2))
ax.text(dim_x - 0.3, p0_y/2, r"$\rho \cos(\gamma)$", fontsize=14, rotation=90, va='center', ha='right')


# -- East Axis Projections --
# y^n (at origin)
ax.plot([0, 0], [-tick_len, tick_len], 'k-', lw=1.5)
ax.text(0, -0.6, r"$y^n$", fontsize=16, ha='center')

# y_0^n (at p0 x-coord)
ax.plot([p0_x, p0_x], [-tick_len, tick_len], 'k-', lw=1.5)
ax.text(p0_x, -0.6, r"$y_0^n$", fontsize=16, ha='center')

# Dimension line (East)
dim_y = -1.5
ax.annotate('', xy=(0, dim_y), xytext=(p0_x, dim_y), arrowprops=dict(arrowstyle='<->', lw=1.2))
ax.text(p0_x/2, dim_y - 0.5, r"$\rho \sin(\gamma)$", fontsize=14, ha='center')

plt.tight_layout()
plt.show()