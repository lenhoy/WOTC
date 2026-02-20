import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def draw_marine_craft_coord_system():
    """
    Generates a publication-quality plot illustrating the transformation 
    from NED (North-East-Down) cartesian coordinates to Polar states.
    """
    
    # ---------------------------------------------------------
    # 1. Setup and Style
    # ---------------------------------------------------------
    # Use LaTeX font rendering for a professional look (if available)
    try:
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"],
        })
    except:
        # Fallback if LaTeX is not installed
        plt.rcParams.update({
            "font.family": "serif",
        })

    fig, ax = plt.subplots(figsize=(8, 6))
    
    # ---------------------------------------------------------
    # 2. Geometry Definitions
    # ---------------------------------------------------------
    # Target (Origin of the polar system)
    target_pos = np.array([0, 0])
    
    # Vessel State (NED)
    # Placed in the first quadrant (North-East relative to target)
    vessel_n = 5.0  # North position
    vessel_e = 4.0  # East position
    psi_deg = 30    # Heading in degrees
    psi = np.deg2rad(psi_deg)
    
    # Polar State Calculations
    rho = np.sqrt(vessel_n**2 + vessel_e**2)
    gamma = np.arctan2(vessel_e, vessel_n) # Angle from North (Clockwise is positive in NED usually, but for standard polar math we track typically from X. 
                                           # In Marine NED, angles are usually measured clockwise from North.
                                           # Let's visualize strictly as Fossen's definition: gamma is angle of position vector.

    # ---------------------------------------------------------
    # 3. Draw Frames and Axes
    # ---------------------------------------------------------
    
    # Draw NED Global Axes (Origin at Target)
    ax.arrow(0, -1, 0, 9, head_width=0.2, head_length=0.3, fc='k', ec='k', linewidth=1.2) # North Axis
    ax.arrow(-1, 0, 9, 0, head_width=0.2, head_length=0.3, fc='k', ec='k', linewidth=1.2) # East Axis
    
    ax.text(0.2, 8.2, r'$x_n$ (North)', fontsize=12)
    ax.text(8.2, 0.2, r'$y_n$ (East)', fontsize=12)
    
    # Mark the Target
    ax.plot(0, 0, 'ko', markersize=6)
    ax.text(-0.8, -0.6, r'Target $(0,0)$', fontsize=10)

    # ---------------------------------------------------------
    # 4. Draw Vessel (Schematic)
    # ---------------------------------------------------------
    # Simple boat shape definition
    l = 1.5 # length scale
    w = 0.6 # width scale
    
    # Boat vertices relative to center (pointing North/Up initially)
    boat_poly = np.array([
        [0, l],          # Bow
        [w, -l/2],       # Starboard Aft
        [w, -l],         # Starboard Transom
        [-w, -l],        # Port Transom
        [-w, -l/2],      # Port Aft
        [0, l]           # Close loop
    ])
    
    # Rotation Matrix for Heading (psi)
    # Note: In plotting (x=East, y=North), a compass heading psi (clockwise from North) 
    # corresponds to a standard rotation of (90 - psi). 
    # However, simpler to just rotate points manually: 
    # x' = x cos(psi) + y sin(psi) ... wait, NED convention:
    # North is X, East is Y. 
    # Let's stick to plot coordinates: X_plot = East, Y_plot = North.
    # A heading of 0 points Up (North). 
    
    R = np.array([
        [np.cos(psi), np.sin(psi)], 
        [-np.sin(psi), np.cos(psi)]
    ]) # This rotates "North" vector towards East for positive psi
    
    # Rotate and translate boat
    # We need to map Boat X (Surge) -> Plot Y (North)
    # Boat Y (Sway) -> Plot X (East)
    # To avoid confusion, let's define boat in plot coords directly:
    # Bow is at (0, l) (North).
    # Heading psi rotates this vector clockwise.
    
    def rotate_point(x, y, angle_rad):
        # Rotate clockwise by angle_rad
        x_new = x * np.cos(angle_rad) + y * np.sin(angle_rad)
        y_new = -x * np.sin(angle_rad) + y * np.cos(angle_rad)
        return x_new, y_new

    rotated_poly = []
    for bx, by in boat_poly:
        # Initial boat points forward along Y_plot. 
        # Standard rotation formula for clockwise from Y axis:
        px = bx * np.cos(psi) + by * np.sin(psi) 
        py = -bx * np.sin(psi) + by * np.cos(psi)
        # Wait, simple check: 
        # If psi=90 (East), Bow (0,1) should become (1,0).
        # x_new = 0 + 1*1 = 1. y_new = 0 + 0 = 0. Correct.
        rotated_poly.append([px + vessel_e, py + vessel_n])
        
    boat = patches.Polygon(rotated_poly, closed=True, facecolor='#d9d9d9', edgecolor='black', zorder=5)
    ax.add_patch(boat)
    
    # Draw COG / Vessel Center
    ax.plot(vessel_e, vessel_n, 'bo', markersize=4, zorder=6)
    
    # ---------------------------------------------------------
    # 5. Draw Polar State Variables
    # ---------------------------------------------------------
    
    # Position Vector (Rho)
    ax.annotate('', xy=(vessel_e, vessel_n), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color='black', lw=1.5, ls='-'))
    
    # Label Rho
    ax.text(vessel_e/2 - 0.5, vessel_n/2 + 0.5, r'$\rho = \sqrt{x^2 + y^2}$', fontsize=12)
    
    # Dashed projections
    ax.plot([vessel_e, vessel_e], [0, vessel_n], 'k--', alpha=0.5, lw=1) # Vertical to East axis
    ax.plot([0, vessel_e], [vessel_n, vessel_n], 'k--', alpha=0.5, lw=1) # Horizontal to North axis
    
    # Label x, y
    ax.text(vessel_e + 0.1, vessel_n/2, r'$y$ (East)', fontsize=10, color='gray')
    ax.text(vessel_e/2, vessel_n + 0.2, r'$x$ (North)', fontsize=10, color='gray')

    # ---------------------------------------------------------
    # 6. Angles (Gamma and Psi)
    # ---------------------------------------------------------
    
    # Gamma (Bearing / Phase angle) - Angle from North to Position Vector
    # Arc centered at origin, from North (90 deg in plot) to Vector angle
    # Vector angle in plot (math) is arctan(N/E). 
    # Let's calculate math degrees for the arc:
    vec_deg = np.degrees(np.arctan2(vessel_n, vessel_e))
    
    gamma_arc = patches.Arc((0, 0), 3, 3, angle=0, theta1=vec_deg, theta2=90, color='red', lw=1.5)
    ax.add_patch(gamma_arc)
    ax.text(0.8, 1.8, r'$\gamma$', fontsize=14, color='red')

    # Psi (Heading)
    # Draw a local North axis at the vessel for reference
    ax.plot([vessel_e, vessel_e], [vessel_n, vessel_n + 2.5], 'k-.', lw=0.8, alpha=0.7)
    
    # Arc for heading (from local North, clockwise to vessel heading)
    # Vessel heading in plot space is (90 - psi_deg)
    psi_arc = patches.Arc((vessel_e, vessel_n), 2.5, 2.5, angle=0, theta1=90-psi_deg, theta2=90, color='blue', lw=1.5)
    ax.add_patch(psi_arc)
    ax.text(vessel_e + 0.3, vessel_n + 1.5, r'$\psi$', fontsize=14, color='blue')

    # Add velocity vector (optional, but good for "motion control")
    # v_mag = 2.0
    # vx = v_mag * np.sin(psi)
    # vy = v_mag * np.cos(psi)
    # ax.arrow(vessel_e, vessel_n, vx, vy, head_width=0.15, fc='b', ec='b')
    # ax.text(vessel_e + vx, vessel_n + vy + 0.2, r'$u$', color='blue')

    # ---------------------------------------------------------
    # 7. Final Formatting
    # ---------------------------------------------------------
    ax.set_aspect('equal')
    ax.set_xlim(-1, 9)
    ax.set_ylim(-1, 9)
    ax.axis('off') # Turn off the actual plot box, stick to our drawn axes
    
    plt.title(r'Coordinate Transformation: NED to Polar State', fontsize=14, y=0.95)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    draw_marine_craft_coord_system()