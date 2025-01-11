import numpy as np

def simulate_cylinder_barrel_hits(d, L, Npoints=100, Ndir=10000, seed=42):
    """
    d    : Diameter of the cylinder
    L    : Length of the cylinder (z-axis from 0 to L)
    Npoints: Number of points along the diameter at z=0
    Ndir   : Number of photon directions to generate for each point
    seed   : Random seed for reproducibility
    
    Returns: A list (length Npoints) of the #barrel_hits for each emission point.
             Also returns the array of x-coordinates for those points.
    """
    np.random.seed(seed)
    
    R = d / 2.0  # radius
    # Create an array of x-coordinates for the Npoints on the diameter
    x_points = np.linspace(-R, R, Npoints)  # y=0, z=0

    # Prepare output
    barrel_hits_per_point = []

    for x0 in x_points:
        # We'll count how many of the Ndir photons from this point reach the barrel
        barrel_count = 0
        
        # Generate random directions in 4 pi
        # Method: generate cos(theta) ~ uniform(-1,1), phi ~ uniform(0,2pi)
        # Then convert to (vx, vy, vz)
        phi = 2.0 * np.pi * np.random.rand(Ndir)
        cos_theta = 2.0 * np.random.rand(Ndir) - 1.0  # in [-1, 1]
        theta = np.arccos(cos_theta)
        
        vx = np.sin(theta) * np.cos(phi)
        vy = np.sin(theta) * np.sin(phi)
        vz = np.cos(theta)  # or just cos_theta

        # For each direction, check intersection
        for i in range(Ndir):
            x_i0 = x0
            vx_i, vy_i, vz_i = vx[i], vy[i], vz[i]
            
            # Distance to barrel (solve (x0 + t vx)^2 + (t vy)^2 = R^2)
            # We have: (x0 + t*vx)^2 + (t*vy)^2 = R^2
            # Expand: x0^2 + 2 x0 t vx + t^2 (vx^2 + vy^2) = R^2
            # => a t^2 + b t + c = 0
            a = vx_i**2 + vy_i**2
            b = 2.0 * x_i0 * vx_i
            c = x_i0**2 - R**2
            
            t_side = None  # We'll keep the smallest positive root if it exists
            if abs(a) > 1e-15:
                disc = b**2 - 4.0*a*c
                if disc >= 0:
                    sqrt_disc = np.sqrt(disc)
                    t1 = (-b + sqrt_disc) / (2.0*a)
                    t2 = (-b - sqrt_disc) / (2.0*a)
                    # We want the smallest positive root
                    possibles = []
                    if t1 > 1e-15:
                        possibles.append(t1)
                    if t2 > 1e-15:
                        possibles.append(t2)
                    if len(possibles) > 0:
                        t_side = min(possibles)  # smallest positive root
                
            # If t_side is valid, we also need to check z in [0, L]
            # z(t) = 0 + t*vz_i
            hits_barrel = False
            if t_side is not None:
                z_side = t_side * vz_i
                if 0.0 <= z_side <= L:
                    hits_barrel = True
            
            # Distance to top endcap
            # we want z(t) = L => t = L/vz_i (must be positive)
            hits_top = False
            t_top = None
            if vz_i > 1e-15:
                t_top = L / vz_i
                # if t_top is > 0, it's a candidate
                # but we just check if it hits the top first
                # We'll see if it's before the side
                if t_side is None:
                    # No side intersection => definitely goes out top if t_top>0
                    hits_top = True
                else:
                    # Compare t_side vs t_top
                    if t_top < t_side:
                        # hits top first
                        hits_top = True
            
            if hits_barrel and not hits_top:
                # Means side intersection is first
                barrel_count += 1
            elif hits_barrel and hits_top:
                # If side intersection t_side < t_top => hits barrel
                if t_side < t_top:
                    barrel_count += 1
                # else => hits top first, no barrel
            else:
                # If we never had a valid side intersection or top intersection
                # it might go backward or out bottom, etc.
                # Or it might not hit side in [0, L].
                pass
        
        barrel_hits_per_point.append(barrel_count)
    
    return x_points, barrel_hits_per_point



def simulate_cylinder_barrel_sectors(d, L, Npoints=100, Ndir=10000, Nsectors=500, seed=42):
    """
    d        : Cylinder diameter
    L        : Cylinder length (z=0 to z=L)
    Npoints  : Number of source points along the diameter at z=0
    Ndir     : Number of photons emitted by each source point
    Nsectors : Number of equally spaced angular sectors on the barrel
    seed     : RNG seed for reproducibility
    
    Returns:
      x_coords               -> array of shape (Npoints,) giving the x-coord for each source point
      total_barrel_hits      -> array of shape (Npoints,) giving total # of photons that hit the barrel
      max_sector_barrel_hits -> array of shape (Npoints,) giving the maximum # of hits in any single sector
    """
    np.random.seed(seed)
    
    R = d / 2.0  # radius
    # The x-coordinates of the Npoints on the diameter
    x_coords = np.linspace(-R, R, Npoints)
    
    total_barrel_hits = np.zeros(Npoints, dtype=int)
    max_sector_barrel_hits = np.zeros(Npoints, dtype=int)
    
    for ipt, x0 in enumerate(x_coords):
        # We'll store which sector each barrel-hit goes to
        # Then we can do a histogram to see how many hits per sector
        sector_hits = np.zeros(Nsectors, dtype=int)
        
        # Generate random directions in 4 pi
        # phi uniform(0,2pi), cos(theta) uniform(-1,1)
        rand_phi = 2.0 * np.pi * np.random.rand(Ndir)
        cos_theta = 2.0 * np.random.rand(Ndir) - 1.0
        theta = np.arccos(cos_theta)  # in [0, pi]
        
        vx = np.sin(theta) * np.cos(rand_phi)
        vy = np.sin(theta) * np.sin(rand_phi)
        vz = cos_theta  # same as np.cos(theta)
        
        # Count how many photons from this point hit the barrel
        barrel_count = 0
        
        for i in range(Ndir):
            # Quadratic for side intersection
            #  (x0 + t vx)^2 + (t vy)^2 = R^2
            a = vx[i]**2 + vy[i]**2
            b = 2.0 * x0 * vx[i]
            c = x0**2 - R**2
            
            t_side = None
            if abs(a) > 1e-14:
                disc = b**2 - 4*a*c
                if disc >= 0:
                    sqrt_disc = np.sqrt(disc)
                    t1 = (-b + sqrt_disc) / (2*a)
                    t2 = (-b - sqrt_disc) / (2*a)
                    candidates = []
                    if t1 > 1e-14:
                        candidates.append(t1)
                    if t2 > 1e-14:
                        candidates.append(t2)
                    if len(candidates) > 0:
                        t_side = min(candidates)  # smallest positive root
            
            hits_barrel = False
            if t_side is not None:
                # Check z(t_side) in [0, L]
                z_side = t_side * vz[i]  # z0=0
                if 0.0 <= z_side <= L:
                    hits_barrel = True
            
            # Check top endcap (z=L), if vz>0
            hits_top = False
            t_top = None
            if vz[i] > 1e-14:
                t_top = L / vz[i]
                # Compare t_side vs t_top
                if t_side is None:
                    # No valid side intersection => goes out top
                    hits_top = True
                else:
                    # If top is reached before side => hits top
                    if t_top < t_side:
                        hits_top = True
            
            if hits_barrel and not hits_top:
                # Definitely side is first
                barrel_count += 1
                # Find the (x_side, y_side)
                x_side = x0 + t_side * vx[i]
                y_side = 0 + t_side * vy[i]  # y0=0
                # Convert to angle in [0, 2pi)
                phi_side = np.arctan2(y_side, x_side)
                if phi_side < 0:
                    phi_side += 2.0*np.pi
                # Determine sector index
                sector_idx = int((phi_side / (2.0*np.pi)) * Nsectors)
                # Make sure it's within [0, Nsectors-1]
                if sector_idx == Nsectors:
                    sector_idx = Nsectors - 1
                sector_hits[sector_idx] += 1
            
            elif hits_barrel and hits_top:
                # side, top both possible => check which is first
                if t_side < t_top:
                    barrel_count += 1
                    x_side = x0 + t_side * vx[i]
                    y_side = t_side * vy[i]
                    phi_side = np.arctan2(y_side, x_side)
                    if phi_side < 0:
                        phi_side += 2.0*np.pi
                    sector_idx = int((phi_side / (2.0*np.pi)) * Nsectors)
                    if sector_idx == Nsectors:
                        sector_idx = Nsectors - 1
                    sector_hits[sector_idx] += 1
                # else => hits top first => no barrel
            
            # if no barrel or top => presumably it goes out the bottom (vz<0),
            # or the side intersection wasn't valid in [0,L]. We don't count it.
        
        # Store the results for this source point
        total_barrel_hits[ipt] = barrel_count
        max_sector_hits = sector_hits.max() if barrel_count > 0 else 0
        max_sector_barrel_hits[ipt] = max_sector_hits
    
    return x_coords, total_barrel_hits, max_sector_barrel_hits



def simulate_photons_in_cylinder(R, Z, Npoints, Ngamma, reemit_prob, seed=123):
    """
    Simulate photons generated uniformly in a cylinder of radius R, length Z.
    Each of the Npoints points emits Ngamma photons isotropically.
    Photons stop when hitting the barrel (in which case we count them) or
    are absorbed after repeatedly hitting the end-caps. The end-caps
    re-emit with probability = reemit_prob, else absorb the photon.

    Returns:
      barrel_hit_count : total number of photons that hit the barrel.
      total_photons    : total photons generated = Npoints * Ngamma.
    """
    np.random.seed(seed)
    
    # 1) Generate Npoints uniformly in cylinder volume
    # r = R * sqrt(U1), theta = 2 pi U2, z = Z * U3
    U1 = np.random.rand(Npoints)
    U2 = np.random.rand(Npoints)
    U3 = np.random.rand(Npoints)
    
    r = R * np.sqrt(U1)
    theta = 2.0 * np.pi * U2
    z_cyl = Z * U3
    
    # Convert to cartesian
    x_cyl = r * np.cos(theta)
    y_cyl = r * np.sin(theta)
    # z_cyl is already in [0, Z]
    
    # We'll track how many photons eventually hit the barrel
    barrel_hit_count = 0
    total_photons = Npoints * Ngamma
    
    # Helper function to generate isotropic direction vector
    def random_direction():
        phi = 2.0 * np.pi * np.random.rand()
        cos_t = 2.0 * np.random.rand() - 1.0
        sin_t = np.sqrt(1.0 - cos_t**2)
        vx = sin_t * np.cos(phi)
        vy = sin_t * np.sin(phi)
        vz = cos_t
        return vx, vy, vz
    
    for i in range(Npoints):
        x0 = x_cyl[i]
        y0 = y_cyl[i]
        z0 = z_cyl[i]
        
        for _ in range(Ngamma):
            # Start a photon at (x0, y0, z0)
            # Generate initial direction
            vx, vy, vz = random_direction()
            
            alive = True
            px, py, pz = x0, y0, z0  # photon position
            
            while alive:
                # Compute intersection with barrel: (px + t vx)^2 + (py + t vy)^2 = R^2
                # Quadratic: a t^2 + b t + c = 0
                a = vx**2 + vy**2
                b = 2.0*(px*vx + py*vy)
                c = px**2 + py**2 - R**2
                
                t_barrel = None
                if abs(a) > 1e-15:
                    disc = b*b - 4.0*a*c
                    if disc >= 0:
                        sqrt_disc = np.sqrt(disc)
                        t1 = (-b + sqrt_disc) / (2.0*a)
                        t2 = (-b - sqrt_disc) / (2.0*a)
                        t_candidates = []
                        if t1 > 1e-15:
                            t_candidates.append(t1)
                        if t2 > 1e-15:
                            t_candidates.append(t2)
                        if len(t_candidates) > 0:
                            t_barrel = min(t_candidates)
                
                # Intersection with top end-cap z=Z if vz>0 => t_top = (Z - pz)/vz
                # Intersection with bottom end-cap z=0 if vz<0 => t_bottom = -pz/vz
                t_top = None
                if vz > 1e-15:
                    dt = (Z - pz)/vz
                    if dt > 1e-15:
                        t_top = dt
                t_bottom = None
                if vz < -1e-15:
                    dt = -pz/vz  # pz>0
                    if dt > 1e-15:
                        t_bottom = dt
                
                # Decide which intersection is smallest positive
                valid_times = []
                labels = []  # 'barrel', 'top', or 'bottom'
                
                if t_barrel is not None:
                    valid_times.append(t_barrel)
                    labels.append('barrel')
                if t_top is not None:
                    valid_times.append(t_top)
                    labels.append('top')
                if t_bottom is not None:
                    valid_times.append(t_bottom)
                    labels.append('bottom')
                
                if len(valid_times) == 0:
                    # No intersections forward => photon escapes outside or backward
                    # We end the photon (doesn't hit barrel)
                    alive = False
                    break
                
                # Find the smallest positive intersection
                idx_min = np.argmin(valid_times)
                t_min = valid_times[idx_min]
                surf_label = labels[idx_min]
                
                # Move photon to intersection
                px_new = px + t_min*vx
                py_new = py + t_min*vy
                pz_new = pz + t_min*vz
                
                if surf_label == 'barrel':
                    # Photon hits the barrel -> we count it, then it's done
                    barrel_hit_count += 1
                    alive = False
                else:
                    # It hits an end-cap (top or bottom)
                    # Absorb with prob (1 - reemit_prob), re-emit with prob reemit_prob
                    if np.random.rand() < reemit_prob:
                        # Re-emit from (px_new, py_new, pz_new)
                        px, py, pz = px_new, py_new, pz_new
                        vx, vy, vz = random_direction()
                        # continue the while loop
                    else:
                        # absorbed => photon terminates
                        alive = False
    
    return barrel_hit_count, total_photons


